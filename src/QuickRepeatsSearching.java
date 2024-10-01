
import java.util.Arrays;
import java.awt.BasicStroke;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.imageio.ImageIO;

public final class QuickRepeatsSearching {

    public void SetSequences(String[] seq, String[] sname) {
        this.seq = seq;
        this.sname = sname;
        nseq = seq.length;
    }

    public void SetFileName(String a) {
        filePath = a;
    }

    public void SetRepeatLen(int kmerln, int minlenseq, int gap) {
        this.kmerln = kmerln;
        if (this.kmerln < 12) {
            this.kmerln = 12;
        }
        this.gap = gap;
        if (this.gap < kmerln) {
            this.gap = kmerln;
        }
        this.minlenseq = minlenseq;
        if (this.minlenseq < kmerln) {
            this.minlenseq = kmerln;
        }
    }

    public void SetFlanks(int i) {
        flanks = i;
    }

    public void SetMaskedPicture(boolean i) {
        MaskedPicture = i;
    }

    public void SetShowSeq(boolean i) {
        SeqShow = i;
    }

    public void SetMasked(boolean i) {
        MaskedShow = i;
    }

    public void SetGFF(boolean i) {
        GFFShow = i;
    }

    public void SetImage(int w, int h) {
        if (w > 0) {
            iwidth = w;
        }
        if (h > 0) {
            iheight = h;
        }
    }

    public void Run() throws IOException {
        for (int i = 0; i < nseq; i++) {
            int l = seq[i].length();
            repeatslen = 0;
            reallen = 0;
            bb = new ArrayList<>();

            if (l > minlenseq) {
                startTime = System.nanoTime();

                LowComplexitySequence m1 = new LowComplexitySequence();
                m1.FindAllSSRs(seq[i], telomers);
                byte[] ssrmsk = m1.MapBytes();
                int[] ssr = m1.IntBlocks();

                bb.add(ssr);

                MaskingSequence ms = new MaskingSequence();

                int[] u = ms.Mask(seq[i], ssrmsk, kmerln, minlenseq);
                repeatslen = ms.RepeatLength();
                reallen = ms.ReallSeqLength();
                repeatslen = (repeatslen * 100) / (l - reallen);
                gaps = (reallen * 100) / l;

                System.out.println("Target sequence length = " + l + " nt");
                System.out.println("Sequence coverage by repeats =" + String.format("%.2f", repeatslen) + "%");
                System.out.println("Sequence gap (bp)=" + (int) reallen + " (" + String.format("%.4f", gaps) + "%)\n");

                if (u.length > 1) {
                    if (MaskedShow) {
                        MaskSave(i, u);
                    }
                    if (MaskedPicture) {
                        PictureMasking(u);
                    } else {
                        System.out.println("Clustering started...");
                        ClusteringMasking(seq[i], u, sim);
                    }
                }

                if (bb != null) {
                    if (GFFShow) {
                        GffSave(i);
                    }
                    PictureSave(i, iwidth, iheight);
                }

            }
        }
    }

    private int ClusteringMasking(String seq, int[] z2, int sim) {

        SequencesClustering sc = new SequencesClustering(seq, z2, sim);
        int[][] d = sc.ResultArray(); // d[j][0] = x1; d[j][1] = length;
        int[] q = sc.Result();        // cluster ID for each block
        int ncl = sc.getNcl();

        if (ncl < 1) {
            return -1;
        }

        for (int j = 1; j <= ncl; j++) {
            ArrayList<Integer> z = new ArrayList<>();
            z.add(0);
            for (int i = 0; i < q.length; i++) {
                if (q[i] == j) {
                    z.add(d[i][0]);
                    z.add(d[i][1]);
                }
            }
            bb.add(z.stream().mapToInt(Integer::intValue).toArray());
        }
        return bb.size();
    }

    private int PictureMasking(int[] x1) {
        int n = x1.length / 2;
        int[][] d = new int[n][2];
        int k = (n / 50 < 10) ? n : n / 50;

        for (int j = 0; j < n; j++) {
            int p = j * 2;
            d[j][0] = x1[p];
            d[j][1] = x1[p + 1] - x1[p];
        }
        Arrays.sort(d, (int[] a, int[] b) -> {
            return Integer.compare(b[1], a[1]);
        });

        for (int j = 0; j < n; j += k) {
            int[] k7 = new int[k + k + 1];
            k7[0] = k + k;
            int u = 0;
            int i = 0;
            int t = j + k;
            if (t > n) {
                t = n;
            }
            for (i = j; i < t; i++) {
                k7[++u] = d[i][0];
                k7[++u] = d[i][1];
            }
            bb.add(k7);
        }
        return bb.size();
    }

    public int[] ArrayTrim(int[] srcArray, int n) {
        int[] destArray = new int[n];
        System.arraycopy(srcArray, 0, destArray, 0, n);
        return destArray;
    }

    private void MaskSave(int n, int[] m) throws IOException {
        String maskedfile = filePath + "_" + (n + 1) + ".msk";
        if (nseq == 1) {
            maskedfile = filePath + ".msk";
        }
        try (FileWriter fileWriter = new FileWriter(maskedfile)) {
            System.out.println("Saving masked file: " + maskedfile);

            byte[] c = seq[n].getBytes();
            for (int j = 0; j < m.length; j += 2) {
                for (int i = m[j]; i <= m[j + 1]; i++) {
                    c[i] = (byte) (c[i] - 32);
                }
            }
            fileWriter.write(">" + sname[n] + " Sequence coverage by repeats = " + String.format("%.2f", repeatslen) + "%\n");
            fileWriter.write(new String(c));
        }
    }

    private void GffSave(int n) throws IOException {
        int k = 0;
        int l = seq[n].length();
        long duration = (System.nanoTime() - startTime) / 1000000000;
        System.out.println("Time taken: " + duration + " seconds\n");

        String reportfile = filePath + "_" + (n + 1) + ".gff";
        if (nseq == 1) {
            reportfile = filePath + ".gff";
        }

        try (FileWriter fileWriter = new FileWriter(reportfile); BufferedWriter bufferedWriter = new BufferedWriter(fileWriter)) {
            System.out.println("Saving report file: " + reportfile);
            StringBuilder sr = new StringBuilder();
            sr.append("kmer=").append(kmerln).append("\n").append("\n").append("Minimal repeat block size=").append(minlenseq).append("\n\n");
            sr.append("Sequence length = ").append(l).append("\n");
            sr.append("Sequence coverage by repeats = ").append(String.format("%.2f", repeatslen)).append("%\n");
            sr.append("Sequence gap (bp)=").append((int) reallen).append(" (").append(String.format("%.4f", gaps)).append("%)\n");
            sr.append("Time taken: ").append(duration).append(" seconds\n\n");
            sr.append("__________________________________________________\n Repeats search for: ").append(filePath).append("//").append(sname[n]).append(" ").append(l).append("bp :\n");
            if (SeqShow) {
                sr.append("Seqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\tSequence").append("\n");
            } else {
                sr.append("Seqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase").append("\n");
            }

            bufferedWriter.write(sr.toString());
            for (int i = 0; i < bb.size(); i++) {
                int[] z7 = bb.get(i);
                k++;
                for (int j = 1; j < z7.length; j += 2) {
                    String s0 = "";
                    int x = z7[j] + Math.abs(z7[j + 1]) - 1;

                    if (SeqShow) {
                        if (x > l) {
                            s0 = seq[n].substring(z7[j]);
                        } else {
                            s0 = seq[n].substring(z7[j], x);
                        }
                        if (flanks > 0) {
                            String s1 = "";
                            String s2 = "";
                            if (z7[j] - flanks > 0) {
                                s1 = seq[n].substring(z7[j] - flanks, z7[j]).toUpperCase();
                            } else {
                                if (z7[j] > 1) {
                                    s1 = seq[n].substring(1, z7[1] - 1).toUpperCase();
                                }
                            }
                            if (x + flanks < l) {
                                s2 = seq[n].substring(x, x + flanks).toUpperCase();
                            } else {
                                if (l - x > 0) {
                                    s2 = seq[n].substring(x, l).toUpperCase();
                                }
                            }
                            s0 = s1 + s0 + s2;
                        }
                        if (z7[j + 1] < 0) {
                            s0 = dna.ComplementDNA2(s0);
                        }
                    }
                    sr = new StringBuilder();
                    if (k == 1) {
                        sr.append(sname[n]).append("\t").append("SSR").append("\t").append(k).append("\t").append(z7[j] + 1).append("\t").append(x + 1).append("\t").append(z7[j + 1]).append("\t").append(".").append("\t").append(".").append("\t").append(s0).append("\n");
                    } else {
                        sr.append(sname[n]).append("\t").append(".").append("\t").append(k).append("\t").append(z7[j] + 1).append("\t").append(x + 1).append("\t").append(z7[j + 1]).append("\t").append(".").append("\t").append(".").append("\t").append(s0).append("\n");
                    }
                    bufferedWriter.write(sr.toString());
                }
            }
        }
    }

    private void PictureSave(int n, int dw, int dh) throws IOException {
        int maxClusters = 250;
        int maxImageDimension = 40_000;
        int minImageWidth = 1000;
        int minImageHeight = 100;
        int stepPadding = 50;
        int defaultStep = 20;

        // Adjust number of clusters `b`
        int b = Math.min(bb.size(), maxClusters); // Maximum of 1000 clusters

        // Adjust `z` (step between clusters) based on `b`
        int z = calculateClusterStep(b, defaultStep);

        // Calculate width and height
        int l = seq[n].length();
        int width = calculateWidth(l, dw, maxImageDimension, minImageWidth);
        int height = calculateHeight(b, z, dh, maxImageDimension, minImageHeight, stepPadding);

        // Calculate dot size
        float dotSize = calculateDotSize(b);

        // Attempt to save the image
        try {
            saveImage(n, l, b, z, width, height, dotSize);
        } catch (IOException e) {
            System.out.println("Error saving the picture: " + e.getMessage());
        }
    }

    private int calculateClusterStep(int b, int defaultStep) {
        if (b > 400) {
            return 7;
        }
        if (b > 250) {
            return 8;
        }
        if (b > 100) {
            return 9;
        }
        if (b > 50) {
            return 10;
        }
        return defaultStep;
    }

    private int calculateWidth(int l, int dw, int maxImageDimension, int minImageWidth) {
        int width = (l < 80_000_000) ? l / 200 : maxImageDimension;
        if (dw > 0) {
            width = dw;
        }
        return Math.max(Math.min(width, maxImageDimension), minImageWidth);
    }

    private int calculateHeight(int b, int z, int dh, int maxImageDimension, int minImageHeight, int stepPadding) {
        int height = b * z + stepPadding; // Adding some padding
        if (dh > 0) {
            height = dh;
        }
        return Math.max(Math.min(height, maxImageDimension), minImageHeight);
    }

    private float calculateDotSize(int b) {
        float dotSize = 10 - (b / 100.0f);
        return Math.max(5.0f, Math.min(dotSize, 7.0f));
    }

    private void saveImage(int n, int l, int b, int z, int width, int height, float dotSize) throws IOException {
        double nucleotidesPerPixel = (double) width / l;
        String pngfile = filePath + "_" + (n + 1) + ".png";
        if (nseq == 1) {
            pngfile = filePath + ".png";
        }

        System.out.println("Saving picture " + (width + 100) + "x" + (height + 100) + " : " + pngfile);

        BufferedImage image = new BufferedImage(width + 100, height + 100, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();

        g2d.setStroke(new BasicStroke(dotSize));
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, width + 100, height + 100);
        g2d.setColor(Color.BLACK);
        drawLinesAndLabels(g2d, l, width);
        drawClusters(g2d, b, z, nucleotidesPerPixel);

        g2d.dispose();
        File outputFile = new File(pngfile);
        ImageIO.write(image, "png", outputFile);
    }

    private void drawLinesAndLabels(Graphics2D g2d, int l, int width) {
        g2d.drawLine(50, 15, width + 50, 15); // top line

        int w = (width + 100) / 10;
        int d = l / 10;
        for (int i = 0; i < 12; i++) {
            g2d.drawLine(i * w + 50, 5, i * w + 50, 15);
            g2d.drawString(String.format("%,d", (1 + i * d)), 55 + i * w, 10);
        }
    }

    private void drawClusters(Graphics2D g2d, int b, int z, double w1) {
        for (int i = 0; i < b; i++) {
            int[] z7 = bb.get(i);

            // Gray lines at height 22
            for (int j = 1; j < z7.length; j += 2) {
                int x1 = 50 + (int) (z7[j] * w1);
                int x2 = 50 + (int) ((z7[j] + z7[j + 1]) * w1);
                g2d.setColor(Color.DARK_GRAY);
                g2d.drawLine(x1, 22, x2, 22); // draw dark gray line
            }

            // Blue lines at each cluster step
            int y = 50 + (i * z);
            for (int j = 1; j < z7.length; j += 2) {
                int x1 = 50 + (int) (z7[j] * w1);
                int x2 = 50 + (int) ((z7[j] + z7[j + 1]) * w1);
                g2d.setColor(Color.BLUE);
                g2d.drawLine(x1, y, x2, y); // draw blue line
            }
        }
    }

    private double gaps = 0;
    private double reallen = 0;
    private double repeatslen = 0;
    private long startTime;
    private int nseq = 0;
    private int iwidth = 0;
    private int iheight = 0;
    private int minlenseq = 50;     // Minimal repeat block size
    private int kmerln = 21;       // kmer=12-21
    private int flanks = 20;
    private int gap = 21;         // gap between repeat blocks, gap=kmer
    private final int sim = 50;
    private final int telomers = 14; // Kmax=11 ->SSR  Kmax=14 -> telomers
    private boolean SeqShow;
    private boolean MaskedShow;
    private boolean MaskedPicture;
    private boolean GFFShow;
    private String filePath;
    private String[] seq;
    private String[] sname;
    private ArrayList<int[]> bb;
}
