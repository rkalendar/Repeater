
import java.util.Arrays;
import java.awt.BasicStroke;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.imageio.ImageIO;

public final class Tandems2 {

    public void SetSequences(String[] seq, String[] sname) {
        this.seq = seq;
        this.sname = sname;
        nseq = seq.length;
    }

    public void SetFileName(String a) {
        filePath = a;
    }

    public void SetRepeatLen(int kmer, int minlen, int minseq, int g) {
        kmerln = kmer;
        gap = g;
        minlenblock = minlen;
        minlenseq = minseq;
        if (kmerln < 12) {
            kmerln = 12;
        }
        if (minlenblock < kmer) {
            minlenblock = kmer;
        }
        if (minlenseq < minlenblock) {
            minlenseq = minlenblock;
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
            reallen = 0;
            repeatslen = 0;

            if (l > minlenseq) {
                startTime = System.nanoTime();
                int[] u = Mask(seq[i], kmerln, minlenblock, minlenseq);
                if (u.length > 1) {
                    repeatslen = (repeatslen * 100) / (l - reallen);
                    gaps = (reallen * 100) / l;

                    System.out.println("Target sequence length = " + l + " nt");
                    System.out.println("Sequence coverage by repeats =" + String.format("%.2f", repeatslen) + "%");
                    System.out.println("Sequence gap (bp)=" + (int) reallen + " (" + String.format("%.4f", gaps) + "%)\n");

                    if (MaskedShow) {
                        MaskSave(i, u);
                    }
                    bb = new ArrayList<>();

                    if (MaskedPicture) {
                        PictureMasking(u);
                    } else {
                        ClusteringMasking(seq[i], u, sim);
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

    private int[] Mask(String seq, int kmer, int minlenblock, int minlenseq) {
        HashMap<String, Integer> map = new HashMap<>();
        String aseq = dna.ComplementDNA(seq);

        int l = seq.length();
        int[] ax = new int[kmer];
        int[] bx = new int[5];

        byte b[] = seq.getBytes();
        byte c[] = aseq.getBytes();

        for (int i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];
            c[i] = tables.dx2[c[i]];
            if (b[i] == 4) {
                reallen++;
            }
        }

        for (int i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (int i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                String s = seq.substring(i - kmer + 1, i + 1);
                if (map.containsKey(s)) {
                    int p = map.get(s) + 1;
                    map.replace(s, p);
                } else {
                    map.put(s, 1);
                }
            }
            bx[b[i + 1 - kmer]]--;
            for (int j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }
        //reverse
        bx = new int[5];
        for (int i = 0; i < kmer - 1; i++) {
            ax[i] = c[i];
            bx[ax[i]]++;
        }
        for (int i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = c[i];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                String s = aseq.substring(i - kmer + 1, i + 1);
                if (map.containsKey(s)) {
                    int p = map.get(s) + 1;
                    map.replace(s, p);
                }
            }
            bx[c[i + 1 - kmer]]--;
            for (int j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }

        int[] u = new int[l];
        bx = new int[5];
        for (int i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (int i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            int x = i - kmer + 1;
            if (bx[4] == 0) {
                String s = seq.substring(x, i + 1);
                int p = map.get(s);
                if (p > 1) {
                    for (int j = x; j < x + kmer; j++) {
                        u[j]++;
                    }
                }
            }
            bx[b[i + 1 - kmer]]--;
            for (int j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }
        //reverse
        bx = new int[5];
        for (int i = 0; i < kmer - 1; i++) {
            ax[i] = c[i];
            bx[ax[i]]++;
        }
        for (int i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = c[i];
            if (bx[4] == 0) {
                int x = i - kmer + 1;
                String s = aseq.substring(x, i + 1);
                if (map.containsKey(s)) {
                    int p = map.get(s);
                    if (p > 1) {
                        for (int j = x; j < x + kmer; j++) {
                            u[j]++;
                        }
                    }
                }
            }
            bx[c[i + 1 - kmer]]--;
            for (int j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }

        // combining
        ArrayList<Integer> z = new ArrayList<>();
        for (int i = 0; i < u.length; i++) {
            if (u[i] > 0) {
                int x1 = i;
                int x2 = i;
                int y = 0;
                int x = i + 1;
                for (int j = x; j < u.length; j++) {
                    if (u[j] > 1) {
                        x2 = j;
                    } else {
                        y++;
                        if (y > gap) {
                            break;
                        }
                    }
                    i = j;
                }

                if (x2 > x1) {
                    if (z.isEmpty()) {
                        z.add(x1);
                        z.add(x2);
                    } else {
                        int h = z.size() - 1;
                        int x4 = z.get(h);
                        if (x4 + minlenblock > x1) {
                            z.set(h, x2);
                        } else {
                            z.add(x1);
                            z.add(x2);
                        }
                    }
                }
            }
        }

        int n = 0;
        for (int i = 0; i < z.size(); i += 2) {
            int x1 = z.get(i);
            int x2 = z.get(i + 1);
            if (x2 - x1 > minlenseq) {
                n += 2;
            }
        }
        u = new int[n];   //x1-x2 block
        n = -1;
        for (int i = 0; i < z.size(); i += 2) {
            int x1 = z.get(i);
            int x2 = z.get(i + 1);
            if (x2 - x1 > minlenseq) { // join blocks at short distance
                u[++n] = x1;
                u[++n] = x2;
                repeatslen = repeatslen + (x2 - x1 + 1);
            }
        }
        return u;
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
            sr.append("kmer=").append(kmerln).append("\n").append("Minimal repeat=").append(minlenblock).append("\n").append("Repeat filter=").append(minlenseq).append("\n\n");
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
                    sr.append(sname[n]).append("\t").append(".").append("\t").append(k).append("\t").append(z7[j] + 1).append("\t").append(x + 1).append("\t").append(z7[j + 1]).append("\t").append(".").append("\t").append(".").append("\t").append(s0).append("\n");
                    bufferedWriter.write(sr.toString());
                }
            }
        }
    }

    private void PictureSave(int n, int dw, int dh) throws IOException {
        int b = bb.size();
        int z = 20;      // step between clusters
        if (b > 1000) {  // maximum 1000 clusters 
            b = 1000;
        }
        if (b > 50) {
            z = 10;
        }
        if (b > 100) {
            z = 9;
        }
        if (b > 250) {
            z = 8;
        }
        if (b > 500) {
            z = 7;
        }

        int l = seq[n].length();
        int width = l < 80_000_000 ? l / 200 : 40_000;
        if (dw > 0) {
            width = dw;
        }
        if (width > 40_000) { //46_340 too big a picture leads to an error   2,147,483,647 pixels (or 46,340 x 46,340 pixels). 
            width = 40_000; //https://jobcardsystems.com/index.php/blog/29-effective-handling-of-big-images-in-java
        }
        if (width < 1000) {
            width = 1000;
        }
        int height = b * z;      //too big a picture leads to an error
        if (dh > 0) {
            height = dh;
        }
        if (height > 40_000) {
            height = 40_000;
        }
        if (height < 100) {
            height = 100;
        }
        height = height + 50;

        float dotSize = 10 - (b / 100);
        if (dotSize < 5) {
            dotSize = 5.0f;
        }
        if (dotSize > 7) {
            dotSize = 7.0f;
        }

        try {
            double w1 = (double) width / l;       // nucleotides per pixel        
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
            g2d.drawLine(50, 15, width + 50, 15);

            int w = (width + 100) / 10;
            int d = l / 10;
            for (int i = 0; i < 12; i++) {
                g2d.drawLine(i * w + 50, 5, i * w + 50, 15);
                g2d.drawString(String.format("%,d", (1 + (int) (i * d))), 55 + i * w, 10);
            }

            for (int i = 0; i < b; i++) {
                int[] z7 = bb.get(i);
                for (int j = 1; j < z7.length; j += 2) {
                    int x1 = 50 + (int) (z7[j] * w1);
                    int x2 = 50 + (int) ((z7[j] + z7[j + 1]) * w1);
                    g2d.setColor(Color.DARK_GRAY);
                    g2d.drawLine(x1, 22, x2, 22); //(int x1, int y1, int x2, int y2)
                }
            }

            for (int i = 0; i < b; i++) {
                int[] z7 = bb.get(i);
                int y = 50 + (i * z);  // y1-y2 line    
                for (int j = 1; j < z7.length; j += 2) {
                    int x1 = 50 + (int) (z7[j] * w1);
                    int x2 = 50 + (int) ((z7[j] + z7[j + 1]) * w1);
                    g2d.setColor(Color.BLUE);
                    g2d.drawLine(x1, y, x2, y); //(int x1, int y1, int x2, int y2)
                }
            }
            g2d.dispose();
            File outputFile = new File(pngfile);
            ImageIO.write(image, "png", outputFile);
        } catch (IOException e) {
            System.out.println("Saving picture is large.");
        }
    }

    private double reallen = 0;
    private double gaps = 0;
    private double repeatslen = 0;
    private long startTime;
    private int nseq = 0;
    private int iwidth = 0;
    private int iheight = 0;
    private int minlenblock = 30;   // repeat length user control  
    private int minlenseq = 50;     // sequence length 
    private int kmerln = 21;
    private int flanks = 20;
    private int gap = 21;
    private final int sim = 70;
    private boolean SeqShow;
    private boolean MaskedShow;
    private boolean MaskedPicture;
    private boolean GFFShow;
    private String filePath;
    private String[] seq;
    private String[] sname;
    private ArrayList<int[]> bb;
}
