import java.awt.BasicStroke;
import java.io.FileWriter;
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

    public void SetKmerLen(int i) {
        kmerln = i;
        if (kmerln < 12) {
            kmerln = 12;
        }
        minlenblock = kmerln + kmerln;
    }

    public void SetMinimalRepeatLen(int i) {
        minlenblock = i;
    }

    public void SetSequenceLen(int i) {
        minlenseq = i;
        if (minlenseq < minlenblock) {
            minlenseq = minlenblock;
        }
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

    public void SetFlanks(int i) {
        flanks = i;
    }

    public void SetShowSeq(boolean i) {
        SeqShow = i;
    }

    public void Run() throws IOException {
        startTime = System.nanoTime();
        for (int i = 0; i < nseq; i++) {
            reallen = 0;
            byte[] u = Mask(seq[i], kmerln);
            MaskSave(i, u);
            ClusteringMasking(seq[i], u, kmerln, minlenblock, minlenseq);
            if (GFFShow) {
                GffSave(i);
            }
            PictureSave(i, iwidth, iheight);
        }
    }

    private byte[] Mask(String seq, int kmer) {
        int l = seq.length();
        int i = 0;
        int j = 0;
        int p = 0;

        String s;
        String aseq = dna.ComplementDNA(seq);
        HashMap<String, Integer> px2 = new HashMap<>();

        int[] ax = new int[kmer];
        int[] bx = new int[5];

        byte b[] = seq.getBytes();
        byte c[] = aseq.getBytes();
        for (i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];
            c[i] = tables.dx2[c[i]];
            if (bx[4] == 0) {
                reallen++;
            }
        }

        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                s = seq.substring(i - kmer + 1, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s) + 1;
                    px2.replace(s, p);
                } else {
                    px2.put(s, 1);
                }
            }
            bx[b[i + 1 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }
        //reverse
        bx = new int[5];
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = c[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = c[i];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                s = aseq.substring(i - kmer + 1, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s) + 1;
                    px2.replace(s, p);
                }
            }
            bx[c[i + 1 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }

        byte[] u = new byte[l];
        bx = new int[5];
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            int x = i - kmer + 1;
            if (bx[4] == 0) {
                s = seq.substring(x, i + 1);
                p = px2.get(s);
                if (p > 1) {
                    for (j = x; j < x + kmer; j++) {
                        u[j] = 1;
                    }
                }
            }
            bx[b[i + 1 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }
        //reverse
        bx = new int[5];
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = c[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = c[i];
            if (bx[4] == 0) {
                int x = i - kmer + 1;
                s = aseq.substring(x, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s);
                    if (p > 1) {
                        for (j = x; j < x + kmer; j++) {
                            u[j] = 2;
                        }
                    }
                }
            }
            bx[c[i + 1 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }
        return u;
    }

    private int ClusteringMasking(String seq, byte[] u, int kmer, int minlenblock, int minlenseq) {
        ArrayList<Integer> z = new ArrayList<>();
        int j = 0;
        for (int i = 0; i < u.length; i++) {
            if (u[i] > 0) {
                int x1 = i;
                int x2 = i;
                for (j = i + 1; j < u.length; j++) {
                    if (u[j] > 0) {
                        x2 = j;
                    } else {
                        break;
                    }
                }
                i = j - 1;
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
        int[] z2 = new int[n];//x1-x2 block
        n = -1;
        for (int i = 0; i < z.size(); i += 2) {
            int x1 = z.get(i);
            int x2 = z.get(i + 1);
            if (x2 - x1 > minlenseq) { // join blocks at short distance
                z2[++n] = x1;
                z2[++n] = x2;
            }
        }

        SequencesClustering sc = new SequencesClustering(seq, z2, 70, true, kmer);
        int[] q = sc.Result(); // cluster ID for each block
        int ncl = sc.getNcl();
        if (q.length < 1) {
            return -1;
        }

        bb = new ArrayList<>();
        for (int f = 1; f < ncl + 1; f++) {
            int[] k7 = new int[q.length + q.length + 1];
            int h = 0;
            int p = 0;
            for (j = 0; j < q.length; j++) {
                p = j * 2;
                if (q[j] == f) {
                    k7[0] = k7[0] + 2;
                    h++;
                    k7[h] = z2[p];
                    h++;
                    k7[h] = z2[p + 1] - z2[p] + 1;
                }
            }
            if (h > 0) {
                k7 = ArrayTrim(k7, k7[0] + 1);
                for (int i = 2; i < k7.length; i += 2) {
                    if (k7[0] < k7[i]) {
                        k7[0] = k7[i];
                    }
                }
                bb.add(k7);
            }
        }

        SortBlocks();
        return bb.size();
    }

    private void SortBlocks() {
        if (bb == null) {
            return;
        }
        int n = bb.size();
        for (int i = 0; i < n - 1; i++) {
            int minIndex = i;
            int[] zi = bb.get(i);
            for (int j = i + 1; j < n; j++) {
                int[] zj = bb.get(j);
                if (zj[0] > zi[0]) {
                    minIndex = j;
                    zi = zj; // Update zi to the new minimum element
                }
            }
            // Swap the elements at i and minIndex
            int[] temp = bb.get(i);
            bb.set(i, bb.get(minIndex));
            bb.set(minIndex, temp);
        }
    }

    public int[] ArrayTrim(int[] srcArray, int n) {
        int[] destArray = new int[n];
        System.arraycopy(srcArray, 0, destArray, 0, n);
        return destArray;
    }

    private void MaskSave(int n, byte[] m) throws IOException {
        if (m == null) {
            return;
        }
        int l = seq[n].length();
        repeatslen = 0;

        String maskedfile = filePath + "_" + (n + 1) + ".msk";
        if (nseq == 1) {
            maskedfile = filePath + ".msk";
        }

        try (FileWriter fileWriter = new FileWriter(maskedfile)) {
            fileWriter.write(">" + sname[n]);
            System.out.println("Saving masked file: " + maskedfile);
            char[] c = seq[n].toCharArray();
            for (int i = 0; i < l; i++) {
                if (m[i] > 0) {              // if (m[i] == 0) {
                    repeatslen++;
                    c[i] = (char) (c[i] - 32);
                }
            }

            double z = (repeatslen * 100) / reallen;
            System.out.println("Sequence coverage by repeats = " + String.format("%.2f", z) + "%");
            if (MaskedShow) {
                fileWriter.write("Sequence coverage by repeats = " + String.format("%.2f", z) + "%\n\n");
                fileWriter.write(new String(c));
            }
        }
    }

    private void GffSave(int n) throws IOException {
        if (bb == null) {
            return;
        }
        int k = 0;
        int l = seq[n].length();
        long duration = (System.nanoTime() - startTime) / 1000000000;

        String reportfile = filePath + "_" + (n + 1) + ".gff";
        if (nseq == 1) {
            reportfile = filePath + ".gff";
        }

        StringBuilder sr = new StringBuilder();
        sr.append("kmer=").append(kmerln).append("\n").append("Minimal repeat=").append(minlenblock).append("\n").append("\n");
        sr.append("__________________________________________________\n Repeats search for: ").append(filePath).append("//").append(sname[n]).append(" ").append(l).append("bp :\n");
        sr.append("Time taken: ").append(duration).append(" seconds\n\n");
        if (SeqShow) {
            sr.append("Seqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\tSequence").append("\n");
        } else {
            sr.append("Seqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\tSequence").append("\n");
        }

        double d = ((l - reallen) * 100) / l;
        double z = (repeatslen * 100) / reallen;
        sr.append("Sequence gap (bp)=").append(l - reallen).append(" (").append(String.format("%.3f", d)).append("%)\n");

        try (FileWriter fileWriter = new FileWriter(reportfile)) {
            System.out.println("Saving report file: " + reportfile);
            fileWriter.write("Sequence coverage by repeats = " + String.format("%.2f", z) + "%\n\n");

            fileWriter.write(sr.toString());
            for (int i = 0; i < bb.size(); i++) {
                int[] z7 = bb.get(i);
                k++;
                for (int j = 1; j < z7.length; j += 2) {

                    String s0 = "";
                    int x = z7[j] + Math.abs(z7[j + 1]) - 1;

                    if (SeqShow) {
                        if (x > seq[n].length()) {
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
                    fileWriter.write(sr.toString());
                }
            }
        }
    }

    private void PictureSave(int n, int dw, int dh) throws IOException {
        if (bb == null) {
            return;
        }
        int b = bb.size();
        int z = 20;          // step between clusters

        if (b > 1000) {
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
        int width = l > 5_000_000 ? 10000 : l / 500; // image=10000x3000 l < 5_000_000 ? l / 250 : 5_000 + (l - 5_000) / 250;

        if (dw > 0) {
            width = dw;
        }
        if (width > 46340) { //too big a picture leads to an error  
            width = 46340; //https://jobcardsystems.com/index.php/blog/29-effective-handling-of-big-images-in-java
        }
        if (width < 1000) {
            width = 1000;
        }
        int height = b * z;      //too big a picture leads to an error
        if (dh > 0) {
            height = dh;
        }
        if (height > 46340) {
            height = 46340;
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

        double w1 = (double) width / l;       // nucleotides per pixel        
        String pngfile = filePath + "_" + (n + 1) + ".png";
        if (nseq == 1) {
            pngfile = filePath + ".png";
        }

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
            int y = 40 + (i * z);  // y1-y2 line    
            for (int j = 1; j < z7.length; j += 2) {
                int x1 = 50 + (int) (z7[j] * w1);
                int x2 = 51 + (int) ((z7[j] + z7[j + 1]) * w1);
                g2d.setColor(Color.BLUE);
                g2d.drawLine(x1, y, x2, y); //(int x1, int y1, int x2, int y2)
            }
        }
        g2d.dispose();
        try {
            System.out.println("Saving picture : " + pngfile);
            File outputFile = new File(pngfile);
            ImageIO.write(image, "png", outputFile);
        } catch (IOException e) {
            System.out.println("Saving picture is large.");
        }
    }

    private String filePath;
    public int nseq;
    public int iwidth = 0;
    public int iheight = 0;
    public double reallen = 0;
    public double repeatslen = 0;
    private String[] seq;
    private String[] sname;
    private int minlenblock = 25;   // repeat length user control  
    private int minlenseq = 50;     // sequence length 
    private int kmerln = 21;
    private int flanks = 20;
    private boolean SeqShow;
    private ArrayList<int[]> bb;
    private long startTime;
    private boolean MaskedShow;
    private boolean GFFShow;
}
