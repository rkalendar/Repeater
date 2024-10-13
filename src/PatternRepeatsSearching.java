
import java.awt.BasicStroke;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.imageio.ImageIO;

public final class PatternRepeatsSearching {

    public void SetSequences(String[] seq, String[] sname) {
        this.seq = seq;
        this.sname = sname;
        nseq = seq.length;
    }

    public void SetFileName(String a) {
        filePath = a;
    }

    public void SetRepeatLen(int kmer, int minlen, int minseq) {
        kmerln = kmer;
        minlenblock = minlen;
        minlenseq = minseq;
        if (kmerln < 12) {
            kmerln = 12;
        }
        if (minlenblock < kmerln) {
            minlenblock = kmerln;
        }
        if (minlenseq < minlenblock) {
            minlenseq = minlenblock;
        }
    }

    public void SetImage(int w, int h) {
        if (w > 0) {
            iwidth = w;
        }
        if (h > 0) {
            iheight = h;
        }
    }

    public void SetNumberUnits(int i) {
        nblocks = i - 1;
        if (nblocks < 1) {
            nblocks = 1;
        }
        if (nblocks > 1000) {
            nblocks = 1000;
        }
    }

    public void SetTelomers(boolean i) {
        if (i) {
            telomers = 14;
        }
    }

    public void SetFlanks(int i) {
        flanks = i;
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

    public void Run() throws IOException {
        startTime = System.nanoTime();
        for (int i = 0; i < nseq; i++) {
            FindAllRepeats(seq[i], kmerln);
            PrintResult(i);
            PictureSave(i, iwidth, iheight);
        }
    }

    public void RunSSR() throws IOException {
        for (int i = 0; i < nseq; i++) {
            LowComplexitySequence m1 = new LowComplexitySequence();
            m1.FindAllSSRs(seq[i], telomers);
            bb = m1.Blocks();
            PrintResult(i);
            PictureSave(i, iwidth, iheight);
        }
    }

    public int[] ArrayExtend(int[] srcArray, int n) {
        int[] destArray = new int[srcArray.length + n];
        System.arraycopy(srcArray, 0, destArray, 0, srcArray.length);
        return destArray;
    }

    public int[] ArrayTrim(int[] srcArray, int n) {
        int[] destArray = new int[n];
        System.arraycopy(srcArray, 0, destArray, 0, n);
        return destArray;
    }

    public byte[] ArrayExtendByte(byte[] srcArray, int n) {
        byte[] destArray = new byte[srcArray.length + n];
        System.arraycopy(srcArray, 0, destArray, 0, srcArray.length);
        return destArray;
    }

    private int FindAllRepeats(String seq, int kmer) {
        int k, n, t, i, h, e, y, z, w, x, q, j, r, p;
        int l = seq.length();
        int g = l + l + 1;

        String s;
        String aseq = dna.ComplementDNA(seq);
        HashMap<String, int[]> px2 = new HashMap<>();

        int[] ax = new int[kmer];
        int[] bx = new int[5];

        byte b[] = seq.getBytes();
        for (i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];
        }
        b = ArrayExtendByte(b, l + 1);
        b[l] = 4;
        for (i = 1; i < l + 1; i++) {
            b[l + i] = tables.cdnat2[b[l - i]];
        }
        numnonn = kmer - 1;
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                numnonn++;
                s = seq.substring(i - kmer + 1, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s)[0] + 1;
                    px2.put(s, new int[]{p, 0});
                } else {
                    px2.put(s, new int[]{1, 0});
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
            ax[i] = b[l + i + 1];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[l + i + 1];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                s = aseq.substring(i - kmer + 1, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s)[0] + 1;
                    px2.put(s, new int[]{p, 0});
                }
            }
            bx[b[l + i + 2 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }

        t = 0;
        n = 0;
        bx = new int[5];
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                s = seq.substring(i - kmer + 1, i + 1);
                p = px2.get(s)[0];
                if (p > nblocks) {
                    t++;
                    n = n + p;
                    px2.get(s)[0] = -p;
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
            ax[i] = b[l + i + 1];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[l + i + 1];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                s = aseq.substring(i - kmer + 1, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s)[0];
                    if (p > nblocks) {
                        t++;
                        n = n + p;
                        px2.get(s)[0] = -p;
                    }
                }
            }
            bx[b[l + i + 2 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }
        if (n == 0 | t == 0) {
            return 0;
        }

        int[] u = new int[n + 1];
        int[][] x1 = new int[t + 1][2];
        z = 1;
        t = 0;
        bx = new int[5];
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                s = seq.substring(i - kmer + 1, i + 1);
                p = px2.get(s)[0];
                if (p < 0) {
                    px2.get(s)[0] = -px2.get(s)[0];
                    if (px2.get(s)[0] > nblocks) {
                        px2.get(s)[1] = z;
                        u[z] = i + 1 - kmer;
                        t++;
                        u[0] = t;
                        x1[t][1] = z;
                        x1[t][0] = px2.get(s)[0];
                        z = z + px2.get(s)[0];
                    }
                } else {
                    if (p > nblocks) {
                        px2.get(s)[1]++;
                        u[px2.get(s)[1]] = i + 1 - kmer;
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
            ax[i] = b[l + i + 1];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[l + i + 1];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                s = aseq.substring(i - kmer + 1, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s)[0];
                    if (p < 0) {
                        px2.get(s)[0] = -px2.get(s)[0];
                        if (px2.get(s)[0] > nblocks) {
                            px2.get(s)[1] = z;
                            u[z] = l + i + 2 - kmer;
                            t++;
                            u[0] = t;
                            x1[t][1] = z;
                            x1[t][0] = px2.get(s)[0];
                            z = z + px2.get(s)[0];
                        }
                    } else {
                        if (p > nblocks) {
                            px2.get(s)[1]++;
                            u[px2.get(s)[1]] = l + i + 2 - kmer;
                        }
                    }
                }
            }
            bx[b[l + i + 2 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }

        if (t == 0) {
            return 0;
        }

        bb = new ArrayList<>();
        byte[] m = new byte[g];
        int[] k7;
        int[] k8;

        for (i = 1; i < t + 1; i++) {
            k7 = new int[x1[i][0] + 1];

            for (k = 1; k <= x1[i][0]; k++) {
                x = u[x1[i][1] + k - 1];
                if (m[x] == 0) {
                    for (r = k + 1; r <= x1[i][0]; r++) {
                        y = u[x1[i][1] + r - 1];
                        h = kmer;
                        e = 0;// no mistmatches
                        p = 2; // fully complement 2 bases

                        for (;;) {
                            if (y > l) {
                                if (y + h > g - 1) {
                                    break;
                                }
                            } else {
                                if (y + h > l - 1) {
                                    break;
                                }
                            }
                            if (x < l) {
                                if (x + h > l - 1) {
                                    break;
                                }
                            } else {
                                if (x + h > g - 1) {
                                    break;
                                }
                            }
                            if (b[x + h] == 4 && b[y + h] == 4) {
                                break;
                            }

                            if (b[x + h] == b[y + h]) {
                                if (p > 0) {
                                    e = 0;
                                }
                                p++;
                            } else {
                                if (e > 2) {
                                    h = h - 3;
                                    break;
                                }
                                e++;
                                p = 0;
                            }
                            h++;
                        }
                        if (h > minlenblock) {
                            k7[0] = 1;
                            if (k7[k] < h) {
                                k7[k] = h;
                            }
                            if (k7[r] < h) {
                                k7[r] = h;
                            }

                        }

                    }
                }
            }

            if (k7[0] == 1) {
                n++;
                w = 0;
                k8 = new int[1 + x1[i][0] + x1[i][0]];
                for (k = 1; k <= x1[i][0]; k++) {

                    if (k7[k] > 0) {
                        k7[0] = k7[0] + k7[k];
                        x = u[x1[i][1] + k - 1];

                        for (r = 0; r < k7[k]; r++) {
                            m[x + r] = 1;
                        }
                        if (k8[0] > 0) {
                            if (k8[k8[0] - 1] < x && k8[k8[0] - 1] + k8[k8[0]] > x) {
                                if (k8[k8[0]] < k7[k]) {
                                    k8[k8[0]] = x + k7[k] - 2 - k8[k8[0] - 1];
                                }
                                x = -1;
                            }
                        }

                        if (x > -1) {
                            k8[0] = k8[0] + 2;
                            w++;
                            k8[w] = x;
                            w++;
                            k8[w] = k7[k];
                        }

                    }
                }

                k8 = ArrayTrim(k8, w + 1);
                bb.add(k8);
            }
        }
        n = bb.size();
        if (n == 0) {
            return 0;
        }

        for (i = 0; i < n; i++) {
            k7 = bb.get(i);
            for (r = 1; r < k7[0]; r += 2) {
                if (k7[r] >= l) {
                    k7[r] = l - (k7[r] - l) + 2 - k7[r + 1];
                    k7[r + 1] = -k7[r + 1];
                }
                w = k7[r] + Math.abs(k7[r + 1]);
                if (w > l) {
                    q = w - l;
                    if (k7[r + 1] > 0) {
                        k7[r + 1] = k7[r + 1] - q;
                    } else {
                        k7[r + 1] = k7[r + 1] + q;
                    }
                }
            }
            bb.set(i, Sorting(k7));
        }

        p = -1;
        for (i = 0; i < n; i++) {
            k7 = bb.get(i);
            if (k7[0] > 0) {
                p++;
                bb.set(p, bb.get(i));
            }
        }
        if (p > -1) {
            bb = new ArrayList<>(bb.subList(0, p + 1));
        }
        return p;

    }

    public int[] Sorting(int[] a) {
        int j, m, h, p, z, n1, n2, k1, k2;
        if (a[0] < 4) {
            return a;
        }
        z = a[0] / 2;
        j = z / 2;
        while (j > 0) {
            h = z - j;
            do {
                p = 0;
                for (m = 1; m <= h; m++) {
                    k1 = m + m - 1;
                    k2 = m + m + j + j - 1;
                    if (a[k1] > a[k2]) {
                        n1 = a[k1];
                        n2 = a[k1 + 1];
                        a[k1] = a[k2];
                        a[k1 + 1] = a[k2 + 1];
                        a[k2] = n1;
                        a[k2 + 1] = n2;
                        p = m;
                    }
                }
                h = p - j;
            } while (p > 0);
            j = j / 2;
        }
        return a;
    }

    public int[] joinTwo(int[] a, int[] b) {
        int j = 0;
        int[] joinedArray = Arrays.copyOf(a, a.length + b.length - 1);
        for (int i = 1; i < b[0]; i += 2) {
            if (b[i] > 0) {
                j++;
                joinedArray[a[0] + j] = b[i];
                j++;
                joinedArray[a[0] + j] = b[i + 1];
            }
        }
        joinedArray[0] = a[0] + j;
        joinedArray = Sorting(joinedArray);
        return Compress(joinedArray);
    }

    public int[] Compress(int[] a) {
        int x1, x2;
        int k = 0;

        for (int i = 1; i < a[0] - 2; i += 2) {
            if (a[i] > 0) {
                x1 = a[i] + Math.abs(a[i + 1]);

                for (int j = i + 2; j < a[0]; j += 2) {
                    if (a[j] > 0) {
                        if (x1 < a[j]) {
                            break;
                        }
                        x2 = a[j] + Math.abs(a[j + 1]);

                        if (a[i] >= a[j]) {
                            if (x2 > a[i]) {
                                if ((a[i + 1] > 0 && a[j + 1] > 0) || (a[i + 1] < 0 && a[j + 1] < 0)) {
                                    if (a[i] > a[j]) {
                                        a[i] = a[j];
                                    }
                                    if (x2 > x1) {
                                        if (a[i + 1] > 0) {
                                            a[i + 1] = x2 - a[i];
                                        } else {
                                            a[i + 1] = -(x2 - a[i]);
                                        }
                                    } else {
                                        if (a[i + 1] > 0) {
                                            a[i + 1] = x1 - a[i];
                                        } else {
                                            a[i + 1] = -(x1 - a[i]);
                                        }
                                    }
                                    a[j + 1] = 0;
                                }
                            }
                        } else if (a[i] < a[j]) {
                            if (x1 > a[j]) {
                                if ((a[i + 1] > 0 && a[j + 1] > 0) || (a[i + 1] < 0 && a[j + 1] < 0)) {
                                    if (x2 > x1) {
                                        if (a[i + 1] > 0) {
                                            a[i + 1] = x2 - a[i];
                                        } else {
                                            a[i + 1] = -(x2 - a[i]);
                                        }
                                    } else {
                                        if (a[i + 1] > 0) {
                                            a[i + 1] = x1 - a[i];
                                        } else {
                                            a[i + 1] = -(x1 - a[i]);
                                        }
                                    }
                                    a[j + 1] = 0;
                                }
                            }
                        }
                    }
                }
            }
        }

        int[] resultArray = new int[a.length];
        for (int i = 1; i < a[0]; i += 2) {
            if (a[i + 1] != 0) {
                k += 2;
                resultArray[k - 1] = a[i];
                resultArray[k] = a[i + 1];
            }
        }
        resultArray[0] = k;
        return Arrays.copyOf(resultArray, (int) (k + 1));
    }

    public boolean isCluster(int[] a, int[] b) {
        int x1, x2, d;
        for (int i = 1; i < a[0]; i += 2) {
            x1 = a[i] + Math.abs(a[i + 1]);
            for (int j = 1; j < b[0]; j += 2) {
                if (x1 < b[j]) {
                    break;
                }
                x2 = b[j] + Math.abs(b[j + 1]);
                d = 0;
                if (a[i] >= b[j]) {
                    if (x1 >= x2) {
                        d = x2 - a[i];
                    } else {
                        d = Math.abs(a[i + 1]);
                    }
                }
                if (a[i] < b[j]) {
                    if (x1 >= x2) {
                        d = x1 - b[j];
                    } else {
                        d = Math.abs(b[j + 1]);
                    }
                }

                if (d > minlenseq) {
                    return true;
                }
            }
        }
        return false;
    }

    private void PrintResult(int n) throws IOException {
        if (bb == null) {
            return;
        }
        int k = 0;
        int l = seq[n].length();
        if (numnonn == 0) {
            numnonn = l;
        }
        double z = numnonn;
        int v = numnonn;

        long duration = (System.nanoTime() - startTime) / 1000000000;

        String reportfile = filePath + "_" + (n + 1) + ".gff";
        String maskedfile = filePath + "_" + (n + 1) + ".msk";
        if (nseq == 1) {
            reportfile = filePath + ".gff";
            maskedfile = filePath + ".msk";
        }

        if (MaskedShow) {
            byte[] m = new byte[l];

            try (FileWriter fileWriter = new FileWriter(maskedfile)) {
                fileWriter.write(">" + sname[n]);
                System.out.println("Saving masked file: " + maskedfile);
                for (int i = 0; i < bb.size(); i++) {
                    int[] z7 = bb.get(i);
                    if (z7[0] > 1) {
                        for (int j = 1; j < z7[0]; j += 2) {
                            int x = z7[j] + Math.abs(z7[j + 1]) - 1;
                            for (int h = z7[j]; h < x; h++) {
                                m[h] = 4;
                            }
                        }
                    }
                }

                byte[] c = seq[n].getBytes();
                for (int i = 0; i < l; i++) {
                    if (m[i] == 0) {
                        z--;
                        c[i] = (byte) (c[i] - 32);
                    }
                }
                z = (z * 100 / v);
                fileWriter.write("Sequence coverage by repeats " + String.format("%.2f", z) + "%\n\n");
                fileWriter.write(new String(c));
                System.out.println("Sequence coverage by repeats " + String.format("%.2f", z) + "%");
            }
        }

        if (GFFShow) {
            StringBuilder sr = new StringBuilder();
            sr.append("kmer=").append(kmerln).append("\n").append("Minimal repeat=").append(minlenblock).append("\n").append("Repeat filter=").append(minlenseq).append("\nQuick analysis is false\n");
            sr.append("__________________________________________________\n Repeats search for: ").append(filePath).append("//").append(sname[n]).append(" ").append(l).append("bp :\n");
            sr.append("Time taken: ").append(duration).append(" seconds\n\n");
            if (SeqShow) {
                sr.append("Seqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\tSequence").append("\n");
            } else {
                sr.append("Seqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\tSequence").append("\n");
            }
            if (v < l) {
                double d = ((l - v) * 100) / l;
                sr.append("Sequence gap (bp)=").append(l - v).append(" (").append(String.format("%.3f", d)).append("%)\n");
            }

            try (FileWriter fileWriter = new FileWriter(reportfile)) {

                System.out.println("Saving report file: " + reportfile);
                if (MaskedShow) {
                    fileWriter.write("Sequence coverage by repeats " + String.format("%.2f", z) + "%\n\n");
                }
                fileWriter.write(sr.toString());
                for (int i = 0; i < bb.size(); i++) {
                    int[] z7 = bb.get(i);
                    if (z7[0] > 1) {
                        k++;
                        for (int j = 1; j < z7[0]; j += 2) {
                            if (Math.abs(z7[j + 1]) > minlenseq) {

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
                                if (z7[j + 1] > 0) {
                                    sr.append(sname[n]).append("\t").append(".").append("\t").append(k).append("\t").append(z7[j] + 1).append("\t").append(x + 1).append("\t").append(z7[j + 1]).append("\t").append(".").append("\t").append("+").append("\t").append(s0).append("\n");
//                                sr.append(k).append("\t").append(z7[j] + 1).append("\t").append(x + 1).append("\t").append(z7[j + 1]).append("\t").append("\t").append(s0).append("\n");
                                    fileWriter.write(sr.toString());
                                } else {
                                    sr.append(sname[n]).append("\t").append(".").append("\t").append(k).append("\t").append(-z7[j + 1]).append("\t").append(x + 1).append("\t").append(z7[j] - 1).append("\t").append(".").append("\t").append("-").append("\t").append(s0).append("\n");
//                               sr.append(k).append("\t").append(x + 1).append("\t").append(-z7[j] - 1).append("\t").append(-z7[j + 1]).append("\t").append(s0).append("\n");
                                    fileWriter.write(sr.toString());
                                }
                            }
                        }
                    }
                }
            }
        }
        /*
GFF format General Feature Format is a format for describing genes and other features associated with DNA, RNA and Protein sequences. GFF lines have nine tab-separated fields:
Generic Feature Format Version 3 (GFF3) https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
1. seqid - Must be a chromosome or scaffold or contig.
2. source - The program that generated this feature.
3. type - "repeat".
4. start - The starting position of the feature in the sequence. The first base is numbered 1.
5. stop - The ending position of the feature (inclusive).
6. score - length 
7. strand - Valid entries include '+', '-', or '.' (for don't know/care).
8. phase - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
9. attributes - All lines with the same group are linked together into a single item.
         */
    }

    private void PictureSave(int n, int dw, int dh) throws IOException {
        int maxClusters = 50;
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
            for (int j = 1; j < z7.length; j += 2) {
                int x1 = 50 + (int) (z7[j] * w1);
                int x2 = 50 + (int) ((z7[j] + z7[j + 1]) * w1);
                g2d.setColor(Color.DARK_GRAY);
                g2d.drawLine(x1, 22, x2, 22); // draw dark gray line
            }

            int y = 50 + (i * z);
            for (int j = 1; j < z7.length; j += 2) {
                if (Math.abs(z7[j + 1]) > minlenseq) {
                    int x1 = 50 + (int) (z7[j] * w1);
                    int x2 = 50 + (int) ((z7[j] + Math.abs(z7[j + 1])) * w1);
                    if (z7[j + 1] > 0) {
                        g2d.setColor(Color.BLUE);
                    } else {
                        g2d.setColor(Color.RED);
                    }
                    g2d.drawLine(x1, y, x2, y); //(int x1, int y1, int x2, int y2)
                }
            }
        }
    }

    private String filePath;
    public int nseq;
    public int iwidth = 0;
    public int iheight = 0;
    private String[] seq;
    private String[] sname;
    private int nblocks = 2;
    private int numnonn = 0;   // calculated non-N bases at the sequence
    private final int ssrlen = 30;
    private int minlenblock = 50;   // initial sequence length user control  
    private int minlenseq = 50;     // sequence length 
    private final int mnblock = 20;
    private int kmerln = 12;
    private int flanks = 20;
    private int telomers = 14; // Kmax=11 ->SSR  Kmax=14 -> telomers
    private boolean SeqShow;
    private boolean MaskedShow;
    private boolean GFFShow;
    private ArrayList<int[]> bb;
    private long startTime;
}
