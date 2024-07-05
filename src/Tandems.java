
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

public final class Tandems {

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
        if (kmerln < 5) {
            kmerln = 5;
        }
        minlenblock = kmerln + kmerln;
        minlenseq = minlenblock;
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

    public void SetTandemLen(int i) {
        minlenblock = i;
        if (minlenblock < mnblock) {
            minlenblock = mnblock;
        }
        minlenseq = minlenblock;
    }

    public void SetSequenceLen(int i) {
        minlenseq = i;
        if (minlenseq < minlenblock) {
            minlenseq = minlenblock;
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
            FindAllSSRs(seq[i], telomers);
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

    private byte[] SimpleRepeatsMasking(byte b[], int ssr) {
        int n1;
        int n2;
        int n3;
        int x;
        int e;
        int lmer = 17;
        int Kmax = ssr;   // Kmax=11 for SSR  Kmax=14 for telomers
        int l = b.length;
        if (l < lmer + 3) {
            return b;
        }

        byte[][] v2 = new byte[5][5];
        byte[][][] v3 = new byte[5][5][5];
        byte[] k2 = new byte[l - lmer];

        for (int i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];
        }

        for (int i = 0; i < lmer; i++) { //max=18
            n1 = b[i];        //16
            n2 = b[i + 1];    //17
            v2[n1][n2]++;
            if (v2[n1][n2] == 1) {
                k2[0]++;
            }
        }

        for (int i = 0; i < lmer - 1; i++) {
            n1 = b[i];       //15
            n2 = b[i + 1];   //16
            n3 = b[i + 2];   //17
            v3[n1][n2][n3]++;
            if (v3[n1][n2][n3] == 1) {
                k2[0]++;
            }
        }

        for (int i = 1; i < l - lmer; i++) {
            k2[i] = k2[i - 1];
            n1 = b[i - 1];
            n2 = b[i];
            n3 = b[i + 1];
            v2[n1][n2]--;
            if (v2[n1][n2] == 0) {
                k2[i]--;
            }

            v3[n1][n2][n3]--;
            if (v3[n1][n2][n3] == 0) {
                k2[i]--;
            }

            n1 = b[i + lmer - 2]; //16
            n2 = b[i + lmer - 1]; //17
            n3 = b[i + lmer];     //18
            v2[n2][n3]++;
            if (v2[n2][n3] == 1) {
                k2[i]++;
            }

            v3[n1][n2][n3]++;
            if (v3[n1][n2][n3] == 1) {
                k2[i]++;
            }
        }

        for (int i = 0; i < l - lmer; i++) {
            if (k2[i] < Kmax && b[i] < 4) {
                x = i;
                for (e = i + 1; e < l - lmer; e++) {
                    if (k2[e] > Kmax || b[e] == 4) {
                        break;
                    }
                }
                i = e - 1;
                if (e - x > ssrlen) {
                    if (x > 0) {
                        x--;
                    }
                    for (int h = x; h < e + lmer; h++) {
                        b[h] = 5;
                    }
                }
            }
        }
        return b;
    }

    private int FindAllSSRs(String seq, int ssr) {
        byte b[] = seq.getBytes();
        b = SimpleRepeatsMasking(b, ssr);
        bb = new ArrayList<>();
        int[] k7;
        int e;
        int h = b.length - 1;
        for (int i = 0; i <= h; i++) {
            if (b[i] == 5) {
                for (e = i + 1; e <= h; e++) {
                    if (b[e] < 5) {
                        break;
                    }
                }
                if (minlenseq < e - i) {
                    k7 = new int[3];
                    k7[0] = 2;
                    k7[1] = i;
                    k7[2] = e - i;
                    bb.add(k7);
                }
                i = e;
            }
        }

        return bb.size();
    }

    private int FindAllRepeats(String seq, int kmer) {
        int l = seq.length();
        int g = l + l + 1;
        int k = 0;
        int n = 0;
        int t = 0;
        int i = 0;
        int h = 0;
        int e = 0;
        int y = 0;
        int z = 0;
        int w = 0;
        int x = 0;
        int q = 0;
        int j = 0;
        int r = 0;
        int p = 0;

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
        double z = numnonn;
        int l = seq[n].length();
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
                char[] c = seq[n].toCharArray();
                for (int i = 0; i < l; i++) {
                    if (m[i] == 0) {
                        z--;
                        c[i] = (char) (c[i] - 32);
                    }
                }
                z = (z * 100 / v);
                fileWriter.write("Sequence coverage by repeats " + String.format("%.2f", z) + "%\n\n");
                fileWriter.write(new String(c));
                System.out.println("Sequence coverage by repeats " + String.format("%.2f", z) + "%");
            }
        }

        StringBuilder sr = new StringBuilder();
        sr.append("kmer=").append(kmerln).append("\n").append("Initial string length filter=").append(minlenblock).append("\n").append("String length filter=").append(minlenseq).append("\n").append("Quick analysis is false\n");
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
        if (bb == null) {
            return;
        }
        int k = 100;
        int b = bb.size();
        int l = seq[n].length();
        int width = l < 5_000_000 ? l / 100 : 5_000 + (l - 5_000) / 200; // image=10000x3000
        if (dw > 0) {
            width = dw;
        }
        if (width > 46340) { //too big a picture leads to an error  
            width = 46340; //https://jobcardsystems.com/index.php/blog/29-effective-handling-of-big-images-in-java
        }
        if (width < 500) {
            width = 500;
        }
        int height = b < 5000 ? k + b : 5100 + (b - 5000) / 200;        //too big a picture leads to an error
        if (dh > 0) {
            height = dh;
        }
        if (height > 46340) {
            height = 46340;
        }
        if (height < 500) {
            height = 500;
        }
        float dotSize = 1 + (b / 500);                                  //7.0f;
        double w1 = (double) width / l;                                 // nucleotides per pixel        
        if (dotSize < 1) {
            dotSize = 1.0f;
        }
        if (dotSize > 10) {
            dotSize = 10.0f;
        }

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
        for (int i = 0; i < 11; i++) {
            g2d.drawLine(i * w + 50, 5, i * w + 50, 20);
            g2d.drawString(String.format("%,d", (1 + (int) (i * d))), 55 + i * w, 10);
        }

        for (int i = 0; i < b; i++) {
            int[] z7 = bb.get(i);
            if (z7[0] > 1) {
                int y = (i + 70);// y1-y2 line                
                for (int j = 1; j < z7[0]; j += 2) {
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
    private int telomers = 11; // Kmax=11 ->SSR  Kmax=14 -> telomers
    private boolean SeqShow;
    private boolean MaskedShow;
    private ArrayList<int[]> bb;
    private long startTime;
}
