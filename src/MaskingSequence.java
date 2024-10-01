
import java.util.ArrayList;
import java.util.HashMap;

public class MaskingSequence {

    private long gapslen = 0;
    private long repeatslen = 0;

    public int[] Mask(String seq, byte[] ssrmsk, int kmer, int minlenseq) {
        int[] u = Masking(seq, ssrmsk, kmer);

        // identification blocks
        ArrayList<Integer> z = new ArrayList<>();

        for (int i = 0; i < u.length; i++) {
            if (u[i] > 1) {
                int x1 = i;
                int x2 = i;
                int x = i + 1;
                if (x1 > 0) {
                    x1 = x1 - 1;
                }

                for (int f = x; f < u.length; f++) {
                    x2 = f;
                    if (u[f] < 2) {
                        i = x2;
                        break;
                    }
                }

                if (x2 > x1) {
                    if (z.isEmpty()) {
                        z.add(x1);
                        z.add(x2);
                    } else {
                        int h = z.size() - 1;
                        int x4 = z.get(h);
                        if (x4 + kmer > x1) {
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
        u = new int[n];   //x1.....x2 block
        n = -1;
        for (int i = 0; i < z.size(); i += 2) {
            int x1 = z.get(i);
            int x2 = z.get(i + 1);
            if (x2 - x1 > minlenseq) {
                u[++n] = x1;         // start x1
                u[++n] = x2 - x1 + 1;// length 
                repeatslen += (x2 - x1 + 1);
            }
        }
        return u;
    }

    private int[] Masking(String seq, byte[] ssrmsk, int kmer) {
        HashMap<String, Integer> map = new HashMap<>();
        String aseq = dna.ComplementDNA(seq);

        int l = seq.length();
        int[] ax = new int[kmer];
        int[] bx = new int[5];

        byte b[] = seq.getBytes();
        byte c[] = aseq.getBytes();

        gapslen = 0;
        for (int i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];
            c[i] = tables.dx2[c[i]];
            if (b[i] == 4) {
                gapslen++;
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
                        if (ssrmsk[j] == 0) {
                            u[j]++;
                        }
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
                        for (int j = l - x - kmer; j < l - x; j++) {
                            if (ssrmsk[j] == 0) {
                                u[j]++;
                            }
                        }
                    }
                }
            }
            bx[c[i + 1 - kmer]]--;
            for (int j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }
        return u;
    }

    public long GapsLength() {
        return gapslen;
    }

    public long RepeatLength() {
        return repeatslen;
    }

}
