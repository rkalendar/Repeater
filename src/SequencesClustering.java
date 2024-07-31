
import java.util.Arrays;
import java.util.HashSet;

public final class SequencesClustering {

    public SequencesClustering(String s, int[] x, int sim, int kmer) {
        if (sim < 60) {
            sim = 60;
        }
        if (kmer < 12) {
            kmer = 12;
        }
        if (kmer > 99) {
            kmer = 99;
        }
        this.sim = sim;
        nseq = x.length / 2;
        ncl = Clustering(s, kmer, x);
    }

    private int Clustering(String seq, int kmer, int[] x1) {
        int n = 0;
        cx = new int[nseq];      // clusters   
        d = new int[nseq][2];
        for (int j = 0; j < nseq; j++) {
            int p = j * 2;
            d[j][0] = x1[p];
            d[j][1] = x1[p + 1] - x1[p];
        }
        Arrays.sort(d, (int[] a, int[] b) -> {
            return Integer.compare(b[1], a[1]);
        });

        for (int j = 0; j < nseq; j++) {
            if (cx[j] == 0) {
                HashSet<String> m2 = new HashSet<>();
                n++;
                cx[j] = n;
                for (int i = 0; i < d[j][1] - kmer + 1; i++) {
                    String s = seq.substring(d[j][0] + i, d[j][0] + i + kmer);
                    if (m2.contains(s)) {
                    } else {
                        m2.add(s);
                    }
                    s = dna.ComplementDNA2(s);
                    if (m2.contains(s)) {
                    } else {
                        m2.add(s);
                    }
                }

// compare            
                for (int k = j + 1; k < nseq; k++) {
                    if (cx[k] == 0) {
                        int r = 0;
                        for (int i = 0; i < d[k][1] - kmer + 1; i++) {
                            String s = seq.substring(d[k][0] + i, d[k][0] + i + kmer);
                            if (m2.contains(s)) {
                                r++;
                            }
                            s = dna.ComplementDNA2(s);//s = a.substring(i, i + kmer);
                            if (m2.contains(s)) {
                                r++;
                            }
                        }
                        int v = (100 * r) / d[k][1];
                        if (v > sim) {
                            cx[k] = n;
                        }
                    }
                }
            }
        }
        return n;
    }

    public int getNcl() {
        return ncl;
    }

    public int[] Result() {
        return cx;
    }

    public int[][] ResultArray() {
        return d;
    }
    private final int nseq;
    private final int sim;
    private final int ncl;
    private int[] cx;
    private int[][] d;
}
