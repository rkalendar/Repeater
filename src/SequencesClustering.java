
import java.util.HashMap;

public final class SequencesClustering {

    public SequencesClustering(String s, int[] x, int sim, boolean t, int kmer) {
        this.seq = s;
        this.x1 = x;
        if (sim < 60) {
            sim = 60;
        }
        if (sim > 90) {
            sim = 90;
        }
        if (kmer < 12) {
            kmer = 12;
        }
        if (kmer > 99) {
            kmer = 99;
        }
        this.sim = sim;
        nseq = x.length / 2;
        ncl = Clustering(kmer);
    }

    public int[] Result() {
        return cx;
    }

    private int Clustering(int kmer) {
        int n = 0;
        int[] lx = new int[nseq];
        HashMap<String, int[]> map = new HashMap<>();
        HashMap<String, int[]> m = new HashMap<>();

        for (int j = 0; j < nseq; j++) {
            int p = j * 2;
            //     String a = dna.ComplementDNA(seq.substring(x1[p], x1[p + 1]));
            lx[j] = x1[p + 1] - x1[p];
            for (int i = 0; i < lx[j] - kmer + 1; i++) {
                String s = seq.substring(x1[p] + i, x1[p] + i + kmer);
                if (m.containsKey(s)) {
                    int[] t = m.get(s);
                    if (t[1] < (j + 1)) {
                        t[0]++;       // number of seq[] 
                        t[1] = j + 1; // index seq[]+1 
                        m.replace(s, t);
                    }
                } else {
                    int[] t = new int[2];
                    t[0] = 1;
                    t[1] = j + 1;
                    m.put(s, t);
                }

                s = dna.ComplementDNA2(s);//s = a.substring(i, i + kmer);
                if (m.containsKey(s)) {
                    int[] t = m.get(s);
                    if (t[1] < (j + 1)) {
                        t[0]++;       // number of seq[] 
                        t[1] = j + 1; // index seq[]+1 
                        m.replace(s, t);
                    }
                } else {
                    int[] t = new int[2];
                    t[0] = 1;
                    t[1] = j + 1;
                    m.put(s, t);
                }
            }

        }

        for (int j = 0; j < nseq; j++) {
            int p = j * 2;
            for (int i = 0; i < lx[j] - kmer + 1; i++) {
                String s = seq.substring(x1[p] + i, x1[p] + i + kmer);
                if (map.containsKey(s)) {
                    int[] y = map.get(s);
                    if (y[0] == j) {
                        y[y[1]]++;
                    } else {
                        y[0] = j;
                        y[y[1] + 1] = j;
                        y[y[1] + 2] = 1;
                    }
                    map.replace(s, y);
                } else {
                    int[] t = m.get(s);
                    int[] y = new int[t[0] + t[0] + 2];
                    y[0] = j; // last seq[]
                    y[1] = 3; // last x[i]
                    y[2] = j; // data 1= seq[]
                    y[3] = 1; // data 2= n->kmer at seq[]
                    map.put(s, y);
                }

                // s = a.substring(i, i + kmer); // memory problem - Exception in thread "main" java.lang.OutOfMemoryError: Java heap space
                s = dna.ComplementDNA2(s);//s = a.substring(i, i + kmer);
                if (map.containsKey(s)) {
                    int[] y = map.get(s);
                    if (y[0] == j) {
                        y[y[1]]++;
                    } else {
                        y[0] = j;
                        y[y[1] + 1] = j;
                        y[y[1] + 2] = 1;
                    }
                    map.replace(s, y);
                } else {
                    int[] t = m.get(s);
                    int[] y = new int[t[0] + t[0] + 2];
                    y[0] = j; // last seq[]
                    y[1] = 3; // last x[i]
                    y[2] = j; // data 1= seq[]
                    y[3] = 1; // data 2= n->kmer at seq[]
                    map.put(s, y);
                }
            }
        }
        m.clear();

        //clustering
        cx = new int[nseq];      // clusters        
        for (int j = 0; j < nseq; j++) {
            int p = j * 2;
            //    String a = dna.ComplementDNA(seq.substring(x1[p], x1[p + 1]));
            int[] x = new int[nseq];
            for (int i = 0; i < lx[j] - kmer + 1; i++) {
                String s = seq.substring(x1[p] + i, x1[p] + i + kmer);
                int[] t = map.get(s);
                for (int h = 2; h < t.length; h = h + 2) {
                    if (t[h] == j) {
                        x[t[h]]++;
                    } else {
                        x[t[h]] = x[t[h]] + t[h + 1];
                    }
                }

                // s = a.substring(i, i + kmer);
                s = dna.ComplementDNA2(s);
                t = map.get(s);
                for (int h = 2; h < t.length; h = h + 2) {
                    if (t[h] == j) {
                        x[t[h]]++;
                    } else {
                        x[t[h]] = x[t[h]] + t[h + 1];
                    }
                }
            }

            for (int h = 0; h < nseq; h++) {
                if (x[h] > 0 && x[j] > 0) {
                    int v = (100 * x[h]) / Math.min(lx[j], lx[h]);
                    if (v > sim) {
                        if (cx[j] == 0) {
                            if (cx[h] > 0) {
                                cx[j] = cx[h];
                            } else {
                                n++;
                                cx[j] = n;
                                cx[h] = n;
                            }
                        } else {
                            if (cx[h] > 0) {
                                if (cx[h] < cx[j]) {
                                    cx[j] = cx[h];
                                } else {
                                    cx[h] = cx[j];
                                }
                            } else {
                                cx[h] = cx[j];
                            }
                        }
                    }
                }
            }

        }
        n++;
        for (int h = 0; h < nseq; h++) {
            if (cx[h] == 0) {
                cx[h] = n;
            }
        }
        return n;
    }

    public int getNcl() {
        return ncl;
    }
    private final int nseq;
    private final int sim;
    private int ncl;
    private int[] cx;
    private int[] x1;
    private String seq;
}
