
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public final class SequencesClustering {

    public SequencesClustering(String s, int[] x, int sim, int kmer, boolean slow) {
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
        if (slow) {
            ncl = ClusteringSlow(s, kmer, x);
        } else {
            ncl = ClusteringPattern(s, x);
        }
    }

    private int ClusteringPattern(String seq, int[] x1) {
        HashMap<String, Integer> pt = new HashMap<>();
        pt.put("aaatt", 0);
        pt.put("aagtt", 1);
        pt.put("ataat", 2);
        pt.put("atgat", 3);
        pt.put("acagt", 4);
        pt.put("acggt", 5);
        pt.put("agact", 6);
        pt.put("aggct", 7);
        pt.put("ttaaa", 8);
        pt.put("ttgaa", 9);
        pt.put("taata", 10);
        pt.put("tagta", 11);
        pt.put("tgaca", 12);
        pt.put("tggca", 13);
        pt.put("tcaga", 14);
        pt.put("tcgga", 15);
        pt.put("ccagg", 16);
        pt.put("ccggg", 17);
        pt.put("caatg", 18);
        pt.put("cagtg", 19);
        pt.put("ctaag", 20);
        pt.put("ctgag", 21);
        pt.put("cgacg", 22);
        pt.put("cggcg", 23);
        pt.put("ggacc", 24);
        pt.put("gggcc", 25);
        pt.put("gaatc", 26);
        pt.put("gagtc", 27);
        pt.put("gtaac", 28);
        pt.put("gtgac", 29);
        pt.put("gcagc", 30);
        pt.put("gcggc", 31);

        d = new int[nseq][2];
        for (int j = 0; j < nseq; j++) {
            int p = j * 2;
            d[j][0] = x1[p];
            d[j][1] = x1[p + 1] - x1[p];
        }
        Arrays.sort(d, (int[] a, int[] b) -> {
            return Integer.compare(b[1], a[1]);
        });
        
        int nkmers = pt.size();
        int kmer = 5;
        int[][] m2 = new int[nseq][nkmers];
        for (int j = 0; j < nseq; j++) {
            for (int i = 0; i < d[j][1] - kmer + 1; i++) {
                String s = seq.substring(d[j][0] + i, d[j][0] + i + kmer);
                if (pt.containsKey(s)) {
                    m2[j][pt.get(s)]++;
                }
                s = dna.ComplementDNA2(s);
                if (pt.containsKey(s)) {
                    m2[j][pt.get(s)]++;
                }
            }
        }

        // compare sequences  
        int n = 0;              // amount of clusters
        cx = new int[nseq];     // clusters           
        for (int i = 0; i < nseq; i++) {
            if (cx[i] == 0) {
                n++;
                cx[i] = n;

                for (int j = i + 1; j < nseq; j++) {
                    if (cx[j] == 0) {
                        int[] m = new int[nkmers + 1];

                        for (int k = 0; k < nkmers; k++) {
                            if (m2[i][k] > 0 && m2[j][k] > 0) {
                                if (m[0] == 0) {
                                    m[0] = 1;
                                    m[m[0]] = k;
                                } else {
                                    m[0]++;
                                    m[m[0]] = k;
                                }
                            }
                        }

                        int v = 0;// pricatical matches 
                        int z = 0;// theoretically maximum possible matches  
                        for (int k = 1; k < 1 + m[0]; k++) {
                            for (int y = k + 1; y < 1 + m[0]; y++) {
                                double di, dj;
                                z++;
                                if (m2[i][m[k]] < m2[i][m[y]]) {
                                    di = (double) m2[i][m[k]] / m2[i][m[y]];
                                    dj = (double) m2[j][m[k]] / m2[j][m[y]];
                                } else {
                                    di = (double) m2[i][m[y]] / m2[i][m[k]];
                                    dj = (double) m2[j][m[y]] / m2[j][m[k]];
                                }
                                if (di == dj) {
                                    v++;
                                }
                                if (di > dj) {
                                    if (di <= dj + (dj * dif)) {
                                        v++;
                                    }
                                }
                                if (di < dj) {
                                    if (dj <= di + (di * dif)) {
                                        v++;
                                    }
                                }
                            }
                            if (v > 0 && v + (v * dif) > z) {
                                cx[j] = n;
                            }
                        }
                    }
                }
            }
        }
        return n;
    }

    private int ClusteringSlow(String seq, int kmer, int[] x1) {
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
                            s = dna.ComplementDNA2(s);
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
    private final double dif = 0.4d;
    private final int sim;
    private final int ncl;
    private int[] cx;
    private int[][] d;
}
