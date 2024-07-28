
public final class SequencesClustering2 {

    public SequencesClustering2(String seq, int[] z, int sim, int minlenseq) {

        if (sim < 60) {
            sim = 60;
        }
        if (sim > 90) {
            sim = 90;
        }
        this.sim = sim;
        nseq = z.length / 2;
        Lclusters = new int[nseq];

        ncl = Clustering12(seq, z);
    }

    private int Clustering12(String seq, int[] z2) {
        int x = 0;
        int l = 0;
        int n = 0;
        int w = 0;
        int z = 0;
        int k = 0;
        int f = 0;
        int h = 0;
        int kmer = 12;

        int[][][][][][][][][][][][][] px = new int[5][5][5][5][5][5][5][5][5][5][5][5][2];

        int[] lx = new int[nseq];
        int[] z3 = new int[kmer];

        for (int j = 0; j < nseq; j++) {
            String taq = seq.substring(z2[j * 2], z2[j * 2 + 1]);
            lx[j] = taq.length();
            l = l + lx[j];
            if (lx[j] > kmer) {
                z3[0] = tables.dx[taq.charAt(0)];
                z3[1] = tables.dx[taq.charAt(1)];
                z3[2] = tables.dx[taq.charAt(2)];
                z3[3] = tables.dx[taq.charAt(3)];
                z3[4] = tables.dx[taq.charAt(4)];
                z3[5] = tables.dx[taq.charAt(5)];
                z3[6] = tables.dx[taq.charAt(6)];
                z3[7] = tables.dx[taq.charAt(7)];
                z3[8] = tables.dx[taq.charAt(8)];
                z3[9] = tables.dx[taq.charAt(9)];
                z3[10] = tables.dx[taq.charAt(10)];
                for (int i = kmer - 1; i < lx[j]; i++) {
                    z3[kmer - 1] = tables.dx[taq.charAt(i)];
                    if (px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] < j + 1) {
                        px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] = j + 1;
                    }
                    px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]++;
                    z3[0] = z3[1];
                    z3[1] = z3[2];
                    z3[2] = z3[3];
                    z3[3] = z3[4];
                    z3[4] = z3[5];
                    z3[5] = z3[6];
                    z3[6] = z3[7];
                    z3[7] = z3[8];
                    z3[8] = z3[9];
                    z3[9] = z3[10];
                    z3[10] = z3[11];
                }
            }
        }

        for (int j = 0; j < nseq; j++) {
            if (lx[j] > kmer) {
                String taq = seq.substring(z2[j * 2], z2[j * 2 + 1]);
                String a2 = dna.ComplementDNA(taq);
                z3[0] = tables.dx[a2.charAt(0)];
                z3[1] = tables.dx[a2.charAt(1)];
                z3[2] = tables.dx[a2.charAt(2)];
                z3[3] = tables.dx[a2.charAt(3)];
                z3[4] = tables.dx[a2.charAt(4)];
                z3[5] = tables.dx[a2.charAt(5)];
                z3[6] = tables.dx[a2.charAt(6)];
                z3[7] = tables.dx[a2.charAt(7)];
                z3[8] = tables.dx[a2.charAt(8)];
                z3[9] = tables.dx[a2.charAt(9)];
                z3[10] = tables.dx[a2.charAt(10)];

                n = nseq + j;
                for (int i = kmer - 1; i < lx[j]; i++) {
                    z3[kmer - 1] = tables.dx[a2.charAt(i)];
                    if (px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] < n) {
                        px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] = n;
                    }
                    px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]++;
                    z3[0] = z3[1];
                    z3[1] = z3[2];
                    z3[2] = z3[3];
                    z3[3] = z3[4];
                    z3[4] = z3[5];
                    z3[5] = z3[6];
                    z3[6] = z3[7];
                    z3[7] = z3[8];
                    z3[8] = z3[9];
                    z3[9] = z3[10];
                    z3[10] = z3[11];
                }
            }
        }

        int n0 = l + l;
        int[] u1 = new int[n0];
        int[] u2 = new int[n0];
        for (int j = 0; j < nseq; j++) {
            if (lx[j] > kmer) {
                String taq = seq.substring(z2[j * 2], z2[j * 2 + 1]);
                z3[0] = tables.dx[taq.charAt(0)];
                z3[1] = tables.dx[taq.charAt(1)];
                z3[2] = tables.dx[taq.charAt(2)];
                z3[3] = tables.dx[taq.charAt(3)];
                z3[4] = tables.dx[taq.charAt(4)];
                z3[5] = tables.dx[taq.charAt(5)];
                z3[6] = tables.dx[taq.charAt(6)];
                z3[7] = tables.dx[taq.charAt(7)];
                z3[8] = tables.dx[taq.charAt(8)];
                z3[9] = tables.dx[taq.charAt(9)];
                z3[10] = tables.dx[taq.charAt(10)];

                for (int i = kmer - 1; i < lx[j]; i++) {
                    z3[kmer - 1] = tables.dx[taq.charAt(i)];
                    if (px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0] > 0) {
                        px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] = w;
                        u1[w] = i - kmer + 1;
                        u2[w] = j;
                        w = w + px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0];
                        px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0] = -1;
                    } else {

                        for (h = px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1]; h < px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] - px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]; h++) {
                            if (j == u2[h]) {
                                h = -1;
                                break;
                            }
                        }
                        if (h > -1) {
                            u1[px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] - px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]] = i - kmer + 1;
                            u2[px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] - px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]] = j;
                            px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]--;
                        }
                    }
                    z3[0] = z3[1];
                    z3[1] = z3[2];
                    z3[2] = z3[3];
                    z3[3] = z3[4];
                    z3[4] = z3[5];
                    z3[5] = z3[6];
                    z3[6] = z3[7];
                    z3[7] = z3[8];
                    z3[8] = z3[9];
                    z3[9] = z3[10];
                    z3[10] = z3[11];
                }
            }
        }
        for (int j = 0; j < nseq; j++) {
            if (lx[j] > kmer) {
                String taq = seq.substring(z2[j * 2], z2[j * 2 + 1]);
                String a2 = dna.ComplementDNA(taq);
                z3[0] = tables.dx[a2.charAt(0)];
                z3[1] = tables.dx[a2.charAt(1)];
                z3[2] = tables.dx[a2.charAt(2)];
                z3[3] = tables.dx[a2.charAt(3)];
                z3[4] = tables.dx[a2.charAt(4)];
                z3[5] = tables.dx[a2.charAt(5)];
                z3[6] = tables.dx[a2.charAt(6)];
                z3[7] = tables.dx[a2.charAt(7)];
                z3[8] = tables.dx[a2.charAt(8)];
                z3[9] = tables.dx[a2.charAt(9)];
                z3[10] = tables.dx[a2.charAt(10)];

                n = nseq + j;
                for (int i = kmer - 1; i < lx[j]; i++) {
                    z3[kmer - 1] = tables.dx[a2.charAt(i)];
                    if (px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0] > 0) {
                        px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] = w;
                        u1[w] = i - kmer + 1;
                        u2[w] = n;
                        w = w + px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0];
                        px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0] = -1;
                    } else {
                        for (h = px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1]; h < px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] - px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]; h++) {
                            if (j == u2[h]) {
                                h = -1;
                                break;
                            }
                        }
                        if (h > -1) {
                            u1[px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] - px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]] = i - kmer + 1;
                            u2[px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] - px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]] = n;
                            px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]--;
                        }
                    }
                    z3[0] = z3[1];
                    z3[1] = z3[2];
                    z3[2] = z3[3];
                    z3[3] = z3[4];
                    z3[4] = z3[5];
                    z3[5] = z3[6];
                    z3[6] = z3[7];
                    z3[7] = z3[8];
                    z3[8] = z3[9];
                    z3[9] = z3[10];
                    z3[10] = z3[11];
                }
            }
        }

        k = 0;
        int p = nseq * nseq * 2;
        int[] u3 = new int[p]; // similarity between seq1 & seq2
        int[] u4 = new int[p]; // i seq1
        int[] u5 = new int[p]; // j seq2
        for (int j = 0; j < nseq - 1; j++) {
            k = k + nseq + nseq - j - 1;
            if (lx[j] > kmer) {
                String taq = seq.substring(z2[j * 2], z2[j * 2 + 1]);
                z3[0] = tables.dx[taq.charAt(0)];
                z3[1] = tables.dx[taq.charAt(1)];
                z3[2] = tables.dx[taq.charAt(2)];
                z3[3] = tables.dx[taq.charAt(3)];
                z3[4] = tables.dx[taq.charAt(4)];
                z3[5] = tables.dx[taq.charAt(5)];
                z3[6] = tables.dx[taq.charAt(6)];
                z3[7] = tables.dx[taq.charAt(7)];
                z3[8] = tables.dx[taq.charAt(8)];
                z3[9] = tables.dx[taq.charAt(9)];
                z3[10] = tables.dx[taq.charAt(10)];

                for (int i = kmer - 1; i < lx[j]; i++) {
                    z3[kmer - 1] = tables.dx[taq.charAt(i)];
                    for (h = 1; h <= -px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][0]; h++) {
                        x = u2[px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] + h - 1];
                        if (u2[px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] + h - 1] == j && u1[px[z3[0]][z3[1]][z3[2]][z3[3]][z3[4]][z3[5]][z3[6]][z3[7]][z3[8]][z3[9]][z3[10]][z3[11]][1] + h - 1] < i - kmer + 1) {
                            break;
                        }
                        if (x > j) {
                            z = x + k - nseq - nseq;
                            u3[z]++;
                            u4[z] = j;
                            u5[z] = x;
                        }
                    }
                    z3[0] = z3[1];
                    z3[1] = z3[2];
                    z3[2] = z3[3];
                    z3[3] = z3[4];
                    z3[4] = z3[5];
                    z3[5] = z3[6];
                    z3[6] = z3[7];
                    z3[7] = z3[8];
                    z3[8] = z3[9];
                    z3[9] = z3[10];
                    z3[10] = z3[11];
                }
            }
        }

        for (int j = 1; j < p; j++) {
            if (u3[j] > 0) {
                int m1 = (u4[j] < nseq) ? u4[j] : u4[j] - nseq;
                int m2 = (u5[j] < nseq) ? u5[j] : u5[j] - nseq;
                int m = (100 * (kmer - 1 + u3[j])) / Math.min(lx[m1], lx[m2]);

                if (m >= sim) {
                    if (Lclusters[m1] + Lclusters[m2] == 0) {
                        r1++;
                        ncl++;
                        Lclusters[m1] = ncl;
                        Lclusters[m2] = ncl;
                    } else {
                        if (Lclusters[m1] > Lclusters[m2] && Lclusters[m2] > 0) {
                            r1--;
                            f = Lclusters[m1];
                            for (int u = 0; u < nseq; u++) {
                                if (Lclusters[u] == f) {
                                    Lclusters[u] = Lclusters[m2];
                                }
                            }
                        } else if (Lclusters[m1] < Lclusters[m2] && Lclusters[m1] > 0) {
                            r1--;
                            f = Lclusters[m2];
                            for (int u = 0; u < nseq; u++) {
                                if (Lclusters[u] == f) {
                                    Lclusters[u] = Lclusters[m1];
                                }
                            }
                        } else if (Lclusters[m1] == 0) {
                            Lclusters[m1] = Lclusters[m2];
                        } else if (Lclusters[m2] == 0) {
                            Lclusters[m2] = Lclusters[m1];
                        }
                    }
                }
            }
        }
        return ncl;
    }

    public int getNcl() {
        return ncl;
    }

    public int getRNcl() {
        return r1;
    }

    public int[] getClusters() {
        return Lclusters;
    }

    private final int nseq;
    private final int[] Lclusters;
    private final int sim;
    private int r1;
    private int ncl;

}
