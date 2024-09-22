
import java.util.ArrayList;

public class LowComplexitySequence {

    private byte[] SimpleRepeatsMasking(byte b[], int telomer) {
        int n1;
        int n2;
        int n3;
        int x;
        int e;
        int lmer = 17;
        int Kmax = telomer;   // Kmax=11 for SSR  Kmax=14 for telomers
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

    public void FindAllSSRs(String seq, int telomer) {
        byte b[] = seq.getBytes();
        mapb = new byte[seq.length()];
        b = SimpleRepeatsMasking(b, telomer);
        ArrayList<int[]> b1 = new ArrayList<>();
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
                    k7 = new int[2];
                    k7[0] = i;
                    k7[1] = e - i;
                    b1.add(k7);
                }
                i = e;
            }
        }

        ibloks = new int[1 + 2 * b1.size()];
        e = 0;
        ibloks[0] = 2 * b1.size();
        for (int i = 0; i < b1.size(); i++) {
            int[] z7 = b1.get(i);
            ibloks[++e] = z7[0];
            ibloks[++e] = z7[1];
            for (int j = z7[0]; j < z7[0] + z7[1]; j++) {
                mapb[j] = 4;
            }
        }
    }

    public ArrayList Blocks() {
        ArrayList<int[]> bb = new ArrayList<>();
        bb.add(ibloks);
        return bb;
    }

    public byte[] MapBytes() {
        return mapb;
    }

    public int[] IntBlocks() {
        return ibloks;
    }
    private byte[] mapb;
    private int[] ibloks;
    private final int ssrlen = 30;
    private final int minlenseq = 50;     // sequence length 
}
