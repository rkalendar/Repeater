
public final class dna {

    public static String AntisenseDNA(String source) {
        byte[] b = source.getBytes();
        for (int i = 0; i < source.length(); i++) {
            if (tables.cdna[b[i]] > 0) {
                if (b[i] < 97) {
                    b[i] = (byte) (b[i] + 32);
                }
                b[i] = tables.cdna[b[i]];
            }
        }
        return new String(b);
    }

    public static String ReverseSeq(String source) {
        return new StringBuilder(source).reverse().toString();
    }

    public static String ReverseSeq(char[] source) {
        StringBuilder s = new StringBuilder();
        return s.append(source).reverse().toString();
    }

    public static String ComplementDNA(String source) {
        return new StringBuilder(AntisenseDNA(source)).reverse().toString();
    }

    public static String ComplementDNA2(String source) {
        byte[] b = source.getBytes();
        int l = source.length();
        int n = l / 2;
        for (int i = 0; i < n; i++) {
            if (tables.cdna[b[i]] > 0) {
                byte t = tables.cdna[b[l - i - 1]];
                b[l - i - 1] = tables.cdna[b[i]];
                b[i] = t;
            }
        }
        if ((l % 2) == 1) {
            if (tables.cdna[b[n]] > 0) {
                b[n] = tables.cdna[b[n]];
            }
        }
        return new String(b);
    }
}
