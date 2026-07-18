
public final class dna {

    public static String ComplementDNA2(String source) {
        // DNA substrings are ASCII; ISO-8859-1 is an exact, allocation-light round trip.
        byte[] b = source.getBytes(java.nio.charset.StandardCharsets.ISO_8859_1);
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
        return new String(b, java.nio.charset.StandardCharsets.ISO_8859_1);
    }
}
