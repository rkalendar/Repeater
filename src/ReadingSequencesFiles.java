public final class ReadingSequencesFiles {

    public ReadingSequencesFiles(byte[] s) {
        source = s;
        ReadingSequences();
    }

    public ReadingSequencesFiles(String s) {
        source = s.getBytes();
        ReadingSequences();
    }

    public String[] getSequences() {
        if (ns == 0) {
            return null;
        }
        return sequence;
    }

    public String[] getNames() {
        if (ns == 0) {
            return null;
        }
        return name_seq;
    }

    public int getNseq() {
        return ns;
    }

    private void ReadingSequences() {
        if (source == null) {
            return;
        }

        int l = source.length;
        int s = 0; // total length
        ns = 0;    // amount fasta sequences

        for (int i = 0; i < l; i++) {
            if (source[i] < 9) {
                ns = 0;
                break;
            }
            if (tables.cdn[source[i]] > 0) {
                s++;
            }
            if (source[i] == 62) {
                ns++;
            }
        }
        if (s == 0 || ns == 0) {
            return;
        }
        name_seq = new String[ns];
        sequence = new String[ns];
        // One reusable scratch buffer. The first pass counted `s` = every valid base
        // byte in the whole file (headers included), so `s` is >= any single sequence's
        // compacted length. Compacting straight from `source` into this buffer removes
        // the two Arrays.copyOfRange allocations per sequence (the raw slice and the
        // compacted slice), leaving only each String's own backing storage.
        final byte[] buf = new byte[s];
        int n = -1;
        int t = 0;

        for (int i = 0; i < l; i++) {
            if (source[i] == 62) {
                if (t > 0) {
                    int x = 0;
                    for (int j = t; j < i - 1; j++) {   // exclusive end i-1, as before
                        byte cc = tables.cdn[source[j]];
                        if (cc > 0) {
                            buf[x++] = cc;
                        }
                    }
                    sequence[n] = new String(buf, 0, x, java.nio.charset.StandardCharsets.ISO_8859_1);
                    lSeqs = lSeqs + x;
                }
                n++;
                for (int j = i + 1; j < l; j++) {
                    if (source[j] == 10 || source[j] == 13) {
                        name_seq[n] = new String(source, i + 1, j - i).trim();
                        i = j;
                        t = j + 1;
                        break;
                    }
                }
            }
        }
        int x = 0;
        for (int j = t; j < l; j++) {
            byte cc = tables.cdn[source[j]];
            if (cc > 0) {
                buf[x++] = cc;
            }
        }
        lSeqs = lSeqs + x;
        sequence[n] = new String(buf, 0, x, java.nio.charset.StandardCharsets.ISO_8859_1);
    }

    public void SetFolder(String p) {
        if (!p.isEmpty()) {
            fld = p + "__";
        }
    }

    public int getLength() {
        return (int) lSeqs;
    }
    private String fld = "";
    private String[] name_seq;
    private String[] sequence;
    private byte[] source = null;
    private double lSeqs = 0;
    private int ns = 0;
}
