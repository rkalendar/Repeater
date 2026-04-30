import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class MaskingPairwiseAlignmentSequence {

    private byte[] mask;
    private long nogapslen = 0;
    private long repeatslen = 0;

    public int[] Mask(String seq, int kmer, int minlenseq) {
        mask = MaskingL(seq, kmer, minlenseq);
        int n = mask.length;

        List<int[]> blocks = new ArrayList<>();

        int i = 0;
        while (i < n) {
            if (mask[i] == 1) {
                int start = i;
                int end = i;
                int highCount = 0;

                while (i < n && mask[i] == 1) {
                    if (mask[i] == 1) {
                        highCount++;
                    }
                    end = i;
                    i++;
                }

                if (highCount > kmer) {
                    int seqStart = start;
                    int seqEnd = Math.min(n - 1, end + kmer - 1);

                    if (blocks.isEmpty()
                            || blocks.get(blocks.size() - 1)[1] < seqStart) {
                        blocks.add(new int[]{seqStart, seqEnd});
                    } else {
                        int[] prev = blocks.get(blocks.size() - 1);
                        prev[1] = Math.max(prev[1], seqEnd);
                    }
                }
            } else {
                i++;
            }
        }

        int totalLen = 0;
        List<int[]> filtered = new ArrayList<>();
        for (int[] block : blocks) {
            int start = block[0];
            int end = block[1];
            int length = end - start + 1;
            if (length > minlenseq) {
                filtered.add(new int[]{start, length});
                totalLen += length;
            }
        }
        this.repeatslen = totalLen;

        int[] result = new int[filtered.size() * 2];
        int idx = 0;
        for (int[] block : filtered) {
            result[idx++] = block[0];
            result[idx++] = block[1];
        }
        return result;
    }

    public byte[] ArrayTrim(byte[] srcArray, int n) {
        byte[] destArray = new byte[n];
        System.arraycopy(srcArray, 0, destArray, 0, n);
        return destArray;
    }

    public byte[] ArrayExtendByte(byte[] srcArray, int n) {
        byte[] destArray = new byte[srcArray.length + n];
        System.arraycopy(srcArray, 0, destArray, 0, srcArray.length);
        return destArray;
    }

// - No String substrings: 2-bit k-mer keys (A/C/G/T -> 0..3), skip windows with code 4 (N/gap).
// - HashMap<Long,int[]>: entry[0]=count (negated while "marked"), entry[1]=rolling index into u[].
// - Fixed (n == 0 || t == 0) and clarified loops.
// - Keeps your reverse-complement layout b[0..l-1], b[l]=4, b[l+1..l+l].
    private byte[] MaskingL(String seq, int kmer, int minlenblock) {
        if (seq == null || kmer <= 0 || minlenblock <= 0) {
            return new byte[0];
        }
        final int l = seq.length();
        if (l == 0 || l < kmer) {
            return new byte[l];
        }

        final int g = l + l + 1;
        final byte[] msk = new byte[g];

        // Build code array: forward [0..l-1], sentinel [l]=4, reverse-complement [l+1..l+l]
        byte[] b = seq.getBytes();           // expect ASCII ACGTN
        for (int i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];  // 0..4
        }
        b = java.util.Arrays.copyOf(b, g);
        b[l] = 4;
        for (int i = 1; i <= l; i++) {
            b[l + i] = tables.cdnat2[b[l - i]];
        }

        // rolling key helpers
        final long keyMask = (kmer == 32) ? -1L : ((1L << (2 * kmer)) - 1L);

        // Count all valid k-mers on both strands
        java.util.HashMap<Long, int[]> km = new java.util.HashMap<>(Math.min(l, 1 << 20));
        int valid = 0;
        long key = 0;

        // Track "no-gaps" windows like your original (if you use this elsewhere)
        nogapslen = kmer - 1;

        // Forward msk
        for (int i = 0; i < l; i++) {
            int c = b[i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    long k = key;
                    int[] v = km.get(k);
                    if (v == null) {
                        km.put(k, new int[]{1, 0});
                    } else {
                        v[0]++;
                    }
                    nogapslen++;
                }
            } else {
                valid = 0;
                key = 0;
            }
        }
        // Reverse msk
        valid = 0;
        key = 0;
        for (int i = 0; i < l; i++) {
            int c = b[l + 1 + i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    long k = key;
                    int[] v = km.get(k);
                    if (v == null) {
                        km.put(k, new int[]{1, 0});
                    } else {
                        v[0]++;
                    }
                }
            } else {
                valid = 0;
                key = 0;
            }
        }

        // Mark repeated k-mers and compute totals
        int distinct = 0, total = 0;
        valid = 0;
        key = 0;
        for (int i = 0; i < l; i++) {               // forward
            int c = b[i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    int[] v = km.get(key);
                    if (v != null && v[0] > 1) {
                        distinct++;
                        total += v[0];
                        v[0] = -v[0];
                    }
                }
            } else {
                valid = 0;
                key = 0;
            }
        }
        valid = 0;
        key = 0;
        for (int i = 0; i < l; i++) {               // reverse
            int c = b[l + 1 + i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    int[] v = km.get(key);
                    if (v != null && v[0] > 1) {
                        distinct++;
                        total += v[0];
                        v[0] = -v[0];
                    }
                }
            } else {
                valid = 0;
                key = 0;
            }
        }

        if (total == 0 || distinct == 0) {
            return java.util.Arrays.copyOf(msk, l);
        }

        // Occurrence buckets
        int[] u = new int[total + 1];       // positions; u[0] == t (# of distinct)
        int[][] x1 = new int[distinct + 1][2]; // [i][0]=count, [i][1]=start index in u
        int z = 1, t = 0;

        // Collect positions (forward)
        valid = 0;
        key = 0;
        for (int i = 0; i < l; i++) {
            int c = b[i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    int pos = i + 1 - kmer;
                    int[] v = km.get(key);
                    if (v != null) {
                        if (v[0] < 0) {                       // first time we open this k-mer bucket
                            v[0] = -v[0];
                            if (v[0] > 1) {
                                v[1] = z;
                                u[z] = pos;
                                t++;
                                u[0] = t;
                                x1[t][0] = v[0];
                                x1[t][1] = z;
                                z += v[0];
                            }
                        } else if (v[0] > 1) {                // subsequent occurrences
                            v[1]++;
                            u[v[1]] = pos;
                        }
                    }
                }
            } else {
                valid = 0;
                key = 0;
            }
        }
        // Collect positions (reverse)
        valid = 0;
        key = 0;
        for (int i = 0; i < l; i++) {
            int c = b[l + 1 + i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    int pos = l + i + 2 - kmer;   // your indexing for RC slice
                    int[] v = km.get(key);
                    if (v != null) {
                        if (v[0] < 0) {
                            v[0] = -v[0];
                            if (v[0] > 1) {
                                v[1] = z;
                                u[z] = pos;
                                t++;
                                u[0] = t;
                                x1[t][0] = v[0];
                                x1[t][1] = z;
                                z += v[0];
                            }
                        } else if (v[0] > 1) {
                            v[1]++;
                            u[v[1]] = pos;
                        }
                    }
                }
            } else {
                valid = 0;
                key = 0;
            }
        }

        if (t == 0) {
            return java.util.Arrays.copyOf(msk, l);
        }

        // Pairwise extension per bucket (kept your original logic)
        for (int bi = 1; bi <= t; bi++) {
            int count = x1[bi][0];
            int start = x1[bi][1];
            for (int a = 1; a <= count; a++) {
                int x = u[start + a - 1];
                for (int bpos = a + 1; bpos <= count; bpos++) {
                    int y = u[start + bpos - 1];
                    if (msk[x] != 0 || msk[y] != 0) {
                        continue;
                    }

                    int h = extendBlock(b, l, g, x, y, kmer);
                    if (h > minlenblock) {
                        for (int r = 0; r < h; r++) {
                            msk[x + r] = 1;
                            if (y < l) {
                                msk[y + r] = 1;
                            } else {
                                msk[l - y + l - r] = 1; // your RC mapping
                            }
                        }
                    }
                }
            }
        }

        int x = -1;          // start index of current block
        int y = -1;          // last index seen in the current/previous block
        int r = 0;           // current gap length (zeros after a block)
        boolean inBlock = false;

        for (int i = 0; i < l; i++) {
            if (msk[i] > 0) {
                if (!inBlock) {           // new block
                    x = i;
                    inBlock = true;
                    r = 0;
                } else {                  // continuing after a gap or same block
                    if (r > 0) {          // we just crossed a gap
                        if (r <= minlenblock) {
                            // fill only the gap: [i - r, i)
                            for (int j = i - r; j < i; j++) {
                                msk[j] = 1;
                            }
                        } else {
                            // long gap: start a fresh block at i
                            x = i;
                        }
                        r = 0;
                    }
                }
                y = i;
            } else {
                if (inBlock) {            // count gap only after a block has started
                    r++;
                }
            }
        }
        // return ArrayTrim(msk, l);
        return java.util.Arrays.copyOf(msk, l);
    }

    private static int extendBlock(byte[] b, int l, int g, int x, int y, int kmer) {
        int h = kmer, e = 0, p = 2;
        for (;;) {
            // bounds (mirror your conditions)
            if ((y > l && y + h > g - 1) || (y <= l && y + h > l - 1)) {
                break;
            }
            if ((x < l && x + h > l - 1) || (x >= l && x + h > g - 1)) {
                break;
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
                    h -= 3;
                    break;
                } // keep original behavior
                e++;
                p = 0;
            }
            h++;
        }
        return h;
    }

    private byte[] Masking(String seq, int kmer, int minlenblock) {
        int k, n, t, i, h, e, y, z, x, j, r, p;
        int l = seq.length();
        int g = l + l + 1;

        final byte[] msk = new byte[g];

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
        nogapslen = kmer - 1;
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                nogapslen++;
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
                if (p > 1) {
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
                    if (p > 1) {
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
            return java.util.Arrays.copyOf(msk, l);
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
                    if (px2.get(s)[0] > 1) {
                        px2.get(s)[1] = z;
                        u[z] = i + 1 - kmer;
                        t++;
                        u[0] = t;
                        x1[t][1] = z;
                        x1[t][0] = px2.get(s)[0];
                        z = z + px2.get(s)[0];
                    }
                } else {
                    if (p > 1) {
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
                        if (px2.get(s)[0] > 1) {
                            px2.get(s)[1] = z;
                            u[z] = l + i + 2 - kmer;
                            t++;
                            u[0] = t;
                            x1[t][1] = z;
                            x1[t][0] = px2.get(s)[0];
                            z = z + px2.get(s)[0];
                        }
                    } else {
                        if (p > 1) {
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
            return java.util.Arrays.copyOf(msk, l);
        }

        for (i = 1; i < t + 1; i++) {
            for (k = 1; k <= x1[i][0]; k++) {
                x = u[x1[i][1] + k - 1];

                for (r = k + 1; r <= x1[i][0]; r++) {
                    y = u[x1[i][1] + r - 1];

                    if (msk[x] == 0 && msk[y] == 0) {
                        h = kmer;
                        e = 0; // no mistmatches
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
                            for (r = 0; r < h; r++) {
                                msk[x + r] = 1;
                                if (y < l) {
                                    msk[y + r] = 1;
                                } else {
                                    msk[l - y + l - r] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        x = -1;          // start index of current block
        y = -1;          // last index seen in the current/previous block
        r = 0;           // current gap length (zeros after a block)
        boolean inBlock = false;

        for (i = 0; i < l; i++) {
            if (msk[i] > 0) {
                if (!inBlock) {           // new block
                    x = i;
                    inBlock = true;
                    r = 0;
                } else {                  // continuing after a gap or same block
                    if (r > 0) {          // we just crossed a gap
                        if (r <= minlenblock) {
                            // fill only the gap: [i - r, i)
                            for (j = i - r; j < i; j++) {
                                msk[j] = 1;
                            }
                        } else {
                            // long gap: start a fresh block at i
                            x = i;
                        }
                        r = 0;
                    }
                }
                y = i;
            } else {
                if (inBlock) {            // count gap only after a block has started
                    r++;
                }
            }
        }
        return java.util.Arrays.copyOf(msk, l);
    }

    public long NoGapsLength() {
        return nogapslen;
    }

    public byte[] getByteMask() {
        return mask;
    }

    public long RepeatLength() {
        return repeatslen;
    }
}
