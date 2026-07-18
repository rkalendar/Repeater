import java.util.ArrayList;
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
                while (i < n && mask[i] == 1) {
                    i++;
                }
                int end = i - 1;              // last index of this run
                int highCount = end - start + 1;

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

// - No String substrings: 2-bit k-mer keys (A/C/G/T -> 0..3), skip windows with code 4 (N/gap).
// - Primitive open-addressing map (KmerMap): key=packed 2-bit k-mer,
//   cnt = occurrence count (negated while "marked"), pos = rolling index into u[].
//   This avoids boxing every rolling key to a Long and allocating an int[] per
//   distinct k-mer, which dominated cost/GC on large genomes.
// - Keeps the reverse-complement layout b[0..l-1], b[l]=4, b[l+1..l+l].
    private byte[] MaskingL(String seq, int kmer, int minlenblock) {
        if (seq == null || kmer <= 0 || minlenblock <= 0) {
            return new byte[0];
        }
        // A 64-bit 2-bit-packed key can represent at most 32 bases; beyond that
        // the rolling key would silently collide. Cap to the representable limit.
        if (kmer > 32) {
            kmer = 32;
        }
        final int l = seq.length();
        if (l == 0 || l < kmer) {
            return new byte[l];
        }

        final int g = l + l + 1;
        final byte[] msk = new byte[g];

        // Build code array: forward [0..l-1], sentinel [l]=4, reverse-complement [l+1..l+l]
        // Sequence is already filtered to ASCII base codes upstream, so
        // ISO-8859-1 is exact and skips the UTF-8 encoder pass on large inputs.
        byte[] b = seq.getBytes(java.nio.charset.StandardCharsets.ISO_8859_1);
        for (int i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];  // 0..4
        }
        b = java.util.Arrays.copyOf(b, g);
        b[l] = 4;
        for (int i = 1; i <= l; i++) {
            b[l + i] = tables.cdnat2[b[l - i]];
        }

        // rolling key helpers
        final long keyMask = (kmer >= 32) ? -1L : ((1L << (2 * kmer)) - 1L);

        // Count all valid k-mers on both strands. Distinct k-mers <= 2*(l-kmer+1);
        // the map grows on demand, so a modest initial capacity is fine.
        final KmerMap km = new KmerMap(Math.min((long) l, 1L << 20));
        int valid = 0;
        long key = 0;

        // Track "no-gaps" windows like the original (used elsewhere via NoGapsLength()).
        nogapslen = kmer - 1;

        // Forward
        for (int i = 0; i < l; i++) {
            int c = b[i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    km.count(key);
                    nogapslen++;
                }
            } else {
                valid = 0;
                key = 0;
            }
        }
        // Reverse
        valid = 0;
        key = 0;
        for (int i = 0; i < l; i++) {
            int c = b[l + 1 + i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    km.count(key);
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
                    int s = km.find(key);
                    if (s >= 0 && km.cnt[s] > 1) {
                        distinct++;
                        total += km.cnt[s];
                        km.cnt[s] = -km.cnt[s];
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
                    int s = km.find(key);
                    if (s >= 0 && km.cnt[s] > 1) {
                        distinct++;
                        total += km.cnt[s];
                        km.cnt[s] = -km.cnt[s];
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
                    int s = km.find(key);
                    if (s >= 0) {
                        if (km.cnt[s] < 0) {                  // first time we open this k-mer bucket
                            km.cnt[s] = -km.cnt[s];
                            if (km.cnt[s] > 1) {
                                km.pos[s] = z;
                                u[z] = pos;
                                t++;
                                u[0] = t;
                                x1[t][0] = km.cnt[s];
                                x1[t][1] = z;
                                z += km.cnt[s];
                            }
                        } else if (km.cnt[s] > 1) {           // subsequent occurrences
                            km.pos[s]++;
                            u[km.pos[s]] = pos;
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
                    int pos = l + i + 2 - kmer;   // indexing for the RC slice
                    int s = km.find(key);
                    if (s >= 0) {
                        if (km.cnt[s] < 0) {
                            km.cnt[s] = -km.cnt[s];
                            if (km.cnt[s] > 1) {
                                km.pos[s] = z;
                                u[z] = pos;
                                t++;
                                u[0] = t;
                                x1[t][0] = km.cnt[s];
                                x1[t][1] = z;
                                z += km.cnt[s];
                            }
                        } else if (km.cnt[s] > 1) {
                            km.pos[s]++;
                            u[km.pos[s]] = pos;
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

        // Pairwise extension per bucket (original logic preserved)
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
                                msk[l - y + l - r] = 1; // RC mapping
                            }
                        }
                    }
                }
            }
        }

        int r = 0;           // current gap length (zeros after a block)
        boolean inBlock = false;

        for (int i = 0; i < l; i++) {
            if (msk[i] > 0) {
                if (!inBlock) {           // new block
                    inBlock = true;
                    r = 0;
                } else {                  // continuing after a gap or same block
                    if (r > 0) {          // we just crossed a gap
                        if (r <= minlenblock) {
                            // fill only the gap: [i - r, i)
                            for (int j = i - r; j < i; j++) {
                                msk[j] = 1;
                            }
                        }
                        r = 0;
                    }
                }
            } else {
                if (inBlock) {            // count gap only after a block has started
                    r++;
                }
            }
        }
        return java.util.Arrays.copyOf(msk, l);
    }

    private static int extendBlock(byte[] b, int l, int g, int x, int y, int kmer) {
        int h = kmer, e = 0, p = 2;
        for (;;) {
            // bounds (mirror the original conditions)
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
