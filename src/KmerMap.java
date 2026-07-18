/**
 * Minimal open-addressing hash map from a packed 2-bit k-mer key ({@code long})
 * to a pair of ints: {@code cnt} (occurrence count, negated as a transient
 * "marked" flag) and {@code pos} (rolling write index into an occurrence array).
 *
 * <p>Purpose-built for the k-mer counting passes in
 * {@link MaskingPairwiseAlignmentSequence} and {@link PatternRepeatsSearching}:
 * it replaces {@code HashMap<Long,int[]>} to avoid boxing every rolling key into
 * a {@code Long} and allocating an {@code int[]} per distinct k-mer, which
 * dominated allocation/GC cost on multi-megabase genomes.</p>
 *
 * <p>An empty slot is encoded by {@code cnt == 0}: a live count is always 1,
 * {@code >= 2}, or {@code <= -2} (it starts at 1 and is only negated while
 * {@code > 1}), so it is never 0 and no separate occupancy array is needed.
 * Keys are only ever inserted via {@link #count}; the marking/collecting phases
 * use {@link #find} only, so the table is stable and slot indices stay valid
 * throughout those phases.</p>
 *
 * <p>Not thread-safe; each analysis creates its own instance.</p>
 */
final class KmerMap {

    long[] key;
    int[] cnt;   // occurrence count; 0 == empty slot
    int[] pos;   // rolling index into the occurrence array (the int[]'s 2nd field)

    private int mask;      // capacity - 1 (capacity is a power of two)
    private int size;
    private int threshold; // grow when size reaches this

    private static final double LOAD = 0.6;
    private static final int MAX_CAP = 1 << 30;

    KmerMap(long expected) {
        int cap = tableSizeFor((long) (expected / LOAD) + 1);
        key = new long[cap];
        cnt = new int[cap];
        pos = new int[cap];
        mask = cap - 1;
        threshold = thresholdFor(cap);
    }

    private static int tableSizeFor(long n) {
        long c = 4;
        while (c < n && c < MAX_CAP) {
            c <<= 1;
        }
        return (int) c;
    }

    private static int thresholdFor(int cap) {
        // At the terminal capacity we can no longer grow, so keep real headroom
        // (~1/8 of the table) and let count() fail fast on overflow instead of
        // filling to 100% and spinning forever. Below MAX_CAP, grow early to keep
        // probe chains short.
        return (cap >= MAX_CAP) ? (cap - (cap >>> 3)) : (int) (cap * LOAD);
    }

    // MurmurHash3 64-bit finalizer: spreads the 2*kmer significant bits.
    private static int spread(long h) {
        h ^= (h >>> 33);
        h *= 0xff51afd7ed558ccdL;
        h ^= (h >>> 33);
        return (int) h;
    }

    /** Record one occurrence: insert with cnt=1, or increment an existing count. */
    void count(long k) {
        int i = spread(k) & mask;
        while (cnt[i] != 0) {
            if (key[i] == k) {
                cnt[i]++;
                return;
            }
            i = (i + 1) & mask;
        }
        key[i] = k;
        cnt[i] = 1;
        pos[i] = 0;
        if (++size >= threshold) {
            grow();
        }
    }

    /** @return slot index for {@code k}, or -1 if the key is not present. */
    int find(long k) {
        int i = spread(k) & mask;
        while (cnt[i] != 0) {
            if (key[i] == k) {
                return i;
            }
            i = (i + 1) & mask;
        }
        return -1;
    }

    private void grow() {
        int oldCap = mask + 1;
        if (oldCap >= MAX_CAP) {
            // Cannot address more than MAX_CAP slots with int-indexed arrays.
            // Reaching here needs > ~0.9 billion distinct k-mers (a ~16 GB index),
            // which also overflows the downstream int position arrays, so such
            // inputs are unsupported by design. Fail fast with a clear error
            // rather than let the probe loop fill the table and hang. Because the
            // threshold leaves headroom, an empty slot always exists, so callers
            // never spin: the insert that tripped the threshold has completed and
            // this exception is thrown before the table can fill.
            throw new IllegalStateException(
                    "k-mer index exceeded maximum capacity (" + MAX_CAP
                    + " slots): too many distinct k-mers for this sequence/kmer size");
        }
        int newCap = oldCap << 1;
        long[] nk = new long[newCap];
        int[] nc = new int[newCap];
        int[] np = new int[newCap];
        int nmask = newCap - 1;
        for (int j = 0; j < oldCap; j++) {
            if (cnt[j] != 0) {
                int i = spread(key[j]) & nmask;
                while (nc[i] != 0) {
                    i = (i + 1) & nmask;
                }
                nk[i] = key[j];
                nc[i] = cnt[j];
                np[i] = pos[j];
            }
        }
        key = nk;
        cnt = nc;
        pos = np;
        mask = nmask;
        threshold = thresholdFor(newCap);
    }
}
