
import java.util.Arrays;
import java.util.HashMap;

public final class SequencesClustering {

    // Constants for limits and thresholds
    private static final int MIN_SIMILARITY = 60;
    private static final int MAX_SIMILARITY = 80;
    private static final int KMER_LENGTH = 4;
    private static final int SPD_THRESHOLD = 21;
    private static final double DIFF_THRESHOLD = 1.35d;

    private int ncl;    // Number of clusters
    private int[] cx;   // Cluster mapping for sequences
    private int[][] d;  // Sequence ranges (start and length)

    private static final String[] KMERS = {"aatt", "acgt", "agct", "ttaa", "tgca", "tcga", "ccgg", "catg", "ctag", "ggcc", "gatc", "gtac"};

    /*
    private static final String[] KMERS = {
        "aatc", "aatg", "aact", "aacg", "aagt", "aagc", "atac", "atag", "attc", "attg",
        "atca", "atct", "atcc", "atcg", "atga", "atgt", "atgc", "atgg", "acat", "acag",
        "acta", "actt", "actc", "actg", "acct", "accg", "acga", "acgc", "acgg", "agat",
        "agac", "agta", "agtt", "agtc", "agtg", "agca", "agcc", "agcg", "aggt", "aggc",
        "taac", "taag", "tatc", "tatg", "taca", "tact", "tacc", "tacg", "taga", "tagt",
        "tagc", "tagg", "ttac", "ttag", "ttca", "ttcg", "ttga", "ttgc", "tcaa", "tcat",
        "tcac", "tcag", "tcta", "tctg", "tcca", "tccg", "tcgt", "tcgc", "tcgg", "tgaa",
        "tgat", "tgac", "tgag", "tgta", "tgtc", "tgct", "tgcc", "tgcg", "tgga", "tggc",
        "caat", "caag", "cata", "catt", "catc", "cact", "cacg", "caga", "cagt", "cagc",
        "cagg", "ctaa", "ctat", "ctac", "ctta", "cttg", "ctca", "ctcg", "ctga", "ctgt",
        "ctgc", "ctgg", "ccat", "ccag", "ccta", "cctg", "ccga", "ccgt", "cgaa", "cgat",
        "cgac", "cgag", "cgta", "cgtt", "cgtc", "cgtg", "cgca", "cgct", "cgga", "cggt",
        "gaat", "gaac", "gata", "gatt", "gatg", "gaca", "gact", "gacc", "gacg", "gagt",
        "gagc", "gtaa", "gtat", "gtag", "gtta", "gttc", "gtca", "gtct", "gtcc", "gtcg",
        "gtga", "gtgc", "gcaa", "gcat", "gcac", "gcag", "gcta", "gctt", "gctc", "gctg",
        "gcca", "gcct", "gcga", "gcgt", "ggat", "ggac", "ggta", "ggtc", "ggca", "ggct"
    };
     */
    public SequencesClustering(String seq, int[] x, int similarity) {
        int nseq = x.length / 2;
        if (nseq < 1) {
            return;  // No sequences to cluster
        }

        similarity = adjustSimilarity(similarity);

        // Populate sequence ranges (start, length)
        d = initializeSequenceRanges(x, nseq);

        // Sort the sequences by length (longest first)
        Arrays.sort(d, (a, b) -> Integer.compare(b[1], a[1]));

        // Perform clustering
        ncl = performClustering(seq, nseq, similarity);
    }

    // Adjusts the similarity to ensure it's within defined bounds
    private int adjustSimilarity(int similarity) {
        if (similarity < MIN_SIMILARITY) {
            return MIN_SIMILARITY;
        } else if (similarity > MAX_SIMILARITY) {
            return MAX_SIMILARITY;
        }
        return similarity;
    }

    // Initialize the sequence ranges array (start, length)
    private int[][] initializeSequenceRanges(int[] x, int nseq) {
        int[][] sequenceRanges = new int[nseq][2];
        for (int j = 0; j < nseq; j++) {
            int startIndex = j * 2;
            sequenceRanges[j][0] = x[startIndex];           // Start
            sequenceRanges[j][1] = x[startIndex + 1] - x[startIndex]; // Length
        }
        return sequenceRanges;
    }

    // Performs the clustering based on sequence similarities
    private int performClustering(String seq, int nseq, int similarity) {
        HashMap<String, Integer> pt = initializeKmerMap();
        int nkmers = pt.size();
        int[][] m2 = computeKmerFrequencies(seq, nseq, pt, nkmers);

        // Initialize clusters
        cx = new int[nseq];
        int clusters = 0;

        for (int i = 0; i < nseq; i++) {
            if (cx[i] == 0) {
                clusters++;
                cx[i] = clusters;

                for (int j = i + 1; j < nseq; j++) {
                    if (cx[j] == 0 && areSequencesSimilar(m2, i, j, nkmers, similarity)) {
                        cx[j] = clusters;
                    }
                }
            }
        }

        return clusters;
    }

    // Method to initialize the HashMap using a loop
    private HashMap<String, Integer> initializeKmerMap() {
        HashMap<String, Integer> pt = new HashMap<>();
        for (int i = 0; i < KMERS.length; i++) {
            pt.put(KMERS[i], i);
        }
        return pt;
    }

    // Computes the k-mer frequencies for each sequence
    private int[][] computeKmerFrequencies(String seq, int nseq, HashMap<String, Integer> pt, int nkmers) {
        int[][] m2 = new int[nseq][nkmers];

        for (int j = 0; j < nseq; j++) {
            int sequenceStart = d[j][0];
            int sequenceLength = d[j][1];

            for (int i = 0; i < sequenceLength - KMER_LENGTH + 1; i++) {
                String s = seq.substring(sequenceStart + i, sequenceStart + i + KMER_LENGTH);

                if (pt.containsKey(s)) {
                    m2[j][pt.get(s)]++;
                }

                String complement = dna.ComplementDNA2(s);
                if (pt.containsKey(complement)) {
                    m2[j][pt.get(complement)]++;
                }
            }
        }

        return m2;
    }

    // Determines if two sequences are similar based on k-mer frequencies
    private boolean areSequencesSimilar(int[][] m2, int i, int j, int nkmers, int similarity) {
        int[] sharedKmers = new int[nkmers + 1];
        int sharedKmerCount = 0;

        // Find shared k-mers
        for (int k = 0; k < nkmers; k++) {
            if (m2[i][k] > 0 && m2[j][k] > 0) {
                sharedKmers[++sharedKmerCount] = k;
            }
        }

        // Calculate similarity based on shared k-mers
        int matchCount = 0;
        int maxMatches = 0;

        for (int k = 1; k <= sharedKmerCount; k++) {
            for (int y = k + 1; y <= sharedKmerCount; y++) {
                maxMatches++;

                int di = Math.min(m2[i][sharedKmers[k]], m2[i][sharedKmers[y]]);
                int dj = Math.min(m2[j][sharedKmers[k]], m2[j][sharedKmers[y]]);

                if (Math.abs(di - dj) <= (int) (Math.min(di, dj) * DIFF_THRESHOLD)) {
                    matchCount++;
                }

                if (maxMatches > SPD_THRESHOLD) {
                    break;
                }
            }
        }

        return matchCount > 1 && (100 * matchCount / maxMatches) > similarity;
    }

    // Accessors
    public int getNcl() {
        return ncl;
    }

    public int[] Result() {
        return cx;
    }

    public int[][] ResultArray() {
        return d;
    }
}
