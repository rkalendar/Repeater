import java.io.IOException;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

public class Repeater {

    public static void main(String[] args) {
        if (args.length > 0) {
            String infile = args[0]; // file path or Folder
            int kmer = 18;
            int minlen = 30;
            int seqlen = 60;
            int gap = kmer;
            int width = 0;
            int hight = 0;
            int flanksshow = 0;
            boolean seqshow = false;
            boolean ssrrun = false;
            boolean maskonly = false;
            int imgformat = PatternRepeatsSearching.FMT_BOTH;
            // By default results are saved to the current working directory.
            String outdir = System.getProperty("user.dir");

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            // Parse the output folder from the raw arguments (so the path keeps
            // its original case) and assemble the option string from the remaining
            // option arguments only. The input path (args[0]) and the out=/outdir=
            // folder value are deliberately excluded: a file path may legitimately
            // contain '=' (e.g. a folder named "kmer=20") which would otherwise be
            // mis-parsed as an option and silently corrupt the analysis.
            StringBuilder optbuf = new StringBuilder();
            for (int a = 1; a < args.length; a++) {
                String la = args[a].toLowerCase();
                if (la.startsWith("out=")) {
                    outdir = args[a].substring(4);
                } else if (la.startsWith("outdir=")) {
                    outdir = args[a].substring(7);
                } else {
                    optbuf.append(la).append(" ");
                }
            }
            String s = optbuf.toString();

            // Image format: format=both|png|svg|none (default both).
            if (s.contains("format=")) {
                int j = s.indexOf("format=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    switch (s.substring(j + 7, x).trim()) {
                        case "png" ->
                            imgformat = PatternRepeatsSearching.FMT_PNG;
                        case "svg" ->
                            imgformat = PatternRepeatsSearching.FMT_SVG;
                        case "none", "no", "off" ->
                            imgformat = PatternRepeatsSearching.FMT_NONE;
                        default ->
                            imgformat = PatternRepeatsSearching.FMT_BOTH;
                    }
                }
            }

            if (s.contains("ssronly")) {
                ssrrun = true;
            }
            if (s.contains("maskonly")) {
                maskonly = true;
            }
            if (s.contains("seqshow")) {
                seqshow = true;
            }
            if (s.contains("flanks=")) {
                flanksshow = 0;
                int j = s.indexOf("flanks=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    flanksshow = StrToInt(s.substring(j + 7, x));
                }
                if (flanksshow < 0) {
                    flanksshow = 0;
                }
                if (flanksshow > 1000) {
                    flanksshow = 1000;
                }
            }

            if (s.contains("kmer=")) {
                int j = s.indexOf("kmer=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    kmer = StrToInt(s.substring(j + 5, x));
                }
                if (kmer < 12) {
                    kmer = 12;
                }
                if (kmer > 32) {
                    // 2-bit packing into a 64-bit key supports at most 32 bases.
                    System.out.println("Note: kmer capped at 32 (maximum supported).");
                    kmer = 32;
                }
            }
            if (s.contains("min=")) {
                int j = s.indexOf("min=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    minlen = StrToInt(s.substring(j + 4, x));
                }
                if (minlen < kmer) {
                    minlen = kmer;
                }
            }
            if (s.contains("sln=")) {
                int j = s.indexOf("sln=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    seqlen = StrToInt(s.substring(j + 4, x));
                }
                if (seqlen < minlen) {
                    seqlen = minlen;
                }
            }
            if (s.contains("gap=")) {
                int j = s.indexOf("gap=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    gap = StrToInt(s.substring(j + 4, x));
                }
            }
            if (s.contains("image=")) { // image=40000x5000
                int j = s.indexOf("image=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    String[] d = s.substring(j + 6, x).split("x");
                    if (d.length > 1) {
                        width = StrToInt(d[0]);
                        hight = StrToInt(d[1]);
                        if (width < 1000) {
                            width = 1000;
                        }
                        if (width > 40000) {
                            width = 40000;
                        }
                        if (hight < 100) {
                            hight = 100;
                        }
                        if (hight == 0 || width == 0) {
                            width = 0;
                            hight = 0;
                        }
                    }
                }
            }

            // Make sure the output directory exists and really is a directory.
            if (outdir != null && !outdir.isEmpty()) {
                File od = new File(outdir);
                if (od.exists()) {
                    if (!od.isDirectory()) {
                        System.err.println("Output path is not a directory: " + outdir + " — using current directory instead.");
                        outdir = System.getProperty("user.dir");
                    }
                } else if (!od.mkdirs()) {
                    System.err.println("Could not create output directory: " + outdir + " — using current directory instead.");
                    outdir = System.getProperty("user.dir");
                }
            }
            System.out.println("Output directory: " + outdir);

            File folder = new File(infile);
            if (folder.exists() && (folder.isDirectory() || folder.isFile())) {
                if (folder.isDirectory()) {
                    File[] files = folder.listFiles();
                    int k = -1;
                    String[] filelist = new String[files.length];
                    for (File file : files) {
                        if (file.isFile()) {
                            filelist[++k] = file.getAbsolutePath();
                        }
                    }
                    for (String nfile : filelist) {
                        if (nfile != null) {
                            try {
                                SaveResult(nfile, kmer, minlen, seqlen, gap, flanksshow, seqshow, ssrrun, width, hight, maskonly, outdir, imgformat);
                            } catch (Exception e) {
                                System.err.println("Failed to open file: " + nfile);
                            }
                        }
                    }
                } else {
                    SaveResult(infile, kmer, minlen, seqlen, gap, flanksshow, seqshow, ssrrun, width, hight, maskonly, outdir, imgformat);
                }
            }

        } else {
            System.out.println("REPEATER2 (2024-2026) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/Repeater\n");
            System.out.println("Basic usage:");
            System.out.println("java -jar \\Repeater2\\dist\\Repeater2.jar <inputfile>/<inputfolderpath> <optional_commands>");
            System.out.println("Common options:");
            System.out.println("kmer=19\tminimal kmer=12 (default 18)");
            System.out.println("min=30\tinitial repeat block length (default min=30), it can be equal to 'kmer'");
            System.out.println("sln=90\trepeat block length (default sln=60), it can be equal to 'kmer'");
            System.out.println("flangs=100\textend the flanks of the repeat with an appropriate length (100 nt) (default flangs=0)");
            System.out.println("image=10000x300\t (by default, the dimensionality of the image is automatically determined)");
            System.out.println("format=both\timage output format: both|png|svg|none (default both)");
            System.out.println("out=<path>\toutput folder for results (default: current directory)");
            System.out.println("-seqshow\textract repeat sequences");
            System.out.println("-ssronly\tanalyzing only the SSR/telomers loci\n");
            System.out.println("-maskonly\tonly generate masked output (skip clustering)\n");
            System.out.println("java -jar \\Repeater2\\dist\\Repeater2.jar <inputfile> -ssronly -seqshow flanks=100");
            System.out.println("java -jar \\Repeater2\\dist\\Repeater2.jar <inputfile> kmer=18 min=30 sln=100 -nomask\n");
            System.out.println("java -jar \\Repeater2\\dist\\Repeater2.jar <inputfile> ssr=true seqshow=true flanks=100");
            System.out.println("Large genome settings:");
            System.out.println("java -jar -Xms16g -Xmx64g \\Repeater2\\dist\\Repeater2.jar <inputfile> kmer=20 sln=100 image=10000x300\n");
            System.out.println("Analysing all files in the directory:");
            System.out.println("java -jar \\Repeater2\\dist\\Repeater2.jar \\Repeater2\\test\\\n");
        }
    }

    public static int StrToInt(String str) {
        StringBuilder r = new StringBuilder();
        int z = 0;
        r.append(0);
        for (int i = 0; i < str.length(); i++) {
            char chr = str.charAt(i);
            if (chr > 47 && chr < 58) {
                r.append(chr);
                z++;
                if (z > 10) {
                    break;
                }
            }
            if (chr == '.' || chr == ',') {
                break;
            }
        }
        // Parse as long and clamp so an out-of-range value (e.g. a 10+ digit
        // number) cannot throw NumberFormatException and abort the program.
        long v;
        try {
            v = Long.parseLong(r.toString());
        } catch (NumberFormatException e) {
            v = Integer.MAX_VALUE;
        }
        return (int) Math.min(v, Integer.MAX_VALUE);
    }

    private static void SaveResult(String infile, int kmer, int minlen, int seqlen, int gap, int flanksshow, boolean seqshow, boolean ssrrun, int width, int hight, boolean maskonly, String outdir, int imgformat) {
        try {
            long startTime = System.nanoTime();
            byte[] binaryArray = Files.readAllBytes(Paths.get(infile));
            ReadingSequencesFiles rf = new ReadingSequencesFiles(binaryArray);
            if (rf.getNseq() == 0) {
                System.out.println("There is no sequence(s).");
                System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
                System.out.println(">seq1\nactacatactacatcactctctctccgcacag\n");
                return;
            }
            System.out.println("Running...");
            System.out.println("kmer=" + kmer);
            System.out.println("Minimal repeat length =" + minlen);
            System.out.println("Repeat block length =" + seqlen);
            System.out.println("Target file: " + infile);
            if (rf.getNseq() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getNseq());
            }
            System.out.println("Shown repeated sequence is " + seqshow);
            if (flanksshow > 0) {
                System.out.println("Flanks around sequence is " + flanksshow);
            }
            PatternRepeatsSearching s2 = new PatternRepeatsSearching();
            s2.SetSequences(rf.getSequences(), rf.getNames());
            s2.SetRepeatLen(kmer, minlen, seqlen);
            s2.SetShowSeq(seqshow);
            s2.SetFlanks(flanksshow);
            s2.SetFileName(infile);
            s2.SetOutputDir(outdir);
            s2.SetImageFormat(imgformat);
            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            if (ssrrun) {
                System.out.println("SSR analysis is running.\n");
                s2.SetTelomers(true);
                s2.RunSSR();
            } else {
                if (maskonly) {
                    s2.RunMask();
                } else {
                    s2.Run();
                }
            }
            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds\n");
        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }
}
