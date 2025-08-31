## 🚀 Quick-Start Guide

### 1. Install Java

* Requires **Java 24 or higher**
* [Download Java here](https://www.oracle.com/java/technologies/downloads/)
* Make sure Java is in your system PATH ([instructions](https://www.java.com/en/download/help/path.html))

---

### 2. Get Repeater2

* Copy `Repeater2.jar` from the `dist` folder to your working directory.

---

### 3. Run a Test Example

```bash
java -jar Repeater2.jar test/4.txt
```

---

### 4. Analyze a Genome Folder

```bash
java -Xms16g -Xmx64g -jar Repeater2.jar D:/Genomes/GRCh38.p14/
```

---

### 5. Useful Options

* `kmer=19` → set k-mer length (default 18)
* `sln=90` → minimal repeat string length (default 90)
* `-seqshow` → extract repeat sequences
* `-ssronly` → analyze only SSR/telomeric loci

Example:

```bash
java -jar Repeater2.jar test/ kmer=18 sln=90 -seqshow
```

---

### 6. Output

* Results saved in **GFF3 format** (tab-delimited text)
* Includes repeat clusters, positions, and optional extracted sequences

---

⚡ In just three steps (install Java → copy `.jar` → run `java -jar …`), you can analyze DNA sequences or entire genomes.

---

Would you like me to make this **the first section** of your README (so users see it right away), or keep it as a separate “Quick-Start” section after the introduction?
