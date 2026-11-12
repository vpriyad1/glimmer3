# Glimmer3 Benchmarking on pO157

How well does an ab initio gene finder recover a curated annotation? This runs Glimmer3 on the *Escherichia coli* O157:H7 str. Sakai plasmid pO157 ([AB011549.2](https://www.ncbi.nlm.nih.gov/nuccore/AB011549.2), 92,721 bp) and scores every prediction against the GenBank CDS annotation, boundary by boundary.

**Finding:** Glimmer3 locates 71 of 78 annotated genes (91%) and gets the stop codon right on every one of them, but picks a different start codon on 20. It also predicts 33 genes the annotation does not contain, and those are short and low-scoring.

---

## Results

| | Count | Share of 78 |
|:--|--:|--:|
| Reference CDS with complete coordinates | 78 | |
| Glimmer3 predictions | 108 | |
| **Exact match** (5′ and 3′ agree) | **51** | 65.4% |
| 3′ agrees, 5′ differs | 20 | 25.6% |
| 5′ agrees, 3′ differs | 0 | 0% |
| Both ends differ | 1 | 1.3% |
| Reference gene with no prediction | 6 | 7.7% |
| Predictions with no reference overlap | 33 | |

Three things stand out.

**Stop codons are easy, start codons are not.** All 71 located genes have the correct 3′ end, and no gene has the right start with the wrong stop. A stop codon is unambiguous; an ORF usually offers several in-frame ATG/GTG/TTG candidates. In 14 of the 20 start-site disagreements Glimmer3 chose a start upstream of the annotated one, making the gene longer.

**The extra predictions look like noise.** The 33 predictions that overlap no annotated gene have a median length of 156 bp and a median score of 2.82, against 636 bp and 7.49 for predictions that do overlap one. A length or score filter would remove most of them.

**The misses are small genes.** All six annotated genes with no prediction are under 400 bp.

---

## Method

```
AB011549.2.fasta ── long-orfs ── extract ── build-icm ── glimmer3 ──┐
                                                                     ├── compare ── summary.json
AB011549.2.gbk ───────────────── parse_gbk.py ──────────────────────┘             detailed.json
```

Glimmer3 is trained on the plasmid's own long ORFs (41 ORFs, 50,799 bp) and then run on the full sequence:

```bash
long-orfs -n -t 1.15 AB011549.2.fasta longorfs.coords
extract -t AB011549.2.fasta longorfs.coords > train.seq
build-icm -r glimmer.icm < train.seq
glimmer3 AB011549.2.fasta glimmer.icm AB011549.2_glimmer
```

Each reference CDS is paired with the same-strand prediction that overlaps it by the most bases, and the pair is classified by its 5′ (start codon) and 3′ (stop codon) coordinates. Coordinates are compared biologically: on the minus strand the 5′ end is the higher coordinate.

| Script | Does |
|:--|:--|
| `scripts/parse_gbk.py` | Extracts CDS features from the GenBank file, including the `join()` gene that spans the plasmid origin |
| `scripts/parse_glimmer.py` | Parses the `.predict` file, including the prediction that wraps around the origin |
| `scripts/compare_ref_vs_glimmer.py` | Pairs, classifies, and writes `summary.json` and `detailed.json` |

Python 3 standard library only.

## Reproducing

```bash
conda create -n glimmer -c conda-forge -c bioconda glimmer -y
conda activate glimmer
./run_pipeline.sh
```

The script downloads the record from NCBI, runs Glimmer3, and prints the summary table.

## Limitations

- The annotation holds 83 CDS features. Five have partial coordinates (`<` or `>`), so they have no annotated start or stop to score and are excluded.
- A GenBank annotation is a reference, not ground truth. Some of the 33 unmatched predictions could be real unannotated genes, and some annotated start sites could be wrong.
- Glimmer3 was trained on 41 ORFs from a 92 kb plasmid. That is a small training set, and plasmids carry horizontally acquired genes whose composition differs from the host's.
- One replicon and one gene finder. A comparison against Prodigal on the same sequence would show how much of the start-site error is specific to Glimmer3.

---

Originally a midterm project for Advanced Practical Computer Concepts for Bioinformatics, M.S. Bioinformatics, Johns Hopkins University. The analysis scripts were rewritten and the pipeline rerun in November 2026.
