# SIREN-neonatal-ecoli-genomics
Genomic analysis pipeline for neonatal E. coli isolates (SIREN study): typing, AMR, and core-genome phylogeny
# SIREN neonatal *E. coli* — genomic analysis pipeline

Commands, parameters and software versions used for the genomic characterisation and
core-genome phylogeny of neonatal *E. coli* isolates in the SIREN study.

All analyses were run on Penn State's Roar Collab HPC.

---

## 1. Isolate selection and genome assemblies

Genomes were included in downstream analyses only if they

1. met the GenomeTrakr quality thresholds for *E. coli* in MicroRunQC — mean read Q score
   ≥ 30, mean coverage ≥ 40×, assembly length 4.5–5.9 Mbp, ≤ 500 contigs — and
2. had a genome assembly available from NCBI.

Ten isolates failed one or both criteria and were excluded from all genomic analyses
(`data/excluded_isolates.tsv`).

Assemblies were retrieved from NCBI rather than re-assembled locally, so that every
downstream result is tied to a publicly available genome.

**Software:** NCBI Datasets CLI v18.38.0

```bash
# BioSample -> assembly accession, then download + unpack, one FASTA per BioSample
./scripts/download_assemblies.sh      # input: biosamples_clean.txt (one SAMN per line)
```

The script writes `fasta/<BioSample>.fna` for each genome and `samn_to_assembly.tsv`
mapping each BioSample to the assembly accession it was downloaded from. GenBank (GCA)
accessions were preferred over RefSeq (GCF); where several versions existed, the highest
was used.

Accessions: `data/isolate_inventory.tsv` (isolate ID → BioSample → assembly accession).

---

## 2. Genomic characterisation

| Analysis | Tool | Version | Database |
|---|---|---|---|
| Sequence type (Achtman scheme) | `mlst` | 2.35.0 | `ecoli_achtman_4` (PubMLST) |
| Phylogroup | ClermonTyping | 24.02 | bundled primers + mash reference |
| O/H serotype | ECTyper | 2.0.0 | bundled (run on GalaxyTrakr) |
| Capsule (K) type | kTYPr | 1.0 (commit `925d2ce`) | bundled HMM profiles |
| Acquired AMR genes | ABRicate | 1.4.0 | NCBI AMRFinderPlus, 8,453 seqs, 2026-11-02 |

```bash
# Sequence types
mlst --scheme ecoli_achtman_4 --threads 8 fasta/*.fna > mlst_results.tsv

# Phylogroups
clermonTyping.sh --fastafile fasta_list.txt --name upenn_all --minimal

# Capsule types
ktypr -i fasta -o ktypr_results -n 8 -v

# Acquired AMR genes
abricate --db ncbi --fofn fasta_list.txt --minid 90 --mincov 80 --threads 8 --nopath \
    > abricate_ncbi.tab
abricate --summary abricate_ncbi.tab > abricate_ncbi_summary.tab
```

Notes on these settings:

- `--minid 90 --mincov 80` — ABRicate's defaults are 80/80. A gene detected on more than
  one contig is counted once.
- ABRicate assigns the closest reference allele by nucleotide identity, so allele
  designations below 100 % identity are approximate.
- ABRicate detects acquired genes only; chromosomal point mutations (*gyrA*, *parC*) were
  not assessed.
- ClermonTyping `--minimal` writes the results table only. The tool's mash-based
  prediction was not used to assign phylogroups.
- kTYPr was run in its default flanking-region mode (30 kb window around *kpsC*).
  Isolates with no *kpsC* detected are reported as having no identifiable group 2/3
  capsule locus rather than as capsule-negative.

---

## 3. Core-genome phylogeny

**Bakta** (annotation) → **Panaroo** (pangenome + core-gene alignment) →
**snp-sites** + **IQ-TREE** (maximum-likelihood tree).

252 genomes: all isolates in the invasive, urine, respiratory and CSF cohorts that met
the inclusion criteria above.

### 3a. Annotation — Bakta 1.12.2

Database: Bakta DB v6.0, `light` type
(DOI [10.5281/zenodo.14916843](https://doi.org/10.5281/zenodo.14916843))

```bash
bakta --db <path>/db-light \
      --output bakta/<BioSample> --prefix <BioSample> --locus-tag <BioSample> \
      --genus Escherichia --species coli \
      --threads 8 --skip-plot --force \
      fasta/<BioSample>.fna
```

Run across all genomes by `scripts/01_bakta.sbatch` (4 genomes in parallel, 8 threads
each); the resulting `<BioSample>.gff3` files are collected into `gff/`. Locus tags are
set to the BioSample accession, so every gene in the pangenome traces back to its genome.

### 3b. Pangenome and core alignment — Panaroo 1.8.1

```bash
panaroo -i gff/*.gff3 -o panaroo_out \
        --clean-mode strict --remove-invalid-genes \
        -a core --aligner mafft --core_threshold 0.98 \
        -t 32
```

| Flag | Value | Effect |
|---|---|---|
| `--clean-mode` | `strict` | Most aggressive filtering of likely annotation errors and contamination; removes gene clusters not supported by the pangenome graph. |
| `--remove-invalid-genes` | on | Drops genes with bad start/stop or internal stops instead of failing. |
| `-a core` | — | Aligns core genes only, not the full pangenome. |
| `--aligner` | `mafft` | Per-gene aligner. |
| `--core_threshold` | `0.98` | A gene is core if present in ≥ 98 % of the 252 genomes. |

3,368 core genes, of which 3,161 were retained after Panaroo's post-alignment entropy
filter (threshold set automatically to 0.0587). Two alignments are written:

- `core_gene_alignment.aln` — all 3,368 core genes (851 MB)
- `core_gene_alignment_filtered.aln` — the 3,161 retained genes (785 MB) ← used for the tree

`gene_presence_absence.csv` holds the accessory genome, one row per gene cluster.

Runtime: 58 min on 32 cores, peak RAM 15.6 GB.

### 3c. Tree — snp-sites 2.5.1 + IQ-TREE 3.1.4

```bash
# Reduce the core alignment to variable sites, and count the constant sites
snp-sites -o core_snps.aln panaroo_out/core_gene_alignment_filtered.aln
FCONST=$(snp-sites -C panaroo_out/core_gene_alignment_filtered.aln)

# Maximum-likelihood tree
iqtree3 -s core_snps.aln -fconst "$FCONST" -m MFP -B 1000 -T 32 --prefix tree/core_snps
```

| Flag | Value | Effect |
|---|---|---|
| `-fconst` | counts from `snp-sites -C` | Adds the invariant sites back as counts, so branch lengths are scaled to substitutions per site across the core genome rather than per SNP. |
| `-m MFP` | ModelFinder Plus | Substitution model selected by BIC. |
| `-B 1000` | — | 1,000 ultrafast bootstrap replicates. |
| `-T 32` | — | Threads. |

2,797,775 constant sites added to 265,074 variable sites (3,062,849 total).
Best-fit model **GTR+F+I+R7** (BIC). Consensus tree log-likelihood −9,062,826.67.
Runtime: 10 h 30 min on 32 cores, peak RAM 12.5 GB.

Two pairs of genomes are identical over the core alignment; both pairs are two isolates
from the same infant (HUP0018-1 / HUP0018-2 and UNEB0001-1 / UNEB0001-2). Both were
retained in the tree.

**Outputs**

| File | Contents |
|---|---|
| `core_snps.treefile` | Maximum-likelihood tree (Newick), tips labelled by BioSample |
| `core_snps.contree` | Consensus tree with ultrafast bootstrap support on branches |
| `core_snps.iqtree` | Full report: model selection, parameters, likelihoods |
| `core_snps.splits.nex` | Split support values |
| `core_snps.mldist` | Pairwise maximum-likelihood distances |

---

## Metadata

`data/tree_metadata.tsv` — one row per tip, keyed by BioSample:

| Column | Notes |
|---|---|
| `BioSample` | matches tip labels in the tree |
| `Isolate` | lab isolate ID |
| `Patient` | study ID; an infant may contribute more than one isolate |
| `Center` | contributing site |
| `Cohort` | Invasive / Urine / Respiratory / CSF (2 isolates appear in both Invasive and CSF) |
| `Source` | specimen type |
| `ST` | sequence type as recorded in the study dataset |
| `ST_mlst` | sequence type from this run of `mlst` v2.35.0 |
| `Phylogroup` | ClermonTyping |
| `K_type` | kTYPr |
| `Serotype` | ECTyper (O:H) |
| `fimH` | FimTyper |

Isolate IDs are study codes; no identifiable patient information is included.

---

## Repository layout

```
scripts/
  download_assemblies.sh     BioSample -> NCBI assembly download
  00_setup.sh                conda environments + Bakta database
  01_bakta.sbatch            annotation (SLURM)
  02_panaroo.sbatch          pangenome + core alignment (SLURM)
  03_tree.sbatch             snp-sites + IQ-TREE (SLURM)
data/
  isolate_inventory.tsv      isolate ID, BioSample, assembly accession
  tree_metadata.tsv          per-tip metadata
  excluded_isolates.tsv      the 10 excluded isolates and why
results/
  core_snps.treefile         ML tree
  core_snps.contree          consensus tree with bootstrap support
  core_snps.iqtree           IQ-TREE report
```

Genome assemblies and the core-gene alignment are not stored here. Assemblies are public
at NCBI under the accessions in `data/isolate_inventory.tsv`; the alignment (785 MB) can
be regenerated from the commands above, or requested directly.

---

## References

- Schwengers O, *et al.* Bakta: rapid and standardized annotation of bacterial genomes
  via alignment-free sequence identification. *Microb Genom.* 2021;7(11):000685.
- Tonkin-Hill G, *et al.* Producing polished prokaryotic pangenomes with the Panaroo
  pipeline. *Genome Biol.* 2020;21:180.
- Page AJ, *et al.* SNP-sites: rapid efficient extraction of SNPs from multi-FASTA
  alignments. *Microb Genom.* 2016;2(4):e000056.
- Minh BQ, *et al.* IQ-TREE 2: new models and efficient methods for phylogenetic
  inference in the genomic era. *Mol Biol Evol.* 2020;37(5):1530–1534.
- Kalyaanamoorthy S, *et al.* ModelFinder: fast model selection for accurate phylogenetic
  estimates. *Nat Methods.* 2017;14:587–589.
- Hoang DT, *et al.* UFBoot2: improving the ultrafast bootstrap approximation.
  *Mol Biol Evol.* 2018;35(2):518–522.
- Beghain J, *et al.* ClermonTyping: an easy-to-use and accurate in silico method for
  *Escherichia* genus strain phylotyping. *Microb Genom.* 2018;4(7):e000192.
- Bessonov K, *et al.* ECTyper: in silico *Escherichia coli* serotype and species
  prediction. *Microb Genom.* 2021;7(12):000728.
- Miravet-Verde S, *et al.* In silico typing maps the natural diversity of *Escherichia
  coli* transporter-dependent capsules. *Nat Microbiol.* 2026;11(5):1217–1232.
- Feldgarden M, *et al.* AMRFinderPlus and the Reference Gene Catalog. *Sci Rep.*
  2021;11:12728.
- Jolley KA, *et al.* Open-access bacterial population genomics: BIGSdb software and the
  PubMLST.org website. *Wellcome Open Res.* 2018;3:124.
- O'Leary NA, *et al.* Exploring and retrieving sequence and metadata for species across
  the tree of life with NCBI Datasets. *Sci Data.* 2024;11:732.
- Seemann T. `mlst`, `abricate`. https://github.com/tseemann
