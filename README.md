# M. tuberculosis Variant Analysis — Lineage & Drug Resistance Profiling

> Adapted from the [Galaxy Training Network M. tuberculosis Variant Analysis tutorial](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/tb-variant-analysis/tutorial.html), converted from the graphical Galaxy workflow to an equivalent command-line pipeline.

Analysis of paired-end whole-genome sequencing data from an *M. tuberculosis*
clinical isolate, to call variants against the TB ancestral reference,
filter out unreliable regions, determine the phylogenetic lineage, and
screen for drug-resistance mutations.

Sample used: **004-2** (paired-end, from the GTN tutorial Zenodo dataset).

---

## Pipeline

```
FASTQ (paired-end)
  → FastQC (QC)
  → fastp (default trimming)
  → Snippy (bwa-mem + freebayes + SnpEff, combined)
  → tb_variant_filter (remove PE/PPE/PGRS repeat regions + low-depth sites)
  → TB-Profiler (lineage + drug resistance, called directly from the BAM)
  → sed + bcftools query + awk (clean annotation table of all filtered variants)
```

---

## 1. Reference & Data

```bash
mkdir -p tb_project/{fastq,ref,results}
cd tb_project

wget "https://zenodo.org/record/3960720/files/004-2_1.fastq.gz" -O fastq/004-2_R1.fastq.gz
wget "https://zenodo.org/record/3960720/files/004-2_2.fastq.gz" -O fastq/004-2_R2.fastq.gz
wget "https://zenodo.org/record/3960720/files/Mycobacterium_tuberculosis_ancestral_reference.gbk" -O ref/mtb_ancestral.gbk
```

**Why the "ancestral" reference?** It's a reconstructed common-ancestor
genome for all MTB lineages (Comas et al.), used to avoid reference bias
toward any single modern strain. The Genbank (`.gbk`) format — rather than
plain FASTA — is required so that Snippy's internal SnpEff step can
annotate variants directly against the gene models embedded in the file.

## 2. Quality Control

```bash
fastqc fastq/004-2_R1.fastq.gz fastq/004-2_R2.fastq.gz -o results/fastqc_before/
```

Result: **%GC = 66%**, matching the expected *M. tuberculosis* GC content
(~65.6%) — confirms the sample is genuinely MTB and not contaminated.
Quality scores were high throughout with the expected gradual drop toward
the 3' end of reads; no adapter contamination detected.

## 3. Trimming (fastp, default settings)

```bash
fastp -i fastq/004-2_R1.fastq.gz -I fastq/004-2_R2.fastq.gz \
      -o fastq/004-2_R1.trimmed.fastq.gz -O fastq/004-2_R2.trimmed.fastq.gz \
      --thread 4 \
      --html results/fastp_report.html --json results/fastp_report.json
```

Since FastQC showed no adapter contamination and generally high quality,
default fastp settings were used without additional aggressive trimming.

## 4. Variant Calling (Snippy)

```bash
snippy --cpus 4 \
       --outdir results/snippy_004-2 \
       --ref ref/mtb_ancestral.gbk \
       --R1 fastq/004-2_R1.trimmed.fastq.gz \
       --R2 fastq/004-2_R2.trimmed.fastq.gz \
       --mincov 10 \
       --minfrac 0.1
```

**Parameters:**
- `--mincov 10`: minimum read depth required to call a site.
- `--minfrac 0.1`: minimum fraction of reads supporting the alternate
  allele (lowered from the default 0.9) — deliberately permissive because
  TB samples can contain **mixed infections or heteroresistant
  subpopulations**, where a resistance mutation may be present in only a
  minority of the bacterial population but is still clinically relevant.

Snippy internally runs `bwa-mem` → `freebayes` → `SnpEff`, producing
`snps.vcf`, `snps.tab`, and `snps.bam` in one step.

**Result: 1,098 raw variants** (matching the GTN tutorial's expected count).

## 5. Variant Filtering (tb_variant_filter)

```bash
tb_variant_filter -R farhat_rlc --min_depth_filter --min_depth 10 \
    results/snippy_004-2/snps.vcf results/004-2_filtered.vcf
```

**Why filter?** *M. tuberculosis* contains highly repetitive gene families
(**PE/PPE/PGRS**) and insertion sequences that cause mapping errors and
generate false-positive variants. `tb_variant_filter`'s `farhat_rlc` list
("Refined Low Confidence" regions, Marin et al.) masks these, combined with
a minimum read-depth filter consistent with the `--mincov 10` used in
Snippy.

**Result: 130 variants masked** (~12%), leaving ~968 high-confidence
variants — matching the GTN tutorial's expected reduction (~131).

## 6. Lineage & Drug Resistance (TB-Profiler)

```bash
tb-profiler profile \
    --bam results/snippy_004-2/snps.bam \
    --prefix 004-2 \
    --dir results/tbprofiler_004-2 \
    --csv --txt
```

TB-Profiler calls variants independently, directly from the BAM, against
a curated drug-resistance mutation database (London School of Hygiene &
Tropical Medicine).

### Result

```
Lineage:         lineage4.4.1.1.1  (Euro-American family)
Drug-resistance: MDR-TB
Median depth:    188x
```

| Drug | Resistance | Gene | Mutation | Confidence |
|---|---|---|---|---|
| Rifampicin / Rifapentine | R | `rpoB` | p.Ser450Leu | Assoc w R (0.99) |
| Isoniazid | R | `katG` | p.Ser315Thr | Assoc w R (1.00) |
| Ethambutol | R | `embB` | p.Met300Val | Assoc w R (1.00) |
| Pyrazinamide | R | `pncA` | p.Tyr34Asp | Assoc w R (1.00) |
| Streptomycin | R | `rpsL` | p.Lys43Arg | Assoc w R (0.99) |
| Fluoroquinolones, injectables (2nd line) | Susceptible | — | none found | — |

**Interpretation:** All five mutations are well-documented, high-confidence
resistance markers (e.g. `rpoB` S450Leu is the single most common
rifampicin-resistance mutation worldwide, located in the RRDR hotspot
region). The sample is resistant to all first-line drugs tested but remains
**susceptible to fluoroquinolones and second-line injectables** — meaning
it qualifies as MDR-TB (and is resistant to pyrazinamide/streptomycin on
top of that) but is **not** XDR-TB, so second-line treatment options are
still available.

## 7. Annotated Variant Table

```bash
# Snippy's SnpEff step (using the Genbank reference) prefixes some
# auto-generated gene IDs with "GENE_" when no explicit /gene qualifier
# exists in the Genbank file — strip it for readability:
sed 's/GENE_//g' results/004-2_filtered.vcf > results/004-2_filtered.clean.vcf

echo -e "CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT" > results/004-2_annotated_summary.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' results/004-2_filtered.clean.vcf | \
awk -F'\t' '{split($5,a,"|"); print $1"\t"$2"\t"$3"\t"$4"\t"a[4]"\t"a[2]"\t"a[3]}' >> results/004-2_annotated_summary.txt
```

~40 `HIGH` impact variants (frameshift/stop-gained) were found across the
genome. **Caveat**: because the ancestral reference represents the deep
common ancestor of all MTB lineages, many apparent frameshifts are
single-base indels in homopolymer tracts that reflect either normal
lineage-specific divergence accumulated since the ancestor, or known
indel-calling artifacts in repetitive runs — not necessarily true
loss-of-function events in this specific isolate. They are listed as
candidates, not confirmed findings. The drug-resistance mutations
identified by TB-Profiler, by contrast, are curated, high-confidence SNPs
independent of this list.

---

## Files

- `results/004-2_annotated_summary.txt` — full annotated variant table (all filtered variants)
- `results/tbprofiler_004-2/results/004-2.txt` — TB-Profiler lineage & resistance report
- `results/004-2_filtered.vcf` — filtered variant calls
- `results/fastp_report.html` — QC/trimming report

---

## Pipeline Summary

```
FASTQ (004-2)
  → FastQC → fastp (default trimming)
  → Snippy (bwa-mem + freebayes + SnpEff vs. ancestral reference)
  → tb_variant_filter (farhat_rlc region mask + min depth 10)
  → TB-Profiler (lineage4.4.1.1.1, MDR-TB: RIF/INH/EMB/PZA/SM resistant)
  → sed + bcftools query (clean annotated variant table)
```
