# FAQ

**RiboFlow** is a Nextflow DSL2 pipeline for ribosome-profiling and RNA-seq data, with STAR-based **genome alignment** alongside the bowtie2 **transcriptome → `.ribo`** path. Entry point: `main.nf`.

### Does it produce `.ribo` files?

Yes, with `transcriptome.run: true` (bowtie2 → `ribopy create` → `ribopy merge` into `all.ribo`). The genome path (`genome.run: true`) emits dedup BAM/BED, stats and optional bigWigs, but no `.ribo`. Run either or both.

### What outputs does it produce?

```
<out>/
  alignments/ribo/{individual,merged,stranded}/   # dedup BAM + BED
  bigwigs/ribo/                                   # if do_bigwig
  fastqc/                                         # if do_fastqc
  ribo/                                           # *.ribo + all.ribo (if transcriptome.run)
  stats/{genome,transcriptome}/{stats,individual_stats}.csv
  rnaseq/{alignments,bigwigs,fastqc,stats}/       # if do_rnaseq, same layout
```

`do_bigwig` and `do_fastqc` default to `false`. Intermediates are cached under `intermediates/` via `storeDir`.

### Does it support UMIs?

Yes: `dedup_method: "umicollapse"` plus `umi_tools_extract_arguments` for your UMI layout. See `examples/example_umi_uniq.yaml`.

### Does it support paired-end RNA-seq?

Yes, on the RNA-seq **genome** path — give `[R1, R2]` per sample (see `examples/example_rnaseq_pe.yaml`). Ribo-seq is single-end only. Rejected at startup: PE with `rnaseq.dedup_method: "position"` (use `umicollapse` or `none`), and PE on the RNA-seq transcriptome path.

### Can I run it without Docker?

Yes. `environment.yaml` is a single conda env with every tool (`umicollapse` from bioconda); use `-profile conda`. It is Linux-only — on macOS use the Docker/Apptainer image.

### How are stats generated?

Per-lane counts come from the cutadapt, bowtie2 and STAR logs plus the qpass/dedup count files, then combined by [`rfc`](https://github.com/ribosomeprofiling/RFCommands) into `stats/*/stats.csv` and `individual_stats.csv`.

### Where do I get references?

You provide them: `input.reference.filter` (bowtie2 rRNA index prefix), `input.reference.genome` (STAR index dir), and for `.ribo` also `transcriptome`, `regions`, `transcript_lengths`. [references_for_riboflow](https://github.com/ribosomeprofiling/references_for_riboflow) works here. For the STAR index either build it yourself with the pinned STAR version, or set `genome_fasta` + `gtf` and let the pipeline run `genomeGenerate` (see README).
