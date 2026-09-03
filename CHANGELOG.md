# Changelog

## 2.0.0

Nextflow DSL2 rewrite of the DSL1 `RiboFlow.groovy`.

### Added

- STAR genome alignment path (ribo-seq): clip → rRNA filter → STAR → quality filter → dedup (`umicollapse` / `position` / `none`)-> merge -> strand-specific bigWigs.
- Genomic RNA-seq alignment can handle both single-end and paired-end inputs.
- UMI-aware deduplication with `umicollapse` (from bioconda; single- and paired-end), alongside `position` and `none`.
- Unique-only and multi-mapper stats modes with per-stage unique/multi/secondary breakdowns.
- Quality-filter settings exposed as YAML keys: `unique_only`, `samtools_filter_arguments`, `samtools_count_arguments`, `mapping_quality_cutoff` (per route).
- Profiles `local`, `conda`, `apptainer`, `docker`, `test`


## 1.0.0

- Added UMI support.

## 0.0.0

- Initial release.
