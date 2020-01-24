# FAQ

1. **Does RiboFlow work with genomic reference?**  
No, RiboFlow is designed to work with transcriptomic references only.
If a gene has alternative isoforms, typically exactly one isoform is included in the reference transcriptome.
Optionally, RiboFlow can map the reads, which didn't map to the transcriptome, using HISAT2 to the genome. Yet these mapping results are only used for diagnostic purposes and they have no effect on the resulting ribo file. Hence thoese reads don't change the downstream analysis.

2. **Where can I find references for organisms other than human?**  
Here is a [repository](https://github.com/ribosomeprofiling/references_for_riboflow)
of officially supported RiboFlow references for various organisms.
We plan to add more organisms in the future.

3. **How can I prepare an alternative transcriptomic reference for RiboFlow?**  
You can take a look at our [notebook](https://github.com/ribosomeprofiling/references_for_riboflow/blob/master/transcriptome/human/v1/scripts/appris_explore.ipynb)
where we explain how we built the human refernece transcriptome for RiboFlow.
[That repository](https://github.com/ribosomeprofiling/references_for_riboflow) also contains the scripts used to generate the reference files. 

4. **I have Unique Molecular Identifiers (UMIs) in my sequencing data. Can RiboFlow use UMIs to deduplicate the reads?**  
While UMIs are NOT supported at the moment, we are working on adding UMI support to RiboFlow in the coming months.

5. **Does RiboFlow provide PCR deduplication?**  
Yes, RiboFlow supports deduplication based on read alignment position. More explicitly, if the two (or more) reads have the same length and they map to the exact same position, those reads will be considered as duplicates and collapsed into one read. Deduplication can be turned on or off by setting the parameter `deduplicate` in the [parameters file](https://github.com/ribosomeprofiling/riboflow/blob/master/project.yaml)

6. **Can I use RNA-Seq data with RiboFlow?**  
Yes, you can. In the  [parameters file](https://github.com/ribosomeprofiling/riboflow/blob/master/project.yaml) file, set `do_rnaseq: true` and under `rnaseq` node (see the example in the parameters file) you can pair ribosome profiling data with RNA-Seq data. Please note that RiboFlow processes RNA-Seq data and ribosome profiling data in a parallel fashion. Hence, RNA-Seq reads are mapped in single end mode. If your data is coming from paired-end sequencing, you can simply provide the reads coming from the first round of sequencing (first file from each pair).

7. **Does RiboFlow work with references where genes don't have 3'UTRs?**  
No, all reference genes must have their 3' UTRs annotated.

8. **In my particular annotation, thre is a set of genes missing 3' or 5' UTRs. How can I incorporate these genes in the RiboFlow reference build?**  
We suggest the following simple workaround for this problem. You can take a fixed number of (say, for example, 50) nucleotides preceeding / proceeding the start / stop site in the genomic sequence. Attach those fragments to your transcript sequence and annotate them as 3' / 5' UTRs. Then you can complile this transcriptome sequence + annotation as a valid RiboFlow reference.
