// STAR genome alignment. Ports `genome_alignment` (RiboFlow.groovy:401-483).
// Optionally emits a transcriptome-coordinate BAM (star.output_transcriptome_bam).
process STAR_ALIGN {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(fastq)
    path genome_dir

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.bam"),                   emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.transcriptome_alignment.bam"), optional: true, emit: transcriptome_bam
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.aligned.fastq.gz"),      emit: aligned
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.unaligned.fastq.gz"),    emit: unaligned
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.log"),                   emit: log
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.secondary.count"),       emit: secondary_count

    script:
    def prefix         = "${meta.id}.${meta.lane}"
    def emit_tx_bam    = Utils.as_bool((params.star ?: [:]).output_transcriptome_bam, false)
    def quant_mode_arg = emit_tx_bam ? '--quantMode TranscriptomeSAM \\\n        ' : ''
    def tx_bam_cmd     = emit_tx_bam \
        ? "mv star_out/Aligned.toTranscriptome.out.bam ${prefix}.transcriptome_alignment.bam" \
        : ''
    def sort_threads   = Math.min(task.cpus as int, 8)
    def star_args      = task.ext.star_args ?: (params.star?.ribo_arguments ?: '')
    // STAR holds the genome index in RAM while it sorts the BAM. Reserve 30 GiB for
    // a human-sized index (measured in benchmarks/fixed_resources.config) and give
    // the sort the rest, floored at 8 GiB. Falls back to 32 GiB when task.memory is unset.
    def sort_ram       = task.memory ? Math.max(task.memory.toBytes() - 32212254720L, 8589934592L) : 34359738368L
    """
    set -o pipefail
    mkdir -p star_out
    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --readFilesIn ${fastq} \\
        --readFilesCommand pigz -dc \\
        --readFilesType Fastx \\
        --genomeDir ${genome_dir} \\
        --genomeLoad NoSharedMemory \\
        ${star_args} \\
        --outSAMtype BAM SortedByCoordinate \\
        --limitBAMsortRAM ${sort_ram} \\
        ${quant_mode_arg}--outSAMattributes All \\
        --outSAMstrandField intronMotif \\
        --outSAMattrRGline ID:${prefix} SM:${meta.id} PL:ILLUMINA \\
        --outReadsUnmapped Fastx \\
        --outFileNamePrefix star_out/

    # STAR emitted this coordinate-sorted (--outSAMtype BAM SortedByCoordinate).
    mv star_out/Aligned.sortedByCoord.out.bam ${prefix}.genome_alignment.bam
    samtools index -@ ${sort_threads} ${prefix}.genome_alignment.bam

    ${tx_bam_cmd}

    # Aligned reads → FASTQ, primary records only (-F 2308) so a multimapper
    # appears once. Kept on disk for inspection / optional FastQC.
    samtools fastq -@ ${sort_threads} -F 2308 ${prefix}.genome_alignment.bam \\
        | pigz -p ${task.cpus} > ${prefix}.genome_alignment.aligned.fastq.gz

    if [ -f star_out/Unmapped.out.mate1 ]; then
        pigz -p ${task.cpus} -c star_out/Unmapped.out.mate1 > ${prefix}.genome_alignment.unaligned.fastq.gz
    else
        echo -n | gzip > ${prefix}.genome_alignment.unaligned.fastq.gz
    fi

    cp star_out/Log.final.out ${prefix}.genome_alignment.log

    samtools view -@ ${task.cpus} -c -f 256 ${prefix}.genome_alignment.bam \\
        > ${prefix}.genome_alignment.secondary.count
    """

    stub:
    def prefix      = "${meta.id}.${meta.lane}"
    def emit_tx_bam = Utils.as_bool((params.star ?: [:]).output_transcriptome_bam, false)
    def tx_stub     = emit_tx_bam ? "touch ${prefix}.transcriptome_alignment.bam" : ''
    """
    touch ${prefix}.genome_alignment.bam
    ${tx_stub}
    echo | gzip -c > ${prefix}.genome_alignment.aligned.fastq.gz
    echo | gzip -c > ${prefix}.genome_alignment.unaligned.fastq.gz
    touch ${prefix}.genome_alignment.log
    echo 0 > ${prefix}.genome_alignment.secondary.count
    """
}
