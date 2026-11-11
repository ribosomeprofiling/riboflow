// STAR genome alignment for the RNA-seq path. SE/PE aware. Mirrors star_align.nf
// but uses params.star.rnaseq_arguments and names outputs with the
// rnaseq.genome_alignment.* prefix (storeDir separation). With
// star.output_transcriptome_bam it also emits the TranscriptomeSAM projection.
process STAR_ALIGN_RNASEQ {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(reads)
    path genome_dir

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.bam"),                 emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.transcriptome_alignment.bam"), optional: true, emit: transcriptome_bam
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.aligned.fastq.gz"),    emit: aligned
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.unaligned.fastq.gz"),  emit: unaligned
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.log"),                 emit: log
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.secondary.count"),     emit: secondary_count

    script:
    def prefix       = "${meta.id}.${meta.lane}"
    def is_pe        = meta.single_end == false
    def reads_in     = is_pe ? "${reads[0]} ${reads[1]}" : "${reads}"
    def star_args    = task.ext.star_args ?: (params.star?.rnaseq_arguments ?: '')
    def emit_tx_bam    = Utils.as_bool((params.star ?: [:]).output_transcriptome_bam, false)
    def quant_mode_arg = emit_tx_bam ? '--quantMode TranscriptomeSAM \\\n        ' : ''
    def tx_bam_cmd     = emit_tx_bam \
        ? "mv star_out/Aligned.toTranscriptome.out.bam ${prefix}.rnaseq.transcriptome_alignment.bam" \
        : ''
    def sort_threads = Math.min(task.cpus as int, 8)
    // See star_align.nf: 30 GiB index reserve, sort gets the rest (floored at 8 GiB).
    def sort_ram     = task.memory ? Math.max(task.memory.toBytes() - 32212254720L, 8589934592L) : 34359738368L
    """
    set -o pipefail
    mkdir -p star_out
    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --readFilesIn ${reads_in} \\
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

    mv star_out/Aligned.sortedByCoord.out.bam ${prefix}.rnaseq.genome_alignment.bam
    ${tx_bam_cmd}
    samtools index -@ ${sort_threads} ${prefix}.rnaseq.genome_alignment.bam

    # Primary records only (-F 2308) so a multimapper appears once.
    samtools fastq -@ ${sort_threads} -F 2308 ${prefix}.rnaseq.genome_alignment.bam \\
        | pigz -p ${task.cpus} > ${prefix}.rnaseq.genome_alignment.aligned.fastq.gz

    if [ -f star_out/Unmapped.out.mate1 ]; then
        pigz -p ${task.cpus} -c star_out/Unmapped.out.mate1 > ${prefix}.rnaseq.genome_alignment.unaligned.fastq.gz
    else
        echo -n | gzip > ${prefix}.rnaseq.genome_alignment.unaligned.fastq.gz
    fi

    cp star_out/Log.final.out ${prefix}.rnaseq.genome_alignment.log

    samtools view -@ ${task.cpus} -c -f 256 ${prefix}.rnaseq.genome_alignment.bam \\
        > ${prefix}.rnaseq.genome_alignment.secondary.count
    """

    stub:
    def prefix      = "${meta.id}.${meta.lane}"
    def emit_tx_bam = Utils.as_bool((params.star ?: [:]).output_transcriptome_bam, false)
    def tx_stub     = emit_tx_bam ? "touch ${prefix}.rnaseq.transcriptome_alignment.bam" : ''
    """
    touch ${prefix}.rnaseq.genome_alignment.bam
    ${tx_stub}
    echo | gzip -c > ${prefix}.rnaseq.genome_alignment.aligned.fastq.gz
    echo | gzip -c > ${prefix}.rnaseq.genome_alignment.unaligned.fastq.gz
    touch ${prefix}.rnaseq.genome_alignment.log
    echo 0 > ${prefix}.rnaseq.genome_alignment.secondary.count
    """
}
