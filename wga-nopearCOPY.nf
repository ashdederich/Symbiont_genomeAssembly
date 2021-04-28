#!/usr/bin/env nextflow

params.reads="/uufs/chpc.utah.edu/common/home/vondolen-group1/genomedata/*_{R1,R2}_001.fastq.gz"
reads_ch=Channel.fromFilePairs(params.reads)
params.unzippedReads="/uufs/chpc.utah.edu/common/home/vondolen-group1/genomedata/*_{R1,R2}_001.fq"
unzippedReads_ch=Channel.fromFilePairs(params.unzippedReads)
params.reference='buchnera'
params.ref_genome="/uufs/chpc.utah.edu/common/home/vondolen-group1/RefGenomes/${params.reference}/*.fna.gz"
params.bwa_amb="/uufs/chpc.utah.edu/common/home/vondolen-group1/RefGenomes/${params.reference}/*.fna.gz.amb"
params.bwa_ann="/uufs/chpc.utah.edu/common/home/vondolen-group1/RefGenomes/${params.reference}/*.fna.gz.ann"
params.bwa_bwt="/uufs/chpc.utah.edu/common/home/vondolen-group1/RefGenomes/${params.reference}/*.fna.gz.bwt"
params.bwa_pac="/uufs/chpc.utah.edu/common/home/vondolen-group1/RefGenomes/${params.reference}/*.fna.gz.pac"
params.bwa_sa="/uufs/chpc.utah.edu/common/home/vondolen-group1/RefGenomes/${params.reference}/*.fna.gz.sa"

process trimmomatic {
    publishDir "/uufs/chpc.utah.edu/common/home/vondolen-group1/genomedata/${sample_id}/InitialAssembly/Trimmed/", mode: 'copy'
    tag "${sample_id}"
    label 'alignment'
    label 'wga' 

    input:
    tuple val(sample_id), path(sample_files) from reads_ch
    
    output:
    tuple val(sample_id), file("${sample_id}_R1_paired.fastq.gz"), file("${sample_id}_R2_paired.fastq.gz") into trimmed_paired_reads_ch
    tuple val(sample_id), file("${sample_id}_R1_unpaired.fastq.gz"), file("${sample_id}_R2_unpaired.fastq.gz") into trimmed_unpaired_reads_ch
    tuple val(sample_id), file("${sample_id}_Trimlog.txt")

    script:
    """
    TrimmomaticPE -threads ${task.cpus} \
    -trimlog ${sample_id}_Trimlog.txt \
    $sample_files \
    ${sample_id}_R1_paired.fastq.gz \
    ${sample_id}_R1_unpaired.fastq.gz \
    ${sample_id}_R2_paired.fastq.gz \
    ${sample_id}_R2_unpaired.fastq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process scaffoldAlignment {
    publishDir "/uufs/chpc.utah.edu/common/home/vondolen-group1/genomedata/${sample_id}/SymbiontAssembly/${params.reference}/Alignment", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label 'alignment'
    label 'wga'

    input:
    tuple val(sample_id), file("${sample_id}_R1_paired.fastq.gz"), file("${sample_id}_R2_paired.fastq.gz") from trimmed_paired_reads_ch
    tuple val(sample_id), file("${sample_id}_R1_unpaired.fastq.gz"), file("${sample_id}_R2_unpaired.fastq.gz") from trimmed_unpaired_reads_ch
    path ref_genome from params.ref_genome
    path amb from params.bwa_amb
    path ann from params.bwa_ann
    path bwt from params.bwa_bwt
    path pac from params.bwa_pac
    path sa from params.bwa_sa
    
    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_aligned.sam") into bwacov_ch

    script:
    """
    bwa mem -p -t ${task.cpus} \
    ${ref_genome} \
    <(cat ${sample_id}_R1_paired.fastq.gz ${sample_id}_R2_paired.fastq.gz) <(cat ${sample_id}_R1_unpaired.fastq.gz ${sample_id}_R2_unpaired.fastq.gz) \
    > ${sample_id}_${params.reference}_aligned.sam 
    """
}

process samToBam {

    publishDir "/uufs/chpc.utah.edu/common/home/vondolen-group1/genomedata/${sample_id}/SymbiontAssembly/${params.reference}/Alignment", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label 'alignment'
    label 'wga'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_aligned.sam") from bwacov_ch

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_aligned_sort.bam"), file("${sample_id}_${params.reference}_aligned_sort.bam.bai"),file("${sample_id}_${params.reference}_aligned.fastq") into samToBam_ch

    script:
    """
    samtools view -@ ${task.cpus} -Sbu -o ${sample_id}_${params.reference}_aligned.bam ${sample_id}_${params.reference}_aligned.sam
    samtools sort -@ ${task.cpus} -o ${sample_id}_${params.reference}_aligned_sort.bam ${sample_id}_${params.reference}_aligned.bam
    samtools index -@ ${task.cpus} ${sample_id}_${params.reference}_aligned_sort.bam 
    samtools fastq -@ ${task.cpus} ${sample_id}_${params.reference}_aligned.bam > ${sample_id}_${params.reference}_aligned.fastq
    """
}

process compressSamToFastq {
    publishDir "/uufs/chpc.utah.edu/common/home/vondolen-group1/genomedata/${sample_id}/SymbiontAssembly/${params.reference}/Alignment", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label 'alignment'
    label 'tabix'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_aligned_sort.bam"), file("${sample_id}_${params.reference}_aligned_sort.bam.bai"),file("${sample_id}_${params.reference}_aligned.fastq") from samToBam_ch

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_aligned.fastq.gz") into samToFastq_ch_gz

    script:
    """
    bgzip -@ ${task.cpus} ${sample_id}_${params.reference}_aligned.fastq
    """
}


process minia{
    tag "${sample_id}"
    label 'alignment'
    tag "${params.reference}"
    label 'GATB'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_aligned.fastq.gz") from samToFastq_ch_gz

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_assembly.fa") into minia_deNovo_ch

    script:
    """
    gatb --kmer-sizes 21,33,55,77,99,127 --nb-cores ${task.cpus} -s ${sample_id}_${params.reference}_aligned.fastq.gz

    mv assembly.fasta ${sample_id}_${params.reference}_assembly.fa
    """
}

process GapFiller {
    publishDir "/uufs/chpc.utah.edu/common/home/vondolen-group1/genomedata/${sample_id}/SymbiontAssembly//${params.reference}/GapFilled", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label 'alignment'
    label 'mtg'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_assembly.fa") from minia_deNovo_ch
    tuple val(sample_id), path(sample_files) from unzippedReads_ch

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_gapFilledout.insertions.fa"),file("${sample_id}_${params.reference}_gapFilledout.gfa"), file("${sample_id}_${params.reference}_gapFilledout.info.txt") into gap_filler_ch

    script:
    """
    cat $sample_files > ${sample_id}_reads.fq

    MindTheGap fill \
    -nb-cores ${task.cpus} \
    -in ${sample_id}_reads.fq \
    -contig ${sample_id}_${params.reference}_assembly.fa \
    -kmer-size 51 \
    -abundance-min 5 \
    -max-nodes 300 \
    -max-length 50000 \
    -out ${sample_id}_${params.reference}_gapFilled
    """
}

process compressGapFiller {
    publishDir "/uufs/chpc.utah.edu/common/home/vondolen-group1/genomedata/${sample_id}/SymbiontAssembly//${params.reference}/GapFilled", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label 'alignment'
    label 'tabix'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_gapFilledout.insertions.fa"),file("${sample_id}_${params.reference}_gapFilledout.gfa") from gap_filler_ch

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_gapFilledout.insertions.fa.gz"),file("${sample_id}_${params.reference}_gapFilledout.gfa.gz") into gap_filler_ch_gz

    script:
    """
    bgzip -@ ${task.cpus} ${sample_id}_${params.reference}_gapFilledout.insertions.fa
    bgzip -@ ${task.cpus} ${sample_id}_${params.reference}_gapFilledout.gfa
    """
}