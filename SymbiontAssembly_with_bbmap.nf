#!/usr/bin/env nextflow

params.reads="/path/to/forward-and-reverse/reads/*_{R1,R2}_001.fastq.gz"
reads_ch=Channel.fromFilePairs(params.reads)
reads_ch.into {reads_trimmomatic; reads_mindTheGap }
params.reference='buchnera'
params.ref_genome="/path/to/reference-genomes/in/own/directory/RefGenomes/${params.reference}/*.fna.gz"
params.workdir="~"
workdir_ch=Channel.fromPath(params.workdir)

process trimmomatic {
    publishDir workdir_ch/"Trimmomatic", mode: 'copy'
    tag "${sample_id}"
    label params.label
    label 'wga' 

    input:
    tuple val(sample_id), path(sample_files) from reads_trimmomatic
    
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

process pear {
    tag "${sample_id}"
    label 'alignment'
    label 'pear'

    input:
    tuple val(sample_id), file("${sample_id}_R1_paired.fastq.gz"), file("${sample_id}_R2_paired.fastq.gz") from trimmed_paired_reads_ch

    output:
    tuple val(sample_id), file("${sample_id}_R1R2_PEARoutput.unassembled.forward.fastq"), file("${sample_id}_R1R2_PEARoutput.unassembled.reverse.fastq") into pear_unassembled
    tuple val(sample_id), file("${sample_id}_R1R2_PEARoutput.assembled.fastq") into pear_ch_assembled
    tuple val(sample_id), file ("${sample_id}_R1R2_PEARoutput.discarded.fastq") into pear_ch_discarded

    script:
    """
    pear -j ${task.cpus} \
    -f ${sample_id}_R1_paired.fastq.gz \
    -r ${sample_id}_R2_paired.fastq.gz \
    -o ${sample_id}_R1R2_PEARoutput
    """
}

process compressPear{
    publishDir workdir_ch/"Trimmomatic/PEAR", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label params.label
    label 'tabix'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mapped.bam") from samToBam_ch
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mapped.fastq") from samToBamfastq_ch

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mapped.fastq.gz") into samToFastq_ch_gz

    script:
    """
    bgzip -@ ${task.cpus} ${sample_id}_${params.reference}_bbmap_mapped.fastq
    """
}

process scaffoldAlignment_bbmap {
    tag "${sample_id}"
    tag "${params.reference}"
    label params.label
    label 'bbmap_container'

    input:
    tuple val(sample_id), file("${sample_id}_R1R2_PEARoutput.assembled.fastq") from pear_ch_assembled
    file ref_genome from ref_ch
    
    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mappedreads.sam") into bbmapcov_ch

    script:
    """
    bbmap.sh in=${sample_id}_R1R2_PEARoutput.assembled.fastq ref=${ref_genome} t=1 out=${sample_id}_${params.reference}_bbmap_mappedreads.sam
    """
}

process samToBam {

    publishDir "/uufs/chpc.utah.edu/common/home/vondolen-group1/genomedata/${sample_id}/SymbiontAssembly/${params.reference}/Alignment", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label params.label
    label 'wga'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mappedreads.sam") from bbmapcov_ch

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mapped.bam") into samToBam_ch
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mapped.fastq") into samToBamfastq_ch

    script:
    """
    samtools view -@ ${task.cpus} -Sbu -o ${sample_id}_${params.reference}_bbmap_mapped.bam ${sample_id}_${params.reference}_bbmap_mappedreads.sam
    samtools fastq -@ ${task.cpus} ${sample_id}_${params.reference}_bbmap_mapped.bam > ${sample_id}_${params.reference}_bbmap_mapped.fastq
    """
}

process compressSamToFastq {
    publishDir "/path/to/forward-and-reverse/reads/${sample_id}/SymbiontAssembly/${params.reference}/Alignment", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label params.label
    label 'tabix'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mapped.bam") from samToBam_ch
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mapped.fastq") from samToBamfastq_ch

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mapped.fastq.gz") into samToFastq_ch_gz

    script:
    """
    bgzip -@ ${task.cpus} ${sample_id}_${params.reference}_bbmap_mapped.fastq
    """
}

process spades {
    publishDir "/mnt/genomedata/Nextflow_assemblies/sequences/${sample_id}/SymbiontAssembly/${params.reference}/Spades", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    cpus = 28
    label 'wga'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_bbmap_mapped.fastq.gz") from samToFastq_ch_gz
    tuple val(sample_id), file("${sample_id}_R1R2_PEARoutput.unassembled.forward.fastq"), file("${sample_id}_R1R2_PEARoutput.unassembled.reverse.fastq") into pear_unassembled
    tuple val(sample_id), file ("${sample_id}_R1R2_PEARoutput.discarded.fastq") into pear_ch_discarded


    output:
    tuple val(sample_id), path("${sample_id}_SPADESoutput/${sample_id}_${params.reference}_scaffolds.fasta") into spades_ch

    script:
    """
    spades -m 120 -t ${task.cpus} -k 21,33,55,77,99,127 \
    -1 ${sample_id}_R1R2_PEARoutput.unassembled.forward.fastq.gz \
    -2 ${sample_id}_R1R2_PEARoutput.unassembled.reverse.fastq.gz \
    --merged ${sample_id}_${params.reference}_bbmap_mapped.fastq.gz \
    -s ${sample_id}_R1R2_PEARoutput.discarded.fastq.gz \
    -o ${sample_id}_SPADESoutput
    mv ${sample_id}_SPADESoutput/scaffolds.fasta ${sample_id}_SPADESoutput/${sample_id}_${params.reference}_scaffolds.fasta
    """
}

process compressSPADES {
    publishDir "/mnt/genomedata/Nextflow_assemblies/sequences/${sample_id}/SymbiontAssembly/${params.reference}/Spades", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label params.label
    label 'wga'

    input:
    tuple val(sample_id), path("${sample_id}_SPADESoutput/${sample_id}_${params.reference}_scaffolds.fasta") from spades_ch

    output:
    tuple val(sample_id), path("${sample_id}_SPADESoutput/${sample_id}_${params.reference}_scaffolds.fasta.gz") into spades_ch_gz

    script:
    """
    bgzip -@ ${task.cpus} ${sample_id}_SPADESoutput/${sample_id}_${params.reference}_scaffolds.fasta
    """
}

process GapFiller {
    tag "${sample_id}"
    tag "${params.reference}"
    label params.label
    label 'mtg'

    input:
    tuple val(sample_id), path("${sample_id}_SPADESoutput/${sample_id}_${params.reference}_scaffolds.fasta.gz") from spades_ch_gz
    tuple val(sample_id), path(sample_files) from reads_mindTheGap

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_gapFilled.insertions.fasta"),file("${sample_id}_${params.reference}_gapFilled.gfa"), file("${sample_id}_${params.reference}_gapFilled.info.txt"), file("${sample_id}_${params.reference}_gapFilled.h5"), file("${sample_id}_${params.reference}_gapFilled_seed_dictionary.fasta") into gap_filler_ch

    script:
    """
    cat $sample_files > ${sample_id}_reads.fastq.gz

    MindTheGap fill \
    -nb-cores ${task.cpus} \
    -in ${sample_id}_reads.fastq.gz \
    -contig ${sample_id}_${params.reference}_assembly.fasta.gz \
    -kmer-size 51 \
    -abundance-min 5 \
    -max-nodes 300 \
    -max-length 50000 \
    -out ${sample_id}_${params.reference}_gapFilled
    """
}

process compressGapFiller {
    publishDir "/path/to/forward-and-reverse/reads/${sample_id}/SymbiontAssembly/${params.reference}/GapFilled", mode: 'copy'
    tag "${sample_id}"
    tag "${params.reference}"
    label params.label
    label 'tabix'

    input:
    tuple val(sample_id), file("${sample_id}_${params.reference}_gapFilled.insertions.fasta"),file("${sample_id}_${params.reference}_gapFilled.gfa"), file("${sample_id}_${params.reference}_gapFilled_seed_dictionary.fasta") from gap_filler_ch

    output:
    tuple val(sample_id), file("${sample_id}_${params.reference}_gapFilled.insertions.fasta.gz"),file("${sample_id}_${params.reference}_gapFilled.gfa.gz"), file("${sample_id}_${params.reference}_gapFilled_seed_dictionary.fasta.gz") into gap_filler_ch_gz

    script:
    """
    bgzip -@ ${task.cpus} ${sample_id}_${params.reference}_gapFilled.insertions.fasta
    bgzip -@ ${task.cpus} ${sample_id}_${params.reference}_gapFilled.gfa
    bgzip -@ ${task.cpus} ${sample_id}_${params.reference}_gapFilled_seed_dictionary.fasta
    """
}