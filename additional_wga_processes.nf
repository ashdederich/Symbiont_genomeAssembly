process blastx {

    publishDir 'genomedata/${sample_id}/SymbiontAssembly/${params.reference}/', mode: 'copy'

    input:
    //change this to the correct name of the ref genome, also create a parameter for it. Do it the same way you did for specifying the referencee geneome in the parameter section up top
    tuple val(sample_id), file "${sample_id}_scaffolds.fasta.gz" from spades_ch
    path db from params.whatevernameis

    output:
    tuple val(sample_id), file("${sample_id}_blastxout.xml") into blastx_ch

    script:
    """
    blastx \
    -query ${sample_id}_scaffolds.fasta.gz \
    -db ${db} \
    -out ${sample_id}_blastxout.xml \
    -outfmt 5
    """
}

process rps {

    publishDir 'genomedata/${sample_id}/SymbiontAssembly/${params.reference}/', mode: 'copy'

    input:
    tuple val(sample_id), file("${sample_id}_blastxout.xml") from blastx_ch

    output:
    tuple val(sample_id), file("${sample_id}_rps_parsed_blast.txt") into parsed_blast_ch

    script:
    """
    rps.py \
    ${sample_id}_blastxout.xml > \
    ${sample_id}_rps_parsed_initblast.txt
    """
}

process pullScaffolds {

    publishDir 'genomedata/${sample_id}/SymbiontAssembly/${params.reference}/', mode: 'copy'
   
    input:
    tuple val(sample_id), file("${sample_id}_scaffolds.fasta.gz") from spades_ch
    tuple val(sample_id), file("${sample_id}_rps_parsed_blast.txt") from parsed_blast_ch

    output:
    tuple val(sample_id), file("${sample_id}_scaffold_candidates.fasta.gz") into scaffold_candidates_ch

    script:
    """
    gunzip ${sample_id}_scaffolds.fasta.gz
    grab_seq.pl ${sample_id}_scaffolds.fasta ${sample_id}_rps_parsed_blast.txt > ${sample_id}_scaffold_candidates.fasta
    bgzip -@ ${task.cpus} ${sample_id}_scaffold_candidates.fasta
    bgzip -@ ${task.cpus} ${sample_id}_scaffolds.fasta
    """
}