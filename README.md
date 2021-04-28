What is this?

This is my Nextflow documentation for assembling symbiont genomes. Information about Nextflow can be found at: https://www.nextflow.io/. All documents needed to run Nextflow will automate and parallelize genome assemblies.

What programs are needed to run this?
- Nextflow
- Singularity OR Docker
- Preferably a job scheduler, such as Slurm, but the local computer can be used as well.

What files are needed for each program?
- Nextflow:
    * A Nextflow .nf file, which is uploaded on my GitHub at https://github.com/ashdederich/Symbiont_genomeAssembly/SymbiontGenomeAssembly.nf
        - The SymbiontGenomeAssembly.nf file paths to forward/reverse reads, reference genome(s), and directory to copy the output files in also needs to be updated to where ~your~ files are.
    * A nextflow.config file, which is also uploaded on my GitHub at: https://github.com/ashdederich/Symbiont_genomeAssembly/nextflow.config
        - The nextflow.config file paths to the Singulariy/Docker containers need to be updated to where ~your~ Singularity/Docker containers are located.
- Singularity OR Docker:
    * Most of the programs needed, such as Trimmomatic, bwa, and samtools can be located in the container I created: https://cloud.sylabs.io/library/ashdederich/default/wholegenomeassembly
    * The tabix container, used to gzip files, can be downloaded from: https://quay.io/repository/biocontainers/tabix?tab=tags
    * The Minia pipeline container for de-novo assembly of contigs at multiple k-mer lengths can be downloaded from: https://hub.docker.com/r/cimendes/gatb-minia-pipeline
    * The MindTheGap container for gap-filling between contigs can be downloaded from: https://quay.io/repository/biocontainers/mindthegap?tab=tags
- Slurm:
    * Slurm, or the job scheduler of choice, only needs to be specified in the nextflow.config file


What should you keep in mind when running this?
- All nextflow files, sequencing files, reference genome(s), and Singularity/Docker containers need to be within the same parent directory.
- You need to edit the Nextflow files (both nextflow.config and the .nf file) to contain the correct paths for your own system.
- I have it set for Nextflow to create an output directory and copy output files to specific folders within that directory.
    * Each directory is named by the specific sample ID of each sequence pair.
    * There are nested folders that I have Nextflow create to organize the file outputs. These can be changed as desired as well.
    * Nextflow names files/folders according to the sample ID and type of symbiont, as specified by me.
- Reference files need to be in their own directory, and each symbiont reference genome needs to be in its own file. Additionally, all reference genomes need to be first indexed by bwa, and these files should remain in the same directory as their corresponding reference genome.
    * You can index these reference genomes by running: '''bwa index [reference_genome.fasta.gz]'''
- You can change the symbiont genome as needed through the command prompt, when calling nextflow.
    * Currently, it is set to 'buchnera' within the .nf file. 
    * However, if you wanted to sequence another symbiont genome, such as Wolbachia, you would change your command prompt to: '''nextflow run SymbiontGenomeAssembly.nf --reference 'wolbachia' '''
    * This will updated a few things:
        - The file and folder names created will be updated to reflect that you are assembling a different symbiont.
        - The location of the reference folder and all necessary reference files will be updated.
    * This is useful if there are several symbiont genomes within the same organism that you wish to assemble and want all files to be located under the same sequence ID.