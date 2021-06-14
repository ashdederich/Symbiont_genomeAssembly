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
    * All the programs needed to complete genome assemblies in this nextflow file can be located in the container I created: https://cloud.sylabs.io/library/ashdederich/default/wholegenomeassembly
        - You can download this container by running this command: 
        ```
        singularity pull library://ashdederich/default/wholegenomeassembly
        ```
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

How does this work?
- Trimmomatic inputs your forward and reverse reads files and trims them based off of quality scores. My specifications that I have included have Trimmomatic trimming the 3 bases at the beginning and end of each read, as sequencers tend to have more sequencing errors at the ends of reads. Trimmomatic also creates a 'read window' of 4 bases and will trim that window if its quality score is at or below 15.
- bwa aligns your raw, trimmed reads to a reference genome, of which you should have downloaded and specified in the path. 
- samtools converts your .sam file output from bwa to a .bam, then to a .fastq for input into SPAdes.
- SPAdes takes the raw reads mapped to your reference genome and creates scaffolds from reads that are overlapping.
- GATB MindTheGap takes your scaffolds from SPAdes, in addition to your raw reads, and uses a DeBruijn graph to fill gaps in between contigs.

---*NOTE*--- For memory efficiency, I keep all files gzipped.

Most importantly, how should you execute all files and assemble a genome?
- Your files should be in some form of hierarchy, as follows:
    * The top directory, which contains ALL further directories and files. We could name this 'genome_assemblies'. This directory will contain:
        - Your reads directory, which contains all forward and reverse reads. This is also the directory that output from each program will be stored.
        - Your nextflow documents directory. This contains:
            * Your Singularity files
            * Your nextflow.nf files
            * Your nextflow.config file
            * Your temporary work directory (that nextflow automatically creates when a job is run)
            * Your bin folder with all extra scripts that you may use.
        - Your reference genome folder. This should have separate folders for each reference genome. Each reference genome directory contained in this directory will have:
            * Your .fasta file for the reference genome.
            * Your indexed files of this reference genome. This is to save time and memory during the assembly proceess. To index your reference genome, go to the correct reference genome directory and type: '''bwa index <reference_genome.fasta>'''
- Once your directories are created, go to your nextflow documents directory and run: '''nextflow run SymbiontGenomeAssembly.nf'''
    * If your nextflow job is ever interrupted, and you would like to pick up from where you left off, you can add the resume tag: '''nextflow run SymbiontGenomeAssembly.nf -resume'''

