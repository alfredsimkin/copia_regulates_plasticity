About:  
This repository is for people who would like to reproduce the computational results of our manuscript "The Retrotransposon Copia regulates structural synaptic plasticity at the Drosophila Neuromuscular Junction".  
  
This repository assumes users who have sudo access to an ubuntu linux distribution (either physical or in a virtual machine) and assumes some knowledge of Unix (e.g. modifying .bashrc files to modify your $PATH variable to point to binary files and/or copying binary files into an existing $PATH folder).  
  
This repository also makes the implicit assumption that users who want to know how any given result was gathered can follow the flow of a snakefile, and have some knowledge of Python scripting. Please write to me if you have questions!  
  
Setup:  
  
0. Download or clone this repository and unzip it on an Ubuntu machine, and cd to it.  
  
1. Obtain snakemake:  
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html  
I followed the instructions, and obtained mambaforge by downloading a shell script for my linux distro from here:  
https://github.com/conda-forge/miniforge#mambaforge  
Then I obtained snakemake with this command:  
mamba create -c conda-forge -c bioconda -n snakemake snakemake  
I tested this on snakemake version 7.19.1  
  
2. Download long sequencing reads and illumina reads, and setup renamed shortcuts to files:  
cd to the downloaded_data folder  
conda activate snakemake  
snakemake -s download_data.smk --cores {your_desired_core_count}  
  
3. Obtain FLAIR from here:  
https://flair.readthedocs.io/en/latest/  
I used FLAIR version 1.7.0 for testing and installed snakemake alongside FLAIR in a single conda environment with this command:  
conda create -n flair -c conda-forge -c bioconda flair snakemake  
  
4. Obtain blat from here:  
http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/  
copy the executable file to your $PATH  
  
5. Obtain the dm6 genome from here (I use wget for this):  
https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz  
Unzip it with gunzip dm6.fa.gz  
  
6. Obtain minimap2:  
install curl (if needed) with (e.g.) sudo apt install curl  
Obtain minimap2 binaries from here:  
curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -  
cd to minimap2-2.17_x64-linux  
copy the minimap2 executable to your $PATH  
  
7. Obtain salmon:  
sudo apt install salmon  
I used salmon 1.4.0  
  
8. Obtain scipy and plotly:  
conda activate snakemake  
conda install scipy  
conda install plotly  
  
9. Obtain kaleido:  
sudo apt install python-is-python3  
sudo apt install pip  
pip install -U kaleido  
  
Reproducing CNS vs. BWM isoform abundances results:  
1. cd to the FLAIR_isoforms folder  
2. conda activate flair  
3. snakemake -s run_flair.smk --cores {your_desired_core_count}  
4. examine FLAIR_isoforms/output_files/mapped_reads.counts.tsv  
  
Reconstructing TE seeds:  
1. cd to the constuct_TE_seeds folder  
2. conda activate snakemake  
3. snakemake -s construct_seeds.smk --cores {your_desired_core_count}  
4. examine construct_TE_seeds/output_files/dm6_TE_seeds.fa  
  
measuring copia chimeric long reads (requires "reconstructing TE seeds" as a precursor):  
1. cd to the chimeric_copia_reads folder  
2. conda activate snakemake  
3. snakemake -s chimeric_copia.smk --cores {your_desired_core_count)  
4. examine these output files:  
chimeric_copia_reads/output_files/CNS_final_summary.tsv - see especially column B of the 7 copia reads that map to refseq genes and RNAi regions, around row 1275  
chimeric_copia_reads/output_files/seed_mapping/CNS_seed_results.tsv (there are 1225 reads that have a >100 nucleotide match to consensus copia, but only 1224 of them have a >100 nucleotide match to the genome)  
  
measuring skews in gene expression following RNAi:  
1. cd to the RNAi_skewed_gene_expression folder  
2. conda activate snakemake  
3. snakemake -s modular_snakefile.smk  
4. examine RNAi_skewed_gene_expression/output_files/summarized_final_stats folder  
  
searching for RNAi regions in mature and immature Refseq sequences (requires "reconstructing TE seeds" and "measuring copia chimeric long reads" as precursors):  
1. cd to the RNAi_refseq_genomic folder  
2. conda activate snakemake  
3. snakemake -s RNAi_refseq_genomic.smk --cores {your_desired_core_count}  
4. examine these output files:  
RNAi_refseq_genomic/output_files/summaries/dm6_genes_containment_stats.tsv  
RNAi_refseq_genomic/output_files/summaries/dm6_genes_plus_introns_containment_stats.tsv  
  
searching for copia-mapping regions of the genome:  
1. cd to the copia_consensus_genome_mappers folder  
2. conda activate snakemake  
3. snakemake -s copia_genomic.smk --cores {your_desired_core_count}  
4. examine these output files:  
copia_consensus_genome_mappers/output_files/region_stats.tsv  
copia_consensus_genome_mappers/output_files/summarized_genomic_copia.tsv  
  

