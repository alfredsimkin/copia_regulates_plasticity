Setup:

1. Obtain snakemake:
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
I followed the instructions, and obtained mambaforge by downloading a shell
script for my linux distro from here:
https://github.com/conda-forge/miniforge#mambaforge
Then I obtained snakemake with this command:
mamba create -c conda-forge -c bioconda -n snakemake snakemake
I tested this on snakemake version 7.19.1

2. Obtain Nanopore long sequencing reads from here:


3. Obtain FLAIR from here:
https://flair.readthedocs.io/en/latest/
I used FLAIR version 1.7.0 for testing and installed snakemake alongside FLAIR
in a single conda environment with this command:
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
curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf -
cd to minimap2-2.24_x64-linux
copy the minimap2 executable to your $PATH

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

measuring copia chimeric long reads:
1. cd to the chimeric_copia_reads folder
2. conda activate snakemake
3. snakemake -s chimeric_copia.smk --cores {your_desired_core_count)
4. examine copia_chimeric_TE_seeds/output_files/CNS_final_summary.tsv

measuring skews in gene expression following RNAi:


searching for RNAi regions in mature and immature Refseq sequences:
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


