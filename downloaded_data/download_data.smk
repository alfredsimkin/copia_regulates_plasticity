configfile: 'input_files/config.yaml'
ruleorder: rename_nanopore_reads>download_nanopore_reads
ruleorder: rename_illumina_reads>download_illumina_reads
working_dir=os.getcwd()

def get_nanopore_name(wildcards):
	return 'output_files/nanopore_reads/'+config['nanopore_names'][wildcards.accession]+'.fastq'

def get_illumina_R1(wildcards):
	return 'output_files/illumina_reads/'+config['illumina_names'][wildcards.accession]+'_1.fastq'

def get_illumina_R2(wildcards):
	return 'output_files/illumina_reads/'+config['illumina_names'][wildcards.accession]+'_2.fastq'

def get_old_name(wildcards):
	print('new is', wildcards.newname)
	return config['links'][wildcards.newname]

rule all:
	input:
		done_linking='output_files/done_linking.txt'

rule download_flair:
	output:
		'flair_installed.txt'
	shell:
		'''
		conda create -n flair -c conda-forge -c bioconda flair snakemake
		touch flair_installed.txt
		'''
		

rule download_dm6:
	output:
		dm6='output_files/dm6.fa'
	shell:
		'''
		wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
		gunzip dm6.fa.gz
		mv dm6.fa output_files
		'''

rule download_sra_toolkit:
	output:
		toolkit_file='output_files/sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump'
	shell:
		'''
		wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz
		tar -xzf sratoolkit.3.0.2-ubuntu64.tar.gz
		cp -R sratoolkit.3.0.2-ubuntu64 output_files
		rm -R sratoolkit.3.0.2-ubuntu64
		rm sratoolkit.3.0.2-ubuntu64.tar.gz
		'''
rule download_nanopore_reads:
	input:
		toolkit='output_files/sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump'
	output:
		file='output_files/nanopore_reads/{accession}.fastq'

	shell:
		'''
		output_files/sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump {wildcards.accession}
		mv {wildcards.accession}.fastq output_files/nanopore_reads
		'''

rule download_illumina_reads:
	input:
		toolkit='output_files/sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump'
	output:
		file1='output_files/illumina_reads/{accession}_1.fastq',
		file2='output_files/illumina_reads/{accession}_2.fastq'
	shell:
		'''
		output_files/sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump {wildcards.accession}
		mv {wildcards.accession}_1.fastq output_files/illumina_reads
		mv {wildcards.accession}_2.fastq output_files/illumina_reads
		'''

rule rename_nanopore_reads:
	input:
		input_name=get_nanopore_name
	output:
		output_name='output_files/nanopore_reads/{accession}.fastq'
	shell:
		'mv {input.input_name} {output.output_name}'

rule rename_illumina_reads:
	input:
		R1_input=get_illumina_R1,
		R2_input=get_illumina_R2
	output:
		R1_output='output_files/illumina_reads/{accession}_1.fastq',
		R2_output='output_files/illumina_reads/{accession}_2.fastq'
	shell:
		'''
		mv {input.R1_input} {output.R1_output}
		mv {input.R2_input} {output.R2_output}		
		'''

rule link_files:
	input:
		yaml_file='input_files/linking.yaml',
		nanopore=expand('output_files/nanopore_reads/{accession}.fastq', accession=config['nanopore_names']),
		illumina=expand('output_files/illumina_reads/{accession}_{number}.fastq', accession=config['illumina_names'], number=[1,2]),
		dm6='output_files/dm6.fa'
	output:
		done_linking='output_files/done_linking.txt'
	script:
		'input_files/scripts/link_files.py'
