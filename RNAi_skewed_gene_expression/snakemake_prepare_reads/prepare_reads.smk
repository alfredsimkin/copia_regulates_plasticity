configfile: 'linking.yaml'

def get_R1_files(wildcards):
	return expand('input_files/fastq_files/'+wildcards.type+'_'+wildcards.subtype+'_{rep}_1.fastq', rep=config['catted_outs'][wildcards.type+'_'+wildcards.subtype])

def get_R2_files(wildcards):
	return expand('input_files/fastq_files/'+wildcards.type+'_'+wildcards.subtype+'_{rep}_2.fastq', rep=config['catted_outs'][wildcards.type+'_'+wildcards.subtype])

rule all:
	input:
		R1_files=expand('output_files/catted_reps/{catted_outs}/{catted_outs}_catted_R1.fastq', catted_outs=config['catted_outs']),
		R2_files=expand('output_files/catted_reps/{catted_outs}/{catted_outs}_catted_R2.fastq', catted_outs=config['catted_outs'])

rule link_seeds:
	input:
		input_seeds='../construct_TE_seeds/output_files/dm6_TE_seeds.fa'
	params:
		relative_path='../../../construct_TE_seeds/output_files/dm6_TE_seeds.fa'
	output:
		link_name='input_files/input_fastas/dm6_TE_seeds.fa'
	shell:
		'ln -s -f -T {params.relative_path} {output.link_name}'

rule link_CNS_gene_mapping:
	input:
		input_CNS='../chimeric_copia_reads/output_files/refseq_reads/CNS_gene_mapping.tsv'
	params:
		relative_path='../../../chimeric_copia_reads/output_files/refseq_reads/CNS_gene_mapping.tsv'
	output:
		link_name='input_files/input_fastas/CNS_gene_mapping.tsv'
	shell:
		'ln -s -f -T {params.relative_path} {output.link_name}'

rule cat_reps:
	input:
		R1_files=get_R1_files,
		R2_files=get_R2_files
	output:
		R1='output_files/catted_reps/{type}_{subtype}/{type}_{subtype}_catted_R1.fastq',
		R2='output_files/catted_reps/{type}_{subtype}/{type}_{subtype}_catted_R2.fastq'
	shell:
		'''
		cat {input.R1_files} >{output.R1}
		cat {input.R2_files} >{output.R2}
		'''
