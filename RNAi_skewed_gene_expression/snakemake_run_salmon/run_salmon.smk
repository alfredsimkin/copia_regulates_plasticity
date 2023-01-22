rule all:
	input:
		out_file=expand('output_files/salmon_output/{distance}/{sample}_in_{anchor}/quant.sf', distance=config['distances'], sample=config['catted_outs'], anchor=config['anchor_genes'])

rule index_genes:
	'''
	generates the salmon index of genes that are proximal to some anchor feature
	'''
	input:
		'output_files/proximal_genes/{distance}/proximal_refseqs_to_{anchor}.fa'
	output:
		directory('output_files/proximal_genes/{distance}/refseq_salmon_index_of_{anchor}')
	shell:
		'salmon index -t {input} -i {output}'

rule run_salmon:
	input:
		R1='output_files/catted_reps/{sample}/{sample}_catted_R1.fastq',
		R2='output_files/catted_reps/{sample}/{sample}_catted_R2.fastq',
		index='output_files/proximal_genes/{distance}/refseq_salmon_index_of_{anchor}'
	output:
		out_file='output_files/salmon_output/{distance}/{sample}_in_{anchor}/quant.sf'
	params:
		out_dir=directory('output_files/salmon_output/{distance}/{sample}_in_{anchor}')
	shell:
		'salmon quant -i {input.index} -l a --validateMappings -1 {input.R1} -2 {input.R2} -o {params.out_dir}'
