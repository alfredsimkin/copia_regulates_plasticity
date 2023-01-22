'''
goal is to create files in output_files/proximal_genes

Some features have goal of expression levels of neighboring reefseq genes:
	coord file of features->overlap with neighboring genes->proximal refseq fa
	TE_sourced_FW
	TE_sourced_gypsy_all
	TE_sourced_copia
	converted_tRNAs
	
Other features have goal of checking expression levels of exact transcripts:
	original_features->feature_fa->proximal_feature fa (poorly labeled as refseq)
	chimeric_refseqs
	dm6_refseq
'''
#configfile: 'config2.yaml'

ruleorder: get_refseq_proximal > overlap_genes
ruleorder: get_tRNA_locations > get_supplemental_locations

rule all:
	input:
		proximals=expand('output_files/proximal_genes/{distance}/proximal_refseqs_to_{anchor}.fa', distance=config['distances'], anchor=config['anchor_genes'])

rule convert_tRNA_table:
	input:
		'input_files/tRNA_table.txt'
	output:
		'output_files/converted_tRNAs.fa'
	script:
		'input_files/scripts/convert_tRNA_table.py'

rule grab_chimeric:
	input:
		fasta_file='input_files/input_fastas/dm6_genes.fa',
		tsv_file='input_files/input_fastas/CNS_gene_mapping.tsv'
	output:
		chimeric_refseq='output_files/proximal_genes/{distance}/proximal_refseqs_to_chimeric_refseqs.fa'
	script:
		'input_files/scripts/grab_chimeric_refseq.py'		

rule get_refseq_proximal:
	input:
		refseq_fasta='input_files/input_fastas/dm6_genes.fa'
	params:
		relative_location='../../../input_files/input_fastas/dm6_genes.fa'
	output:
		proximal_fasta='output_files/proximal_genes/{distance}/proximal_refseqs_to_dm6_refseq.fa'
	shell:
		'ln -s -T {params.relative_location} {output.proximal_fasta}'

rule get_refseq_locations:
	input:
		refseq_fasta='input_files/input_fastas/dm6_genes.fa'
	output:
		refseq_coords='output_files/coord_files/dm6_coords.txt'
	script:
		'input_files/scripts/get_refseq_locations.py'

rule get_tRNA_locations:
	input:
		anchor_fasta='output_files/converted_tRNAs.fa',
	output:
		anchor_coords='output_files/coord_files/converted_tRNAs_coords.txt',
	script:
		'input_files/scripts/get_anchor_locations.py'

rule get_supplemental_locations:
	'''
	similar purpose to get_tRNA_locations but works on input TEs with slightly
	different format and from a larger input fasta file that needs filtering.
	'''
	input:
		all_seeds='input_files/input_fastas/dm6_TE_seeds.fa',
		input_yaml='input_files/supplemental_coords.yaml'
	output:
		coord_file='output_files/coord_files/TE_sourced_{anchor}_coords.txt'
	script:
		'input_files/scripts/get_supplemental_locations.py'

rule overlap_genes:
	input:
		anchor_coords='output_files/coord_files/{anchor}_coords.txt',
		refseq_coords='output_files/coord_files/dm6_coords.txt',
		refseq_fasta='input_files/input_fastas/dm6_genes.fa'
	output:
		proximals='output_files/proximal_genes/{distance}/proximal_refseqs_to_{anchor}.fa'
	script:
		'input_files/scripts/overlap_genes.py'
