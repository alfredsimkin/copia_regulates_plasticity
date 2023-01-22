'''
Asks following questions:

1. How many nucleotides of Copia are in the genome?
2. How many distinct regions of the genome map to Copia?
'''

rule all:
	input:
		'output_files/summarized_genomic_copia.tsv'

rule find_copia_nucs:
	'''
	blats copia consensus sequence against genome
	'''
	input:
		query='input_files/input_refs/copia_consensus.fa',
		database='input_files/input_refs/dm6.fa'
	output:
		copia_nucs='output_files/all_copia_in_genome.psl'
	shell:
		'blat {input.database} {input.query} {output.copia_nucs}'

rule get_total_nuc_stats:
	'''
	calculates how many copia nucleotides are in the genome and from how many
	distinct regions
	'''
	input:
		genomic_copia_psl='output_files/all_copia_in_genome.psl'
	output:
		converted_coords='output_files/genomic_coords.txt',
		summary='output_files/summarized_genomic_copia.tsv',
		region_stats='output_files/region_stats.tsv'
	script:
		'input_files/scripts/summarize_genomic_copia.py'
