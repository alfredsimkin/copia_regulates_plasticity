rule all:
	input:
		off_target_folder=directory('output_files/off-target/unaligned'),
		status_summary='output_files/summaries/dm6_genes_containment_stats.tsv'

rule link_seeds:
	input:
		input_seeds='../construct_TE_seeds/output_files/dm6_TE_seeds.fa'
	params:
		relative_path='../../../construct_TE_seeds/output_files/dm6_TE_seeds.fa'
	output:
		link_name='input_files/fasta_files/dm6_TE_seeds.fa'
	shell:
		'ln -s -f -T {params.relative_path} {output.link_name}'

rule get_regions:
	'''
	serves two functions:
	1. combines overlapping transcripts into non-redundant regions
	2. outputs non-redundant coordinate ranges for assessing RNAi overlap with
	genes - format is chrom, start, end, size. These are not fully sorted.
	
	Question: should I really be separating out '+' and '-'? Arguably, if an
	RNAi hit overlaps both '+' and '-' regions that overlap each other, this
	should be counted once, not twice.
	'''
	input:
		gene_file='input_files/fasta_files/dm6_genes_plus_introns.fa'
	output:
		stats_file='output_files/summaries/dm6_genes_plus_introns_region_stats.tsv',
		detailed_regions='output_files/summaries/dm6_genes_plus_introns_region_stats_detailed.tsv',
		refseq_list='output_files/coord_lists/dm6_genes_plus_introns_refseq_coords.tsv'
	script:
		'input_files/scripts/parse_regions.py'

rule align_RNAi:
	input:
		database='input_files/fasta_files/{ref}.fa',
		query='input_files/fasta_files/RNAi.fa'
	output:
		psl_file='output_files/psl_files/RNAi_in_{ref}.psl'
	shell:
		'blat {input.database} {input.query} {output.psl_file}'

rule check_mature:
	'''
	simple script to check whether there are any hits to mature mRNAs in the
	RNAi queries
	'''
	input:
		psl_file='output_files/psl_files/RNAi_in_dm6_genes.psl'
	output:
		status_summary='output_files/summaries/dm6_genes_containment_stats.tsv'
	script:
		'input_files/scripts/check_mature.py'

rule parse_genomic_lists:
	'''
	converts RNAi regions (psl) and Copia regions (fasta) into tabular lists
	'''
	input:
		copia_seeds='input_files/fasta_files/dm6_TE_seeds.fa',
		RNAi_hits='output_files/psl_files/RNAi_in_dm6.psl'
	output:
		copia_list='output_files/coord_lists/copia_coords.tsv',
		RNAi_list='output_files/coord_lists/RNAi_coords.tsv'
	script:
		'input_files/scripts/parse_lists.py'

rule check_containment:
	'''
	checks if RNAi region is inside a gene region, inside a copia region, inside
	both, or inside neither.
	'''
	input:
		copia_list='output_files/coord_lists/copia_coords.tsv',
		RNAi_list='output_files/coord_lists/RNAi_coords.tsv',
		refseq_list='output_files/coord_lists/dm6_genes_plus_introns_refseq_coords.tsv'
	output:
		containment_stats='output_files/summaries/dm6_genes_plus_introns_containment_stats.tsv'
	script:
		'input_files/scripts/check_containments.py'

rule check_off_target:
	'''
	For any 'off target' regions found by 'check_containment' (inside gene
	region but not copia seed), retrieves all transcripts that overlap.
	'''
	input:
		refseq_file='input_files/fasta_files/dm6_genes_plus_introns.fa',
		containment_stats='output_files/summaries/dm6_genes_plus_introns_containment_stats.tsv',
		consensus_seqs='input_files/fasta_files/TE_consensus.fa'
	output:
		off_target_folder=directory('output_files/off-target/unaligned')
	script:
		'input_files/scripts/off_target_transcripts.py'
