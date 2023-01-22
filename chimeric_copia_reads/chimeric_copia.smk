'''
converts input reads to fasta format and searches reads against copia (blat or
minimap2). searches these reads against the genome, and parses them into reads
that map only to copia and are not chimeric vs. reads that map to copia and
somewhere else (chimeras). Annotates where these 'somewhere else' reads map to.

- handle ties
- gypsy
- up-down of chimeric NM_numbers

Not essential but interesting:

Of the non-seed regions that copia mapping reads map to, how many are novel
relative to original 38? How many are novel relative to 85 copia-mapping regions
of the genome (from another snakemake pipeline)?
'''
configfile: 'input_files/config.yaml'

rule all:
	input:
		final_summary='output_files/CNS_final_summary.tsv'
		#alignment=expand('output_files/RNAi_alignments/CNS_{number}_RNAi_alignments.tsv', number=[0,1,2])

rule link_seeds:
	input:
		input_seeds='../construct_TE_seeds/output_files/dm6_TE_seeds.fa'
	params:
		relative_path='../../../construct_TE_seeds/output_files/dm6_TE_seeds.fa'
	output:
		link_name='input_files/reference_fastas/dm6_TE_seeds.fa'
	shell:
		'ln -s -f -T {params.relative_path} {output.link_name}'

rule convert_to_fasta:
	input:
		fastq='input_files/nanopore_fastq/CNS_{number}.fastq'
	output:
		fasta='output_files/nanopore_rep_fasta/CNS_{number}.fasta'
	run:
		fasta_file=open(output.fasta, 'w')
		for line_number, line in enumerate(open(input.fastq)):
			if line_number%4==0 and '@' in line:
				fasta_file.write('>'+line[1:])
			elif line_number%4==1:
				fasta_file.write(line)

rule cat_reps:
	input:
		original0='output_files/nanopore_rep_fasta/CNS_0.fasta',
		original1='output_files/nanopore_rep_fasta/CNS_1.fasta',
		original2='output_files/nanopore_rep_fasta/CNS_2.fasta'
	output:
		catted='output_files/nanopore_fasta/CNS.fasta'
	shell:
		'cat {input.original0} {input.original1} {input.original2} >{output.catted}'

rule search_copia:
	'''
	search for copia consensus sequence in long reads
	'''
	input:
		query='output_files/nanopore_fasta/CNS.fasta',
		database=config['ref_copia']
	output:
		alignment='output_files/copia_alignments/CNS_copia_alignments.tsv'
	shell:
		'minimap2 {input.database} {input.query} -o {output.alignment}'

rule grab_copia_reads:
	'''
	retrieve the subset of reads that map to copia
	'''
	input:
		reads='output_files/nanopore_fasta/CNS.fasta',
		alignment='output_files/copia_alignments/CNS_copia_alignments.tsv'
	params:
		search_term='COPIA_DM_I'
	output:
		reads='output_files/copia_reads/CNS_copia_mapping.fasta',
		database_matches='output_files/copia_reads/CNS_copia_mapping.tsv'
	script:
		'input_files/scripts/grab_search_term.py'

rule search_genes:
	'''
	search for the subset of copia mapping reads that map to Refseq genes
	'''
	input:
		query='output_files/copia_reads/CNS_copia_mapping.fasta',
		database=config['ref_genes']
	output:
		alignment='output_files/refseq_alignments/CNS_gene_alignments.tsv'
	shell:
		'minimap2 {input.database} {input.query} -o {output.alignment}'

rule grab_refseq_reads:
	'''
	retrieve the subset of copia mapping reads that map to Refseq sequences
	'''
	input:
		reads='output_files/copia_reads/CNS_copia_mapping.fasta',
		alignment='output_files/refseq_alignments/CNS_gene_alignments.tsv'
	params:
		search_term='dm6_ncbiRefSeq'
	output:
		reads='output_files/refseq_reads/CNS_gene_mapping.fasta',
		database_matches='output_files/refseq_reads/CNS_gene_mapping.tsv'
	script:
		'input_files/scripts/grab_search_term.py'

rule search_RNAi:
	'''
	search for the subset of Refseq and copia mapping reads that map to RNAi
	regions
	'''
	input:
		query='output_files/refseq_reads/CNS_gene_mapping.fasta',
		database=config['ref_RNAi']
	output:
		alignment='output_files/RNAi_alignments/CNS_RNAi_alignments.tsv'
	shell:
		'minimap2 {input.database} {input.query} -o {output.alignment}'

rule grab_RNAi_reads:
	'''
	retrieve the subset of Refseq and copia mapping reads that map to RNAi
	regions
	'''
	input:
		reads='output_files/refseq_reads/CNS_gene_mapping.fasta',
		alignment='output_files/RNAi_alignments/CNS_RNAi_alignments.tsv'
	params:
		search_term='Copia'
	output:
		reads='output_files/RNAi_reads/CNS_RNAi_mapping.fasta',
		database_matches='output_files/RNAi_reads/CNS_RNAi_mapping.tsv'
	script:
		'input_files/scripts/grab_search_term.py'

rule search_genome:
	'''
	search the genome for reads that map to copia
	'''
	input:
		query='output_files/copia_reads/CNS_copia_mapping.fasta',
		database=config['ref_genome']
	output:
		alignment='output_files/genome_alignments/CNS_genome_alignments.tsv'
	shell:
		'minimap2 {input.database} {input.query} -o {output.alignment}'

rule map_to_seeds:
	'''
	takes genomic mapped reads and checks which are best aligned to known copia
	seed locations vs. which are not.
	'''
	input:
		alignment='output_files/genome_alignments/CNS_genome_alignments.tsv',
		seeds='input_files/reference_fastas/dm6_TE_seeds.fa',
		reads='output_files/copia_reads/CNS_copia_mapping.fasta'
	output:
		seed_summary='output_files/seed_mapping/CNS_seed_results.tsv',
		non_seed_reads='output_files/seed_mapping/CNS_non_seed_reads.fasta',
		seed_reads='output_files/seed_mapping/CNS_seed_reads.fasta'
	script:
		'input_files/scripts/map_to_seeds.py'

rule search_RNAi_non_seeds:
	'''
	checks to see which non-seed mapping reads contain RNAi regions
	'''
	input:
		query='output_files/seed_mapping/CNS_non_seed_reads.fasta',
		database=config['ref_RNAi']
	output:
		alignment='output_files/non_seed_RNAi_alignments/CNS_non_seed_RNAi_alignments.tsv'
	shell:
		'minimap2 {input.database} {input.query} -o {output.alignment}'

rule search_RNAi_seeds:
	'''
	checks to see which seed mapping reads contain RNAi regions
	'''
	input:
		query='output_files/seed_mapping/CNS_seed_reads.fasta',
		database=config['ref_RNAi']
	output:
		alignment='output_files/seed_RNAi_alignments/CNS_seed_RNAi_alignments.tsv'
	shell:
		'minimap2 {input.database} {input.query} -o {output.alignment}'

rule grab_seed_RNAi_reads:
	input:
		reads='output_files/seed_mapping/CNS_seed_reads.fasta',
		alignment='output_files/seed_RNAi_alignments/CNS_seed_RNAi_alignments.tsv'
	params:
		search_term='Copia'
	output:
		reads='output_files/seed_RNAi_reads/CNS_seed_RNAi_mapping.fasta',
		database_matches='output_files/seed_RNAi_reads/CNS_seed_RNAi_mapping.tsv'
	script:
		'input_files/scripts/grab_search_term.py'

rule grab_non_seed_RNAi_reads:
	input:
		reads='output_files/seed_mapping/CNS_non_seed_reads.fasta',
		alignment='output_files/non_seed_RNAi_alignments/CNS_non_seed_RNAi_alignments.tsv'
	params:
		search_term='Copia'
	output:
		reads='output_files/non_seed_RNAi_reads/CNS_non_seed_RNAi_mapping.fasta',
		database_matches='output_files/non_seed_RNAi_reads/CNS_non_seed_RNAi_mapping.tsv'
	script:
		'input_files/scripts/grab_search_term.py'


rule summarize_output:
	input:
		copia_reads='output_files/copia_reads/CNS_copia_mapping.tsv',
		gene_reads='output_files/refseq_reads/CNS_gene_mapping.tsv',
		RNAi_reads='output_files/RNAi_reads/CNS_RNAi_mapping.tsv',
		non_seed_RNAi_reads='output_files/non_seed_RNAi_reads/CNS_non_seed_RNAi_mapping.tsv',
		seed_RNAi_reads='output_files/seed_RNAi_reads/CNS_seed_RNAi_mapping.tsv'
	output:
		final_summary='output_files/CNS_final_summary.tsv'
	script:
		'input_files/scripts/summarize_searches.py'
