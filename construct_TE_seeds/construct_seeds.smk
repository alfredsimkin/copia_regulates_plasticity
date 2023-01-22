rule all:
	input: final_file='output_files/dm6_TE_seeds.fa'
		
rule blat_genome:
	input:
		dm6_genome='input_files/fasta_files/dm6.fa',
		TE_consensus='input_files/fasta_files/TE_consensus.fa'
	output:
		genomic_TE_psl='output_files/dm6_TE_hits.pslx'
	shell:
		'blat {input.dm6_genome} {input.TE_consensus} -out=pslx {output.genomic_TE_psl}'

rule reconstruct_seeds:
	input:
		blat_file='output_files/dm6_TE_hits.pslx',
		dm6_genome='input_files/fasta_files/dm6.fa'
	output:
		all_TE_seeds='output_files/dm6_TE_seeds.fa'
	script:
		'input_files/scripts/reconstruct_seeds.py'
