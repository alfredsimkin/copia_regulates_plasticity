'''
aligns individual fastq files against dm6 using flair
'''

rule all:
	input:
		quantify='output_files/mapped_reads.counts.tsv'

rule quantify_abundances:
	input:
		manifest='input_files/reads_manifest.tsv',
		fasta='input_files/copia_isoforms.fa'
	params:
		prefix='output_files/mapped_reads'
	output:
		quantify='output_files/mapped_reads.counts.tsv'
	threads: 100
	shell:
		'''
		flair quantify -r {input.manifest} -i {input.fasta} -t {threads} -o {params.prefix}
		'''
