configfile: 'config2.yaml'

module prepare_reads:
	snakefile: "snakemake_prepare_reads/prepare_reads.smk"
	config: config
use rule * from prepare_reads as prepare_reads_*

module find_proximal:
	snakefile: 'snakemake_find_proximal_genes/find_proximal.smk'
	config: config
use rule * from find_proximal as find_proximal_*

module run_salmon:
	snakefile: 'snakemake_run_salmon/run_salmon.smk'
	config: config
use rule * from run_salmon as run_salmon_*

module analyze_salmon:
	snakefile: 'snakemake_analyze_salmon/analyze_salmon.smk'
	config: config
use rule * from analyze_salmon as analyze_salmon_*

rule all:
	input:
		rules.prepare_reads_all.input,
		rules.find_proximal_all.input,
		rules.run_salmon_all.input,
		rules.analyze_salmon_all.input
	default_target: True
