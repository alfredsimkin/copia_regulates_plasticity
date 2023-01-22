rule all:
	input:
		plotly_graph=expand('output_files/summarized_final_stats/{distance}/{distance}_{smaller}-{bigger}_stacked_bar.pdf', smaller=config['smaller'], bigger=config['bigger'], distance=config['distances']),

rule compare_pairs:
	'''
	input smaller number first. E.g. to see genes that decrease to 1/2 or
	increase to 2/1, input 1 as smaller number and 2 as bigger number
	'''
	input:
		control='output_files/salmon_output/{distance}/{control}_in_{anchor}/quant.sf',
		experimental='output_files/salmon_output/{distance}/{experimental}_in_{anchor}/quant.sf'
	params:
		fold_change_smaller=config['smaller'],
		fold_change_bigger=config['bigger']
	output:
		stat_file='output_files/final_stats/{distance}/{control}_vs_{experimental}_in_{anchor}_'+str(config['smaller'])+'-'+str(config['bigger'])+'.txt'
	script:
		'input_files/scripts/compare_pairs.py'

rule summarize_comparisons:
	'''
	take text file info from all pair comparisons and tabulate it
	'''
	input:
		input_files=expand('output_files/final_stats/{distance}/{control}_vs_{experimental}_in_{anchor}_{smaller}-{bigger}.txt', control=config['controls'], experimental=config['experimentals'], anchor=config['anchor_genes'], smaller=config['smaller'], bigger=config['bigger'], allow_missing=True)
	output:
		tsv_file=expand('output_files/summarized_final_stats/{distance}/{distance}_{smaller}-{bigger}_summary.tsv', smaller=config['smaller'], bigger=config['bigger'], allow_missing=True)
	script:
		'input_files/scripts/summarize_stats.py'

rule graph_comparisons:
	input:
		summary_file=expand('output_files/summarized_final_stats/{distance}/{distance}_{smaller}-{bigger}_summary.tsv', smaller=config['smaller'], bigger=config['bigger'], allow_missing=True)
	output:
		plotly_graph=expand('output_files/summarized_final_stats/{distance}/{distance}_{smaller}-{bigger}_stacked_bar.pdf', smaller=config['smaller'], bigger=config['bigger'], allow_missing=True)
	script:
		'input_files/scripts/graph_stats.py'
