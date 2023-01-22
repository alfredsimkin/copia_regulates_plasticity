import plotly.express as px
import pandas as pd
summary_path=str(snakemake.input['summary_file'])
graph_path=snakemake.output['plotly_graph']

def populate_dict(summary_path):
	chart_dict={'dataset':[], 'fraction':[], 'direction':[]}
	for line_number, line in enumerate(open(summary_path)):
		line=line.strip().split('\t')
		if line_number>0 and len(line)>2:			
			RNAi=line[1].split('RNAi_copia')[1]
			anchor='_'.join(line[2].split('_')[:-1])
			up=int(line[7])
			down=int(line[8])
			chart_dict['dataset'].append(f'{RNAi}_{anchor}')
			chart_dict['dataset'].append(f'{RNAi}_{anchor}')
			chart_dict['fraction'].append(down/(up+down))
			chart_dict['direction'].append('down')
			chart_dict['fraction'].append(up/(up+down))
			chart_dict['direction'].append('up')
	return chart_dict

chart_dict=populate_dict(summary_path)
chart_df=pd.DataFrame(chart_dict)
fig = px.bar(chart_df, x="dataset", y="fraction", color="direction", title="genes with > 2 fold change in expression after RNAi")
fig.write_html(str(graph_path)+'.html')
fig.write_image(str(graph_path))
