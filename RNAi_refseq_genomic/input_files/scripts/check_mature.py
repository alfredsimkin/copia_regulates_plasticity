input_psl=snakemake.input['psl_file']
output_summary=open(snakemake.output['status_summary'], 'w')
hits=0
for line in open(input_psl):
	line=line.strip().split()
	if len(line)>3 and line[0].isdigit():
		hits+=1
if hits>0:
	print('**********mature mRNAs found containing RNAi hits********')
	exit()
else:
	print('no mature mRNAs found containing any RNAi hits')
output_summary.write(f'mature mRNA genes that have RNAi hits:\t{hits}\n')
