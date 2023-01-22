infiles=list(snakemake.input['input_files'])
outname=str(snakemake.output['tsv_file']).strip()
outfile=open(outname, 'w')
outfile.write(f'control\texperimental\tanchors\tdec_thresh\tinc_thresh\tcommon\tunchanged\tup\tdown\tP(down>=obs)\n')
print(infiles)
for infile in infiles:
	values=[]
	for line in open(infile):
		if ':' in line:
			values.append(line.strip().split(': ')[1])
	file_name=infile.split('/')[-1]
	control, other=file_name.replace('.txt', '').split('_vs_')
	exp, anchor=other.split('_in_')
	outfile.write('\t'.join([control, exp, anchor]+values)+'\n')
