'''
gets coordinates of anchor genes (a small list of genes of interest) and of all
fasta genes and outputs them to coordinate files
'''

anchor_fasta=snakemake.input['anchor_fasta']
anchor_coords=snakemake.output['anchor_coords']

def read_fasta(fasta_file):
	seq, name_list, seq_list, seq_dict='', [], [], {}
	for line in open(fasta_file):
		line=line.strip()
		if '>' in line:
			name_list.append(line[1:])
			if len(seq)>0:
				seq_list.append(seq)
				seq=''
		else:
			seq=seq+line
	seq_list.append(seq)
#	for seq_number, name in enumerate(name_list):
#		seq_dict[name]=seq_list[seq_number]
	return [[name, seq_list[name_number]] for name_number, name in enumerate(name_list)]


def get_anchor_locations(fasta_name):
	fasta_list=read_fasta(fasta_name)
	titles=[gene[0] for gene in fasta_list]
	final_titles=[]
	location_list=[]
	for title in titles:
		final_titles.append(title)
	for title in final_titles:
		chrom, start, end=title.split('_')[1:4]
		location_list.append([chrom, start, end, title])
	return location_list

def print_locations(location_list, output_name):
	output_file=open(output_name, 'w')
	for location_line in location_list:
		output_file.write('\t'.join(location_line)+'\n')

anchor_locs=get_anchor_locations(anchor_fasta)
print_locations(anchor_locs, anchor_coords)
