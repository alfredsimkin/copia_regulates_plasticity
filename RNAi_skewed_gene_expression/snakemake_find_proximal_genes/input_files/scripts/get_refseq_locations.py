'''
gets coordinates of anchor genes (a small list of genes of interest) and of all
fasta genes and outputs them to coordinate files
'''

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

refseq_fasta=snakemake.input['refseq_fasta']
refseq_coords=snakemake.output['refseq_coords']

def get_refseq_locations(fasta_name):
	fasta_list=read_fasta(fasta_name)
	titles=[gene[0] for gene in fasta_list]
	location_list=[]
	for title in titles:
		name='_'.join(title.split(' ')[0].split('_')[-2:])
		coords=title.split(' ')[1].split('range=')[1]
		chrom, endpoints=coords.split(':')
		start, end=endpoints.split('-')
		location_list.append([chrom, start, end, name])
	return location_list

def print_locations(location_list, output_name):
	output_file=open(output_name, 'w')
	for location_line in location_list:
		output_file.write('\t'.join(location_line)+'\n')

refseq_locs=get_refseq_locations(refseq_fasta)
print_locations(refseq_locs, refseq_coords)
