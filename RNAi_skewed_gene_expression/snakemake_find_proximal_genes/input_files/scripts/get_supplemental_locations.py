import yaml
yaml_dict=yaml.load(open(snakemake.input['input_yaml']), Loader=yaml.SafeLoader)
all_seeds=snakemake.input['all_seeds']

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

seed_names=[item[0] for item in (read_fasta(all_seeds))]
coord_file=open(snakemake.output['coord_file'], 'w')
anchor=snakemake.wildcards['anchor']

search_term=yaml_dict['input_anchors'][anchor]
for seed in seed_names:
	if search_term in seed:
		name, location=seed.split('_chr')
		location_list=('chr'+location).split('_')
		start, end, strand=location_list[-3:]
		chrom='_'.join(location_list[:-3])
		coord_file.write('\t'.join([chrom,start,end, seed])+'\n')
