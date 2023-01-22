'''
finds refseq genes that are within 'distance' units of copia genes
'''

anchor_coords=snakemake.input['anchor_coords']
refseq_coords=snakemake.input['refseq_coords']
refseq_fasta=snakemake.input['refseq_fasta']
#distance=int(snakemake.params['distance'])
distance=int(snakemake.wildcards['distance'])
output=str(snakemake.output)

def print_fasta(fasta_list, outfile, mode='w', line_chars=60):
	"""
	this program prints fasta paired lists to fasta format
	"""
	output=open(outfile, mode)
	for sequence in fasta_list:
		output.write(">"+sequence[0]+"\n")
		for char_number, char in enumerate(sequence[1]):
			output.write(char)
			if char_number%line_chars==line_chars-1 or char_number==(len(sequence[1]))-1:
				output.write("\n")
	output.close()

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

def make_dict(coord_file):
	location_dict={}
	for coord in open(coord_file):
		chrom, start, end, name=coord.strip().split('\t')
		start, end=int(start), int(end)
		if chrom not in location_dict:
			location_dict[chrom]=[]
		location_dict[chrom].append([start, end, name])
	for chrom in location_dict:
		location_dict[chrom]=sorted(location_dict[chrom])
	return location_dict

def run_comparison(anchor_dict, refseq_dict, distance):
	overlappers=[]
	for chrom in anchor_dict:
		for anchor in anchor_dict[chrom]:
			anchor_start, anchor_end, anchor_name=anchor
			if chrom in refseq_dict:
				for gene_number, gene in enumerate(refseq_dict[chrom]):
					gene_start, gene_end, gene_name=gene
					if abs(gene_start-anchor_end)<distance or abs(gene_end-anchor_start)<distance:
						overlappers.append([chrom]+gene)
	return overlappers

def extract_overlappers(overlappers, refseq_fasta, output):
	gene_set=set([])
	fasta_list=[]
	refseq_dict=dict(read_fasta(refseq_fasta))
	for overlapper in overlappers:
		overlapper_chrom, overlapper_start, overlapper_end, overlapper_gene=overlapper
		for gene_name in refseq_dict:
			exact_name='_'.join(gene_name.split(' ')[0].split('_')[-2:])
#			print(overlapper_gene, exact_name)
			if overlapper_gene==exact_name:
				if gene_name not in gene_set:
					fasta_list.append([gene_name, refseq_dict[gene_name]])
				else:
					print(overlapper_gene, 'is duplicated!!')
				gene_set.add(gene_name)
	print_fasta(fasta_list, output)
				
anchor_dict=make_dict(anchor_coords)
refseq_dict=make_dict(refseq_coords)
overlappers=run_comparison(anchor_dict, refseq_dict, distance)
print(len(overlappers))
extract_overlappers(overlappers, refseq_fasta, output)
