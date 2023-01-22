'''
this program extracts 'full length' sequences using blat hits and a genome as
inputs, where 'full length' is defined as a hit that starts at the beginning of
a query and ends at the end of the query, and has query gaps and subject gaps
that both are less than 1,000 letters (assumes genome and query are both
unspliced).
'''
def read_fasta(fasta_file):
	'''
	this function converts fasta files into a nested list of paired
	[name, sequence] lists
	'''
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

def print_fasta(fasta_list, outfile, mode='w', line_chars=60):
	"""
	this function prints fasta paired lists to fasta format
	"""
	output=open(outfile, mode)
	for sequence in fasta_list:
		output.write(">"+sequence[0]+"\n")
		for char_number, char in enumerate(sequence[1]):
			output.write(char)
			if char_number%line_chars==line_chars-1 or char_number==(len(sequence[1]))-1:
				output.write("\n")

blat_file=snakemake.input['blat_file']
genome_reference=snakemake.input['dm6_genome']
output_name=snakemake.output['all_TE_seeds']
genome_fasta=dict(read_fasta(genome_reference))

output_list=[]
for line_number, line in enumerate(open(blat_file)):
	if line_number>4:
		split_line=line.split()
		if int(split_line[11])==0 and int(split_line[12])==int(split_line[10]) and int(split_line[5])<1000 and int(split_line[7])<1000:
			TE_name, chrom, start, end, strand=split_line[9], split_line[13], split_line[15], split_line[16], split_line[8]
			strand='strand='+strand
			name='_'.join([TE_name, chrom, start, end, strand])
			seq=genome_fasta[chrom][int(start):int(end)].upper()
			output_list.append([name, seq])

print_fasta(output_list, output_name)
