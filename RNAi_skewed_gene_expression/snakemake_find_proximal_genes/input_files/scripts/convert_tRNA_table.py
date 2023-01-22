input_name=str(snakemake.input)
output_name=str(snakemake.output)

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

paired_list=[]
for gene_number, gene in enumerate(open(input_name)):
	if gene_number>0:
		gene=gene.strip().split()
		gene_name=gene[0]
		gene_coords=gene[5]
		gene_chrom, gene_coords=gene_coords.replace(',','').split(':')
		gene_start, gene_end=map(int, gene_coords.split('..'))
		gene_seq=gene[-1]
		new_title='_'.join(map(str, [gene_name, 'chr'+gene_chrom, gene_start, gene_end]))
		paired_list.append([new_title, gene_seq])

print_fasta(paired_list, output_name)
