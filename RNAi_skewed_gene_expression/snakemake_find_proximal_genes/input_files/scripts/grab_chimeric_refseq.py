'''
grabbing Refseq transcripts from tsv file and pulling matching genes from refseq
fasta file.
'''

input_tsv=snakemake.input['tsv_file']
input_fasta=snakemake.input['fasta_file']

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


def read_fasta_generator(fasta_name):
	'''
	Works like read_fasta but works as a generator instead of as a list. This
	allows large fasta files to be read in line by line without creating a large
	memory usage. Copied from biopython via a stack overflow forum.
	'''
	import gzip
	if fasta_name.endswith('.gz'):
		fp=gzip.open(fasta_name, mode='rt')
	else:
		fp=open(fasta_name)
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line[1:], []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

def slurp_seqs(file_name_list):
	'''
	the goal of this program is to take a list of sequence files (could be fasta
	files or fastq files, gzipped or not gzipped) and to yield a single nested
	list suitable for sending to fastq or fasta format or for iterating through.
	'''
	for file_name in file_name_list:
		altered_name=file_name
		print(altered_name)
		if file_name.endswith('.gz'):
			altered_name=file_name[:-3]
		if altered_name.endswith('.fasta') or altered_name.endswith('.fa') or altered_name.endswith('.fst'):
			new_list=read_fasta_generator(file_name)
		for paired_seq in new_list:
			yield paired_seq

output_fasta=snakemake.output['chimeric_refseq']

matching_seqs=set([])
for line in open(input_tsv):
	line=line.strip().split()
	if len(line)>1 and 'NM' in line[1]:
		NM_names=set(line[1].split(','))
		matching_seqs=matching_seqs|NM_names

seqs=slurp_seqs([input_fasta])
good_list=[]
for name, seq in seqs:
	gene_name=name.split(' ')[0]
	if gene_name in matching_seqs:
		good_list.append([name, seq])
print_fasta(good_list, output_fasta)
