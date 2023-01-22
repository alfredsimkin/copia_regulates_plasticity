input_alignment=snakemake.input['alignment']
input_reads=snakemake.input['reads']
output_reads=snakemake.output['reads']
search_term=snakemake.params['search_term']
output_matches=open(snakemake.output['database_matches'], 'w')

print('input alignment is', input_alignment)

score_c, read_c, database_c=9,0,5

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
		if altered_name.endswith('.fasta') or altered_name.endswith('.fa') or altered_name.endswith('.fst'):
			new_list=read_fasta_generator(file_name)
		for paired_seq in new_list:
			yield paired_seq

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

matching_reads={}
for line in open(input_alignment):
	line=line.strip().split()
	if len(line)>3 and line[score_c].isdigit() and search_term in line[database_c]:
		read_name=line[read_c]
		if int(line[score_c])>100:
			matching_reads.setdefault(read_name, []).append(line[database_c])

seqs=slurp_seqs([input_reads])
good_list=[]
for name, seq in seqs:
	read_name=name.split(' ')[0]
	if read_name in matching_reads:
		good_list.append([name, seq])
print_fasta(good_list, output_reads)

for read in matching_reads:
	output_matches.write(read+'\t'+','.join(matching_reads[read])+'\n')
