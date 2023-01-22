input_alignment=snakemake.input['alignment']
seeds=snakemake.input['seeds']
input_reads=snakemake.input['reads']
non_seed_read_file=snakemake.output['non_seed_reads']
seed_read_file=snakemake.output['seed_reads']
final_summary=open(snakemake.output['seed_summary'], 'w')

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

def parse_seeds(seed_names):
	parsed_list=[]
	for seed in seed_names:
		chrom, start, end=seed.split('_')[-4:-1]
		parsed_list.append([chrom, int(start), int(end)])
	return parsed_list

def get_best_hits(input_alignment):
	score_c, read_c, database_c, start_c, end_c=9,0,5,7,8
	matching_reads={}
	high_scoring_reads=set([])
	for line in open(input_alignment):
		line=line.strip().split()
		if len(line)>3 and line[score_c].isdigit():
			read_name=line[read_c]
			chrom=line[database_c]
			start=int(line[start_c])
			end=int(line[end_c])
			score=int(line[score_c])
			if score>100:
				high_scoring_reads.add(read_name)
				if read_name not in matching_reads or score>matching_reads[read_name][0]:
					matching_reads[read_name]=[score, chrom, start, end]
	print('high scoring reads are', len(high_scoring_reads))
	return matching_reads

def place_reads(matching_reads, parsed_seeds):
	seed_reads, non_seed_reads=[],[]
	for read_name in matching_reads:
		score, read_chrom, read_start, read_end=matching_reads[read_name]
		in_seed=False
		for seed_chrom, seed_start, seed_end in parsed_seeds:
			if read_chrom==seed_chrom and seed_start<read_start and seed_end>read_end:
				in_seed=True
				break
		if in_seed:
			seed_reads.append(read_name)
		else:
			non_seed_reads.append(read_name)
	return seed_reads, non_seed_reads

def summarize_results(aligned_reads, seed_reads, non_seed_reads, final_summary):
	final_summary.write(f'total number of reads in seeds is:\t{len(seed_reads)}\n')
	final_summary.write(f'total number of reads not in seeds is:\t{len(non_seed_reads)}\n')
	final_summary.write(f'below are the genomic mappings of all reads mapped to seeds:\tscore\tchrom\tstart\tend\n')
	for read in seed_reads:
		score, chrom, start, end=aligned_reads[read]
		final_summary.write(f'{read}\t{score}\t{chrom}\t{start}\t{end}\n')
	final_summary.write(f'below are the genomic mappings of all reads that did not map to seeds:\n')
	for read in non_seed_reads:
		score, chrom, start, end=aligned_reads[read]
		final_summary.write(f'{read}\t{score}\t{chrom}\t{start}\t{end}\n')

def print_matchers(input_reads, matching_reads, output_file):
	seqs=slurp_seqs([input_reads])
	good_list=[]
	for name, seq in seqs:
		read_name=name.split(' ')[0]
		if read_name in matching_reads:
			good_list.append([name, seq])
	print_fasta(good_list, output_file)

seed_list=list(read_fasta_generator(seeds))
seed_names=[name for name, seq in seed_list if 'COPIA_DM' in name]
print('input is:', len(input_alignment))
aligned_reads=get_best_hits(input_alignment)
parsed_seeds=parse_seeds(seed_names)
print('aligned is:', len(aligned_reads))
print('parsed is:', len(parsed_seeds))
seed_reads, non_seed_reads=place_reads(aligned_reads, parsed_seeds)
summarize_results(aligned_reads, seed_reads, non_seed_reads, final_summary)
print_matchers(input_reads, non_seed_reads, non_seed_read_file)
print_matchers(input_reads, seed_reads, seed_read_file)
