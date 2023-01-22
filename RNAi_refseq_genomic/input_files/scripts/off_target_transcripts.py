import subprocess

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

def revcom(seq,nuc='DNA'):
	if nuc=='DNA':
		complement={'N':'N','n':'n','A':'T','C':'G','G':'C','T':'A','a':'t','t':'a','c':'g','g':'c', 'U':'A', 'u':'a', '-':'-'}
	else:
		complement={'N':'N','n':'n','A':'U','C':'G','G':'C','U':'A','a':'u','u':'a','c':'g','g':'c','-':'-'}
	return ''.join(reversed([complement[base] for base in seq]))

refseq_file=snakemake.input['refseq_file']
containment_stats=snakemake.input['containment_stats']
off_target_folder=snakemake.output['off_target_folder']
consensus_seqs=snakemake.input['consensus_seqs']
subprocess.call(['mkdir', off_target_folder])
consensus_seqs=read_fasta(consensus_seqs)
copia_consensus=[[name, seq] for name, seq in consensus_seqs if 'COPIA' in name]
fasta_dict=dict(read_fasta(refseq_file))

regions=set([])
recording=False
for line in open(containment_stats):
	if recording:
		regions.add(tuple(eval(line.strip().split('\t')[0])))
	if line.startswith('refseq_region'):
		recording=True
for gene in fasta_dict:
	name=gene.split(' ')[0].split('RefSeq_')[1]
	location=gene.split(' ')[1].split('=')[1]
	t_chrom, coords=location.split(':')
	t_start, t_end=map(int, coords.split('-'))
	for region in regions:
		r_chrom, r_start, r_end=region[:3]
		if t_chrom==r_chrom and t_start<=r_end and t_end>=r_start:
			seq=fasta_dict[gene]
			paired_list=[[name+'_'+location, seq]]+copia_consensus
			print('paired is', paired_list[1])
			print_fasta(paired_list, off_target_folder+'/'+name+'.fa')
			paired_list=[[name+'_'+location+'_revcom', revcom(seq)]]+copia_consensus
			print_fasta(paired_list, off_target_folder+'/'+name+'_revcom.fa')
