'''
takes RNAi hits from the genome and copia seed coordinates and parses them into
lists of the format - chrom start end size. These lists will later (in another
script) be compared with each other and with refseq gene coordinates (computed
elsewhere) to determine how many RNAi matches occur inside genes, how many RNAi
matches occur inside copia seeds, and how many RNAi matches occur inside both.
'''

copia_fasta=snakemake.input['copia_seeds']
RNAi_psl=snakemake.input['RNAi_hits']
RNAi_list=snakemake.output['RNAi_list']
copia_list=snakemake.output['copia_list']

def parse_copia(copia_fasta, copia_list):
	outfile=open(copia_list, 'w')
	for line in open(copia_fasta):
		if 'COPIA_DM' in line:
			name=line[1:-1]
			gene=line.strip().split('_')
			chrom, start, end=gene[3:6]
			start, end=int(start), int(end)
			outfile.write(f'{chrom}\t{start}\t{end}\t{name}\n')
	outfile.close()

def parse_RNAi(RNAi_psl, RNAi_list):
	outfile=open(RNAi_list, 'w')
	for line in open(RNAi_psl):
		line=line.strip().split()
		if len(line)>3 and line[0].isdigit():
			chrom=line[13]
			start=int(line[15])
			end=int(line[16])
			name=line[9]
			outfile.write(f'{chrom}\t{start}\t{end}\t{name}\n')
	outfile.close()

parse_copia(copia_fasta, copia_list)
parse_RNAi(RNAi_psl, RNAi_list)
