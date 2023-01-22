copia_reads=snakemake.input['copia_reads']
gene_reads=snakemake.input['gene_reads']
RNAi_reads=snakemake.input['RNAi_reads']
non_seed_RNAi_reads=snakemake.input['non_seed_RNAi_reads']
seed_RNAi_reads=snakemake.input['seed_RNAi_reads']
outfile=open(snakemake.output['final_summary'], 'w')

def write_reads(message, dict_type):
	outfile.write(f'{message}\n')
	for read in dict_type:
		outfile.write(f'{read}\t{dict_type[read]}\n')

copia_dict=dict([line.strip().split() for line in open(copia_reads)])
gene_dict=dict([line.strip().split() for line in open(gene_reads)])
RNAi_dict=dict([line.strip().split() for line in open(RNAi_reads)])
non_seed_RNAi_dict=dict([line.strip().split() for line in open(non_seed_RNAi_reads)])
seed_RNAi_dict=dict([line.strip().split() for line in open(seed_RNAi_reads)])

outfile.write(f'copia reads: {len(copia_dict)}\n')
outfile.write(f'copia reads that map to refseq genes: {len(gene_dict)}\n')
outfile.write(f'copia reads that map to refseq genes and RNAi regions: {len(RNAi_dict)}\n')
outfile.write(f'copia reads that do not map to known seeds but do map to RNAi regions: {len(non_seed_RNAi_dict)}\n')
outfile.write(f'copia reads that map to known seeds and map to RNAi regions: {len(seed_RNAi_dict)}\n')

write_reads('copia reads are below:', copia_dict)
write_reads('copia reads that also map to refseq genes are below:', gene_dict)
write_reads('copia reads that also map to refseq genes and RNAi regions are below:', RNAi_dict)
write_reads('copia reads that do not map to known Copia seeds and do map to RNAi regions are below:', non_seed_RNAi_dict)
write_reads('copia reads that map to known Copia seeds and map to RNAi regions are below:', seed_RNAi_dict)
