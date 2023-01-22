'''
examines 3 input files of genomic coordinates:
RNAi matches in genome
Refseq regions of genome
copia seed regions in genome

question is how many RNAi regions are in refseq genes, how many RNAi regions are
in copia seeds, and how many RNAi regions are in both
'''

copia_list=snakemake.input['copia_list']
RNAi_list=snakemake.input['RNAi_list']
refseq_list=snakemake.input['refseq_list']
outfile=open(snakemake.output['containment_stats'], 'w')

def process_coords(input_list):
	coords=[]
	for entry in open(input_list):
		chrom, start, end, descriptor=entry.strip().split()
		start, end=list(map(int, [start, end]))
		coords.append([chrom, start, end, descriptor])
	coords.sort()
	return coords

copia_processed=process_coords(copia_list)
RNAi_processed=process_coords(RNAi_list)
refseq_processed=process_coords(refseq_list)
in_both,refseq_only,copia_only,in_neither=0,0,0,0
off_targets=[]
for RNAi_hit in RNAi_processed:
	RNAi_chrom, RNAi_start, RNAi_end, RNAi_size=RNAi_hit
	in_copia=False
	for copia_seed in copia_processed:
		copia_chrom, copia_start, copia_end, copia_name=copia_seed
		if RNAi_chrom==copia_chrom and copia_start<RNAi_start and copia_end>RNAi_end:
			in_copia=True
			break
	in_refseq=False
	for refseq_region in refseq_processed:
		refseq_chrom, refseq_start, refseq_end, refseq_size=refseq_region
		if RNAi_chrom==refseq_chrom and refseq_start<RNAi_start and refseq_end>RNAi_end:
			in_refseq=True
			break
	if in_refseq and in_copia:
		in_both+=1
	elif in_refseq and not in_copia:
		off_targets.append(f'{refseq_region}\t{RNAi_hit}')
		print(refseq_region, RNAi_hit)
		refseq_only+=1
	elif in_copia and not in_refseq:
		copia_only+=1
	else:
		in_neither+=1
outfile.write(f'total RNAi hits in genome:\t{len(RNAi_processed)}\n')
outfile.write(f'in copia seed and in immature mRNA:\t{in_both}\n')
outfile.write(f'in copia seed but not in immature mRNA:\t{copia_only}\n')
outfile.write(f'in immature mRNA but not in copia seed:\t{refseq_only}\n')
outfile.write(f'not in copia seed or in immature mRNA:\t{in_neither}\n')
outfile.write(f'regions in immature mRNA missing from Copia:\n')
outfile.write(f'refseq_region\tRNAi_hit\n')
for off_target in off_targets:
	outfile.write(f'{off_target}\n')
