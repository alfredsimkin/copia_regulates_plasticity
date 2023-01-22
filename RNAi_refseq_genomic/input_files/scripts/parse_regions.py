'''
much of the code for this is borrowed from summarize_genomic_copia.py
Modified slightly to use coordinate ranges from a fasta file instead of from a
psl file, and to treat '+' and '-' strands of the genome as separate from each
other.
'''

input_fasta=snakemake.input['gene_file']
summary_out=open(snakemake.output['stats_file'], 'w') 
converted_coords=open(snakemake.output['detailed_regions'], 'w')
region_stats=open(snakemake.output['refseq_list'], 'w')

def make_lists(input_fasta):
	coords=[]
	for line in open(input_fasta):
		if line.startswith('>'):
			line=line.strip().split(' ')
			region=line[1].split('=')[1]
			chrom, span=region.split(':')
			start, end=list(map(int, span.split('-')))
			strand=line[4].split('=')[1]
			coords.append([chrom, strand, start, end])
	coords.sort()
	return coords

coords=make_lists(input_fasta)
previous_chrom, previous_strand, max_end='junk', 'junk', 0
region_count,nuc_count=0,0
for entry in coords:
	chrom, strand, start, end=entry
	print(chrom, start, end, max_end)
	if chrom!=previous_chrom or strand!=previous_strand or start>max_end+100:
		if previous_chrom!='junk':
			region_stats.write(f'{previous_chrom}\t{region_start}\t{max_end}\t{max_end-region_start}\n')
		region_count+=1
		region_start=start
	overlapping=False
	if chrom==previous_chrom and strand==previous_strand and start<max_end:
		overlapping=True
	if overlapping and end>max_end:
		nuc_count+=(end-max_end)
	elif not overlapping:
		nuc_count+=(end-start)
	converted_coords.write(f'{chrom}\t{start}\t{end}\t{region_count}\t{nuc_count}\n')
	if (chrom!=previous_chrom or strand!=previous_strand) or end>max_end:
		max_end=end
	previous_chrom, previous_strand=chrom, strand
region_stats.write(f'{previous_chrom}\t{region_start}\t{max_end}\t{max_end-region_start}\n')

summary_out.write(f'total_regions:\t{region_count}\n')
summary_out.write(f'total_nucs:\t{nuc_count}\n')
