'''
outputs two summaries:
1. how many distinct regions of the genome align to consensus copia sequence
2. how many nucleotides of the genome are captured in copia alignments

To accomplish this, blat alignments are split up by exon and converted to
genomic coordinates. Genomic coordinates are sorted by chromosome and by start
coordinate. Regions are defined as start coordinates>100 nucleotides after
previous end coordinate. Nucleotides are defined by subtracting end from start
(if nonoverlapping relative to previous entry) or end from previous end
(if overlapping)

'''

input_psl=snakemake.input['genomic_copia_psl']
summary_out=open(snakemake.output['summary'], 'w')
converted_coords=open(snakemake.output['converted_coords'], 'w')
region_stats=open(snakemake.output['region_stats'], 'w')
coords=[]
for line in open(input_psl):
	line=line.strip().split()
	if len(line)>3 and line[0].isdigit():
		tname=line[13]
#		print(line[18], line[20])
		sizes=list(map(int, line[18].strip(',').split(',')))
		tstarts=list(map(int, line[20].strip(',').split(',')))
		for start_number, start in enumerate(tstarts):
			coords.append([tname, start, start+sizes[start_number]])
coords.sort()

previous_chrom, max_end='junk', 0
region_count,nuc_count=0,0
for entry in coords:
	chrom, start, end=entry
	print(chrom, start, end, max_end)
	if chrom!=previous_chrom or start>max_end+100:
		if previous_chrom!='junk':
			region_stats.write(f'{previous_chrom}\t{region_start}\t{max_end}\t{max_end-region_start}\n')
		region_count+=1
		region_start=start
	overlapping=False
	if chrom==previous_chrom and start<max_end:
		overlapping=True
	if overlapping and end>max_end:
		nuc_count+=(end-max_end)
	elif not overlapping:
		nuc_count+=(end-start)
	converted_coords.write(f'{chrom}\t{start}\t{end}\t{region_count}\t{nuc_count}\n')
	if chrom==previous_chrom and end>max_end:
		max_end=end
	elif chrom!=previous_chrom:
		max_end=end
	previous_chrom=chrom
region_stats.write(f'{previous_chrom}\t{region_start}\t{max_end}\t{max_end-region_start}\n')

summary_out.write(f'total_regions:\t{region_count}\n')
summary_out.write(f'total_nucs:\t{nuc_count}\n')
