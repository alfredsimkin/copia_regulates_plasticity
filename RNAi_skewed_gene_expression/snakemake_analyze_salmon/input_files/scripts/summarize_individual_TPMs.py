TPM_path=snakemake.params['TPM_path']
control_list=snakemake.params['control_list']
experimental_list=snakemake.params['experimental_list']
individual_summary=open(snakemake.output['individual_summary'], 'w')

def populate_dict(title_list, sample_list):
	sample_dict={}
	for sample in sample_list:
		title_list.append(f'{sample}_TPM')
		sample_dict.setdefault(sample, {})
		sample_path=TPM_path.replace('|placeholder|', sample)+'/quant.sf'
		for line_number, line in enumerate(open(sample_path)):
			line=line.strip().split('\t')
			if line_number>0:
				sample_dict[sample].setdefault(line[0], float(line[-2]))
	return sample_dict, title_list

control_dict, title_list=populate_dict(['transcript'], control_list)
experimental_dict, title_list=populate_dict(title_list, experimental_list)

for control in control_dict:
	for experimental in experimental_dict:
		title_list.append(f'{experimental}/{control}')
individual_summary.write('\t'.join(map(str, title_list))+'\n')

print(control_dict)
print(experimental_dict)
for transcript in control_dict[control_list[0]]:
	print(transcript)
	out_line=[transcript]
	for control in control_dict:
		out_line.append(control_dict[control][transcript])
	for experimental in experimental_dict:
		out_line.append(experimental_dict[experimental][transcript])
	for control in control_dict:
		control_value=control_dict[control][transcript]
		for experimental in experimental_dict:
			experimental_value=experimental_dict[experimental][transcript]
			if control_value==0:
				out_line.append('Und')
			else:
				out_line.append(experimental_value/control_value)
	individual_summary.write('\t'.join(map(str, out_line))+'\n')
