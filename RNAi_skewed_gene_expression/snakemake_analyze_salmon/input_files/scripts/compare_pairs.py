import scipy.stats
control_file=snakemake.input['control']
experimental_file=snakemake.input['experimental']
output_name=snakemake.output['stat_file']
fold_change_num=int(snakemake.params['fold_change_smaller'])
fold_change_denom=int(snakemake.params['fold_change_bigger'])
dec_thresh=fold_change_num/fold_change_denom
inc_thresh=fold_change_denom/fold_change_num

stat_file=open(output_name, 'w')

def make_exp_dict(input_file):
	exp_dict={}
	for line_number, line in enumerate(open(input_file)):
		if line_number>0:
			line=line.strip().split()
			exp_dict[line[0]]=float(line[3])+1 #adding a pseudocount of 1
	return exp_dict

control_dict=make_exp_dict(control_file)
experimental_dict=make_exp_dict(experimental_file)

greater_count, lesser_count, equal_count, total_count=0,0,0,0
for gene in control_dict:
	if gene in experimental_dict:
		if experimental_dict[gene]/control_dict[gene]>inc_thresh:
			greater_count+=1
		elif experimental_dict[gene]/control_dict[gene]<dec_thresh:
			lesser_count+=1
		else:
			equal_count+=1
		total_count+=1

stat_file.write(f'threshold for decrease in expression is: {dec_thresh}\n')
stat_file.write(f'threshold for increase is: {inc_thresh}\n')
stat_file.write(f'total genes found in both is: {total_count}\n')
stat_file.write('genes that were unchanged in experimental condition is: '+str(equal_count)+'\n')
stat_file.write('genes that went up in experimental condition is: '+str(greater_count)+'\n')
stat_file.write('genes that went down in experimental condition is: '+str(lesser_count)+'\n')

fewer_prob=scipy.stats.binom.cdf(lesser_count-1, lesser_count+greater_count, 0.5)
equal_or_greater=1-fewer_prob
stat_file.write('assuming an equal probability of an increase in expression as a decrease of expression\n')
stat_file.write('and that the only options are an increase in expression or a decrease of expression\n')
stat_file.write('(ignoring genes that do not change), then the probability of seeing '+str(lesser_count)+'\n')
stat_file.write('or more genes go down, in a total of '+str(lesser_count+greater_count)+' genes is: '+str(equal_or_greater)+'\n')
stat_file.close()
