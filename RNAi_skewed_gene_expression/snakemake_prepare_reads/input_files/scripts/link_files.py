import yaml
import subprocess
import custom2
yaml_dict=yaml.load(open(snakemake.input['input_yaml']), Loader=yaml.SafeLoader)
print('starting linking script')
for output_folder in yaml_dict['replicates']:
	infiles=[]
	input_pattern=yaml_dict['replicates'][output_folder]
	subprocess.call(['mkdir', '-p', 'output_files/fastq_links/'+output_folder])
	split_pattern=input_pattern.split('/')
	folder='/'.join(split_pattern[:-1])
	pattern=split_pattern[-1]
	roots=custom2.file_walker(folder)
	for root, dirs, files in roots:
		for file in files:
			if pattern in file:
				infiles.append(root+'/'+file)
	rep=output_folder.split('_')[-1]
#	print('output_folder is', output_folder)
	for path in infiles:
		desired_part='_'.join(path.split('_')[-3:])
#		print('file is: output_files/fastq_links/'+output_folder+'/'+desired_part)
		subprocess.call(['ln', '-s', '-T', path, 'output_files/fastq_links/'+output_folder+'/'+desired_part])

