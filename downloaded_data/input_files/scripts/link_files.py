import yaml
import subprocess
import os
root=os.getcwd()
yaml_dict=yaml.load(open(snakemake.input['yaml_file']), Loader=yaml.SafeLoader)
status_path=snakemake.output['done_linking']
final_file=open(status_path, 'w')

for output_path in yaml_dict:
	input_path=root+'/'+yaml_dict[output_path]
	subprocess.call(['ln', '-s', '-f', '-T', input_path, output_path])
	final_file.write(f'finished linking {input_path} to {output_path}\n')
