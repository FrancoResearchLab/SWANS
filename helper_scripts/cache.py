import os
import sys
reference_type, reference = '', ''

config = sys.argv[1]

def get_config(search_config):

	# configs to exclude from lowercasing
	exclusion_list = (
		'RPATH',
		'CELLRANGER_REFERENCE',
		'TRANSFERDATA_REF_FILE',
		'REGRESSION_FILE',
		'USER_GENE_FILE',
	)

	with open('configs/prelim_configs.yaml') as input_file:
		for line in input_file.readlines():
			line = line.lstrip()
			if line.startswith(search_config):
				reference = line.split(search_config)[1].replace('\n', '')
				reference = reference.replace(':', '').lstrip().rstrip()
				# if the config is not in the exclusion_list, lowercase it
				if not line.startswith(exclusion_list):
					reference = reference.lower() 
	return(reference)

selection = get_config(config)
print(selection)
