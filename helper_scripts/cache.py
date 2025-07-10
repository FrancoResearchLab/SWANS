import os
import sys
reference_type, reference = '', ''

config = sys.argv[1]

def get_config(search_config):
	with open('configs/prelim_configs.yaml') as input_file:
		for line in input_file.readlines():
			line = line.lstrip()
			if line.startswith(search_config):
				reference = line.split(search_config)[1].replace('\n', '')
				reference = reference.replace(':', '').lstrip().rstrip()
				# if the config is not the RPATH, lower it
				if not line.startswith('RPATH'):
					reference = reference.lower() # guarantees everything is lowercase
				#if search_config.startswith('PROJECT'):
				#	reference = reference.lower()
	return(reference)

selection = get_config(config)
print(selection)
