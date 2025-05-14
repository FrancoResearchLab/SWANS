'''
Author:	E. Reichenberger
Date:		10.17.2019		

Purpose: Re-usable script to determine which samples will be analyzed. If samples.sample_list is present, get the files from there, if not, do a greedy glob search. 

Required: 
'''

import glob
import os.path
import os
import sys
from pathlib import Path   

project = sys.argv[1]
samplefile_list = 'samples.sample_list'
starting_data = sys.argv[2] #starting data (cellranger or matrix) passed as argument
run_cellranger = sys.argv[3]

def get_samples(project, starting_data):

	SAMPLE_LIST = [] #list variables
	sample_path_dic = {}

	if os.path.isfile(samplefile_list) == True:
		with open(samplefile_list, 'r') as input_file:
			input_file.readline()
			for line in input_file.readlines():
					line = line.lstrip().rstrip()
					if len(line) > 1: #grrr....check for empty lines
						sLine = line.split('\t')
						sample = sLine[0]
						condition = sLine[1]
						path = sLine[2] #+ '/outs/' #add this dir as o

						SAMPLE_LIST.append(sample)
						if sample not in sample_path_dic:
							sample_path_dic[sample] = path
							
		with open('helper_scripts/make_symbolic_links.sh', 'w') as output_file:
			for s in SAMPLE_LIST:
				path_list = ''

				if starting_data.lower() == 'fastq' and run_cellranger == 'y':
					path_list = 'data/endpoints/' + project + '/' +  s  + '/fastq/'
				if starting_data.lower() == 'cellranger':
					path_list = 'data/endpoints/' + project + '/' +  s  + '/cellranger/'
				if starting_data.lower() == 'matrix':
					path_list = 'data/endpoints/' + project + '/' +  s  + '/matrix/'

				out_line = 'cd ' + path_list + '\n'
				output_file.write(out_line)

				# add final / if it's not there (so it actually copies files, not the dir!)
				# verbose, will tighten later
				if sample_path_dic[s].endswith('/'):
					out_line = 'ln -s ' + sample_path_dic[s] + '* .\n'

				if not sample_path_dic[s].endswith('/'):
					out_line = 'ln -s ' + sample_path_dic[s] + '/* .\n'

				output_file.write(out_line)
				out_line = 'cd - \n\n'
				output_file.write(out_line)

get_samples(project, starting_data)
