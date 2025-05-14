'''
Author:	E. Reichenberger
Date:		10.17.2019		

Purpose: Re-usable script to determine which samples will be analyzed. If samples.sample_list is present, get the files from there, if not, do a greedy glob search. 

Required: 
'''

import glob
import os.path
import os

#files
samplefile_list = 'samples.sample_list'

def get_samples(project):
	raw_matrix_path = ''

	#paths
	raw_matrix_path = 'data/endpoints/' + project + '/'

	SAMPLE_LIST, number_list = [], [] #list variables

	#if file exists
	#if os.path.isdir(samplefile_list) == True: #this is wrong.
	if os.path.isfile(samplefile_list) == True:
		# print(samplefile_list)
		with open(samplefile_list, 'r') as input_file:
			input_file.readline()
			for line in input_file.readlines():
					line = line.lstrip().rstrip()
					if len(line) > 1: #grrr....check for empty lines
						#SAMPLE_LIST.append(line)
						sLine = line.split('\t')
						sample = sLine[0]
						condition = sLine[1]
						path = sLine[2]
						SAMPLE_LIST.append(sample)

	#if file dn exists
	elif os.path.isfile(samplefile_list) == False:
		SAMPLE_LIST = os.listdir(raw_matrix_path)
		#SAMPLE_LIST = os.path.basename(raw_matrix_path)

	SAMPLE_LIST = list(set(SAMPLE_LIST))
	
	if 'analysis' in SAMPLE_LIST:
		SAMPLE_LIST.remove('analysis')

	# HUMAN SORT
	for S in SAMPLE_LIST:
		S = S.replace(raw_matrix_path, '')
		rem = ''
		for s in S:
			if s.isdigit():
				rem = rem + s
				#rem = int(rem)
		if rem not in number_list:
			number_list.append(rem)
 
	#SAMPLE_LIST = [x for _,x in sorted(zip(number_list, SAMPLE_LIST))] This does not work for samples names A1, A2, ....

	if os.path.isfile(samplefile_list) == False:
		with open(samplefile_list, 'w') as output_file:
			for S in SAMPLE_LIST:
				output_file.write(S + '\n')

	return(SAMPLE_LIST)
