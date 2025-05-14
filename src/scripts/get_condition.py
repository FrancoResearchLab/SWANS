def get_condition(query_sample, filename):
	dic = {}

	with open(filename, 'r') as input_file:
		input_file.readline()

		for line in input_file.readlines():
			line = line.replace('\n', '')
			sLine = line.split('\t')
			sample = sLine[0]
			condition = sLine[1]

			if sample not in dic:
				dic[sample] = ''
			dic[sample] = condition

		cond = dic[query_sample]

		return(cond)
