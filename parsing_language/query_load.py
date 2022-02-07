import os

def parseSequence(filename):
	with open(filename) as f:
		sequence = "".join(line.strip() for line in f)
	return sequence

def executeParse(sequence, parse_version):
	cmd = 'gfortran -o ' + parse_version +'.exe ' +  parse_version + '.f && ./' + parse_version + '.exe ' + sequence + ' > prot1_result.txt'
	os.system(cmd)

def parseResult(filename):
	with open(filename) as f:
		sequence = f.readlines()
	sequence_result = list(sequence[1])
	return sequence_result

def query_load(run_points = False, sequence_file = "prot1.txt", parse_version = "Parse_v2", result_file = "prot1_result.txt"):

	if run_points:
		sequence = parseSequence(sequence_file)
		executeParse(sequence, parse_version)
	sequence = parseResult(result_file)
	return sequence



