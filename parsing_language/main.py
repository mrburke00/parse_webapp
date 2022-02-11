from query_load import query_load, parseSequence
from itertools import groupby
from operator import itemgetter
import numpy as np
import re
import json
import operator

def loadQueryRules(query_rules_file):
	with open(query_rules_file) as json_file:
		data = json.load(json_file)
	return data

def slidingWindow(sequence, window_size):
    if len(sequence) <= window_size:
        return sequence
    for i in range(len(sequence) - window_size + 1):
        yield sequence[i:i+window_size]

def typeCheckQuery(query, query_rules):
	# WILL HAVE TO CHANGE FOR MULTIPLE QUERIES
	query_args = {}
	subject = query.split('??')[0].strip('[ ]').lower()
	target = query.split('??')[1].strip('[ ]').lower()
	options = query.split('??')[2].strip('[ ]').lower()
	query = subject + '??' + target

	try: subject in query_rules.keys()
	except: raise ValueError("Subject of query is not supported.")
	query_args['subject'] = subject

	letter_target = str(list(target)[0])
	try: letter_target in query_rules[subject].keys()
	except: raise ValueError("Target of query is not supported.")

	eq_delimiters = query_rules[subject][letter_target]['eq_delims']
	regexPattern = '|'.join(map(re.escape, eq_delimiters))
	try : eq_index = re.search(regexPattern, query)
	except : raise ValueError("Equality of query is not supported.")
	eq = query[eq_index.start():eq_index.end()]

	threshold = float(re.findall('\d*\.?\d+',target)[0])
	if not (query_rules[subject][letter_target].get('min') is None):
		try : threshold >= float(query_rules[subject][letter_target]['min']) and threshold <= float(query_rules[subject][letter_target]['max'])
		except : raise ValueError("Threshold of query is not supported.")
	query_args['target'] = [letter_target, eq, threshold]

	query_options = query_rules[subject]['options']
	required_options = [key for key,value in query_options.items() if value['required'] == 'True']

	for req_op in required_options:
		if req_op.lower() in options.lower():
			continue
		else:
			raise ValueError("Missing Required Option: ", req_op)
	opt_args = {}
	for opt in options.split(':'):
		try : eq_index = re.search(regexPattern, opt)
		except : raise ValueError("Equality of query is not supported.")
		eq = opt[eq_index.start():eq_index.end()]
		
		opt_args[opt[0:eq_index.start()]] =  [eq, opt[eq_index.end():]]

	query_args['options'] = opt_args
	return query_args
'''
###
Function to extract regions based on some percent content threshold
###
'''
def queryPercentile(query_args):
	subject = query_args['subject']
	sequence = subj_sequences[subject]
	target_letter, target_eq, target_thresh = query_args['target']
	target_thresh = float(target_thresh)
	options = query_args['options']

	N = int(options['n'][1])
	win_gen = slidingWindow(sequence, N)

	poss_regions = []
	for index, window in enumerate(win_gen):
		p_count = window.count(target_letter.upper())
		if perc_ops[target_eq](p_count/N,target_thresh) : 
			poss_regions.append(index) 

	if options['merge']:
		consec_regions = [] 
		# extract consecutive indexes 
		for k, g in groupby(enumerate(poss_regions), lambda ix : ix[0] - ix[1]):
			consec_regions.append(list(map(itemgetter(1), g)))
		
		ps_regions = []
		for r in range(0,len(consec_regions),2):
			# merge 
			region_len = len(sequence[consec_regions[r][-1]:consec_regions[r+1][0]])
			
			if perc_ops[target_eq](sequence[consec_regions[r][-1]:consec_regions[r+1][0]].count(target_letter.upper())/region_len,target_thresh):
				ps_regions.append(list(range(consec_regions[r][0], consec_regions[r+1][-1])))
			# don't merge
			else:
				ps_regions.append(consec_regions[r])
				ps_regions.append(consec_regions[r+1])
	else:
		ps_regions = poss_regions
	return ps_regions[0]

'''
###
Function executes pattern for "at least N" letters in row
###
'''
def atLeastPattern(query_args, window):
	subject = query_args['subject']
	target_letter, target_eq, target_thresh = query_args['target']
	target_thresh = int(target_thresh)
	options = query_args['options']
	target_pattern = ''.join([target_letter]*4)
	if target_pattern in window:
		return True
	return False
'''
###
Function to extract regions based on some sequence content (i.e. four A|D in a row)
###
'''	
def queryPattern(query_args):
	subject = query_args['subject']
	sequence = subj_sequences[subject]
	target_letter, target_eq, target_thresh = query_args['target']
	target_thresh = float(target_thresh)
	options = query_args['options']

	N = int(options['n'][1])
	win_gen = slidingWindow(sequence, N)

	print(subject, target_letter, target_eq, target_thresh, options)
	print(target_eq)
	poss_regions = []
	for index, window in enumerate(win_gen):
		if patt_ops[target_eq](query_args, window):
			poss_regions.append(index) 
	return poss_regions

def executeQuery(query_args):
	subject = query_args['subject']
	target_letter, target_eq, target_thresh = query_args['target']
	target_thresh = float(target_thresh)
	options = query_args['options']
	# check if equality is in our operator definition, if exists then percentile query
	if target_eq in perc_ops.keys():
		regions = queryPercentile(query_args)
	# if its not then means we have some sort of pattern query
	else:
		regions = queryPattern(query_args)
		# need to make new function to handle protein sequences
	#print(regions)
	return regions

def constructQuery(query, query_rules):
	queries_list = query.split('>>') ## need to not hardcode this
	query_args = []
	for query in queries_list:
		try: query_arg = typeCheckQuery(query, query_rules)
		except: raise ValueError("Query Failed")
		query_args.append(query_arg)
	print("+",query_args)
	return query_args
'''
# Function to pipe mutiple query results into one
'''
def projectRegions(regions):
	# need to have buffer option
	return list(set.intersection(*map(set,regions)))
 




#percentile query operations
perc_ops = { ">=": operator.ge, "=": operator.eq, "<=": operator.ge }
#pattern query opertaions
patt_ops = { "&=": atLeastPattern}
#subjects 
subj_sequences = { "pfd" : None, "protein" : None}
def main():

	subj_sequences['pfd'] = query_load()
	subj_sequences['protein'] = list(parseSequence("prot1.txt"))

	#query = '[PFD??P>=0.9??N=50:MERGE=True]'
	query = '[PFD??P>=0.9??N=50:MERGE=True]>>[Protein??A&=4??N=50]'
	#query = '[Protein??A&=4??N=50]'	
	
	query_rules = loadQueryRules('query_rules.json')
	queries = constructQuery(query, query_rules)
	regions = []
	for query in queries:
		print(query)
		regions.append(executeQuery(query))

	##FOR TESTING##
	regions[1] = list(range(300, 335))

	if len(regions) > 1:
		final_region = projectRegions(regions)

	print(final_region)



#[*][p>=50][*]
#[p>=50][*]
#[*][sss][p>50]
#[dsfasdf] | [adfasdaaaa]


if __name__ == "__main__":
    main()