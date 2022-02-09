from query_load import query_load
from itertools import groupby
from operator import itemgetter
import numpy as np
import re
import json
import operator

ops = { ">=": operator.ge, "=": operator.eq, "<=": operator.ge }


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

def queryPercentile(sequence, query_args):
	subject = query_args['subject']
	target_letter, target_eq, target_thresh = query_args['target']
	target_thresh = float(target_thresh)
	options = query_args['options']

	N = int(options['n'][1])
	win_gen = slidingWindow(sequence, N)

	poss_regions = []
	for index, window in enumerate(win_gen):
		p_count = window.count(target_letter.upper())
		if ops[target_eq](p_count/N,target_thresh) : 
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
			
			if ops[target_eq](sequence[consec_regions[r][-1]:consec_regions[r+1][0]].count(target_letter.upper())/region_len,target_thresh):
				ps_regions.append(list(range(consec_regions[r][0], consec_regions[r+1][-1])))
			# don't merge
			else:
				ps_regions.append(consec_regions[r])
				ps_regions.append(consec_regions[r+1])
	else:
		ps_regions = poss_regions
	return ps_regions

def executeQuery(query_args, sequence):
	subject = query_args['subject']
	target_letter, target_eq, target_thresh = query_args['target']
	target_thresh = float(target_thresh)
	options = query_args['options']
	# perform some sort of threshold metric
	if target_eq in ops.keys():
		regions = queryPercentile(sequence, query_args)
	# perform some sort of sequence or pattern metric
	else:
		print('issue')
		# need to make new function to handle protein sequences

	return regions



def constructQuery(query, query_rules, sequence):
	try: query_args = typeCheckQuery(query, query_rules)
	except: raise ValueError("Query Failed")
	print(query_args)

	#ptions = options.split(':')
	#print(options)
	#for req_op in required_options:
	#	break

		#	targets = ['P', 'F', 'D']
		#	target = next(t for t in targets if t in params)
		#	threshold = params[eq_index.end():]
		#	queries['eq'] = eq
		#	queries['target'] = target
		#	queries['threshold'] = threshold


	'''	
	queries = {}
		options = query.split('?')[1:]
		params = query.split('?')[0]

		if options:
			for option in options:
				queries[option.split('=')[0]] = option.split('=')[1]

		if params:
			eq_delimiters = [">=", "=>", "=", "<=", "=<"]
			regexPattern = '|'.join(map(re.escape, eq_delimiters))
			eq_index = re.search(regexPattern, params)
			eq = params[eq_index.start():eq_index.end()]
			targets = ['P', 'F', 'D']
			target = next(t for t in targets if t in params)
			threshold = params[eq_index.end():]
			queries['eq'] = eq
			queries['target'] = target
			queries['threshold'] = threshold
		return queries			
	'''

'''def extractPS(sequence, queries):
	threshold = float(queries['threshold'])
	N = int(queries['N'])

	win_gen = slidingWindow(sequence, N)

	poss_regions = []
	for index, window in enumerate(win_gen):
		p_count = window.count(queries['target'])
		if p_count/N >= threshold: 
			poss_regions.append(index) 
			print(index, index + 50)

	print(poss_regions)
	consec_regions = [] 
	# extract consecutive indexes
	for k, g in groupby(enumerate(poss_regions), lambda ix : ix[0] - ix[1]):
		consec_regions.append(list(map(itemgetter(1), g)))
	

	ps_regions = []
	for r in range(0,len(consec_regions),2):
		# merge 
		if sequence[consec_regions[r][-1]:consec_regions[r+1][0]].count(queries['target']) >= threshold:
			ps_regions.append(list(range(consec_regions[r][0], consec_regions[r+1][-1])))
		# don't merge
		else:
			ps_regions.append(consec_regions[r])
			ps_regions.append(consec_regions[r+1])
	return ps_regions'''
#[*][p>=50][*]
#[p>=50][*]
#[*][sss][p>50]
#[dsfasdf] | [adfasdaaaa]

def main():

	sequence = query_load()
	query = '[PFD??P>=0.9??N=50:MERGE=True]'
	query_rules = loadQueryRules('query_rules.json')
	queries = constructQuery(query, query_rules, sequence)
	#ps_regions = extractPS(sequence, queries)
	#print(ps_regions)

if __name__ == "__main__":
    main()