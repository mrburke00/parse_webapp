from query_load import query_load
from itertools import groupby
from operator import itemgetter
import numpy as np
import re

def slidingWindow(sequence, window_size):
    if len(sequence) <= window_size:
        return sequence
    for i in range(len(sequence) - window_size + 1):
        yield sequence[i:i+window_size]

def extractPS(sequence, queries):
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
	return ps_regions

def constructQuery(query):
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
				

#[*][p>=50][*]
#[p>=50][*]
#[*][sss][p>50]
#[dsfasdf] | [adfasdaaaa]
def main():

	sequence = query_load()
	query = 'P>=0.9?N=50'
	queries = constructQuery(query)
	ps_regions = extractPS(sequence, queries)
	print(ps_regions)

if __name__ == "__main__":
    main()