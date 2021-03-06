
import os
import sys
import random
parse = os.path.abspath('parse_src')
sys.path.insert(1,parse)
import fortran_parse
import json
from io import StringIO
'''
Remove old artifact files with each new query 
'''
def clean():
	os.system('rm -r residue_level_rmodel*')
	os.system('rm -r /cmd2web/src/web_src/static/js/residue_level_rmodel_new*')

'''
For non file sequence queries, but can still handle multiple protein sequences
Parse each input sequence and send to Parse.f script
Get each return and send to fortran_parse.process_parse() display 
'''
''' ##broken##
def fortran_wrap(temp):
	clean()
	cmds = cmd.split()
	f_call = cmds[0:6]
	cmds = cmds[6:]
	json_list = []
	i = 1
	for cmd in cmds:
		if '>' in cmd:
			name = cmd.strip()
		else:
			out_file_name = '/tmp/' + str(random.randint(0,sys.maxsize)) + '.out'
			f_cmd = " ".join(f_call)
			f_cmd += " " + cmd
			f = open(out_file_name, 'w')
			os.system(f_cmd + ' > ' + out_file_name)
			f.close()
			temp = []
			out = []
			for line in open(out_file_name, 'r'):
				print(line)
				if line.strip() != '':
					temp.append(line.strip())
			os.system('mv residue_level_rmodel.csv residue_level_rmodel_'+str(i)+'.csv')
			out = fortran_parse.process_parse(temp,i, name)
			json_list.append(out)
			i = i + 1
	return json.dumps(json_list)
'''
def fortran_wrap_single(temp):
	os.system('mv residue_level_rmodel.csv residue_level_rmodel_1.csv')
	out = fortran_parse.process_parse(temp,1,'sample1')
	clean()
	return json.dumps(out)

def fortran_wrap_file_single(cmd, filename = None):
	clean()
	cmds = cmd.split()
	f_call = cmds[0:6]
	cmds = cmds[6:]
	json_list = []
	k = 1
	tmp = str(cmds[0])
	queries = tmp.split('\\n')

	for i in range(0,len(queries),2):
		name = queries[i].strip()
		sequence = queries[i+1].strip()
		out_file_name = '/tmp/' + str(random.randint(0,sys.maxsize)) + '.out'
		f_cmd = " ".join(f_call)
		f_cmd += " " + sequence.strip()
		f = open(out_file_name, 'w')
		print(f_cmd)
		os.system(f_cmd + ' > ' + out_file_name)
		f.close()
		temp = []
		out = []
		for line in open(out_file_name, 'r'):
			if line.strip() != '':
				temp.append(line.strip())
		os.system('mv residue_level_rmodel.csv residue_level_rmodel_'+str(k)+'.csv')
		out = fortran_parse.process_parse(temp,k,name)
		json_list.append(out)
		k = k + 1
	return json.dumps(json_list)


'''
For file specific protein sequence input
Parse each sequence in file and send to Parse.f script
Get each return and send to fortran_parse.process_parse() for web display 
'''
''' ##broken##
def fortran_wrap_file(cmd, filename = None):
	clean()
	cmds = cmd.split()
	f_call = cmds[0:6]
	cmds = cmds[6:]
	json_list = []
	k = 1
	if filename is None:
		filename = cmds[0]
	with open(filename, 'r') as infile:
		lines = [line for line in infile]
		for i in range(0,len(lines),2):
			name = lines[i].strip()
			sequence = lines[i+1].strip()
			out_file_name = '/tmp/' + str(random.randint(0,sys.maxsize)) + '.out'
			f_cmd = " ".join(f_call)
			f_cmd += " " + sequence.strip()
			f = open(out_file_name, 'w')
			os.system(f_cmd + ' > ' + out_file_name)
			f.close()
			temp = []
			out = []
			for line in open(out_file_name, 'r'):
				if line.strip() != '':
					temp.append(line.strip())
			os.system('mv residue_level_rmodel.csv residue_level_rmodel_'+str(k)+'.csv')
			out = fortran_parse.process_parse(temp,k,name)
			json_list.append(out)
			k = k + 1
	return json.dumps(json_list)
'''

#cmd = 'gfortran -o /Users/DBurke/Documents/Layerlab/parse_webapp/cmd2web/src/web_src/static/js/Parse.exe /Users/DBurke/Documents/Layerlab/parse_webapp/cmd2web/src/web_src/static/js/Parse.f && /Users/DBurke/Documents/Layerlab/parse_webapp/cmd2web/src/web_src/static/js/./Parse.exe test_seq.txt'
#fortran_wrap_file(cmd)
