
import os
import sys
import random
parse = os.path.abspath('parse_src')
sys.path.insert(1,parse)
import fortran_parse
import json

def clean():
	os.system('rm -r residue_level_rmodel*')
	os.system('rm -r /cmd2web/src/web_src/static/js/residue_level_rmodel_new*')

def fortran_wrap(cmd):
	clean()
	cmds = cmd.split()
	f_call = cmds[0:6]
	cmds = cmds[6:]
	json_list = []
	i = 1
	for cmd in cmds:
		if cmd != '>':
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
			out = fortran_parse.process_parse(temp,i)
			json_list.append(out)
			i = i + 1
		#return json.dumps(out)
		else:
			continue
	return json.dumps(json_list)

def fortran_wrap_file(cmd):
	clean()
	cmds = cmd.split()
	f_call = cmds[0:6]
	cmds = cmds[6:]
	json_list = []
	k = 1
	with open(cmds[0], 'r') as infile:
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


#cmd = 'gfortran -o /Users/DBurke/Documents/Layerlab/parse_webapp/cmd2web/src/web_src/static/js/Parse.exe /Users/DBurke/Documents/Layerlab/parse_webapp/cmd2web/src/web_src/static/js/Parse.f && /Users/DBurke/Documents/Layerlab/parse_webapp/cmd2web/src/web_src/static/js/./Parse.exe test_seq.txt'
#fortran_wrap_file(cmd)
