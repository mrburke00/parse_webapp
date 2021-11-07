import pandas as pd
import os

def manipulate_csv(count):
    parse = os.path.abspath('residue_level_rmodel_'+str(count)+'.csv')
    df = pd.read_csv(parse, header = None)
    df.rename(columns={0: 'col1', 1:'col2', 2: 'col3', 3: 'col4'}, inplace=True)
    parent = os.path.abspath(os.getcwd())
    df.to_csv(parent + '/cmd2web/src/web_src/static/js/residue_level_rmodel_new_'+str(count)+'.csv', index=False) # save to new csv file



def process_parse(temp, count, name):

    out = {}
    out['result'] = temp[0]
    line = temp[1].split()
    out['sample_name'] = name
    out['sequence_length'] = line[2]

    res = [temp.index(i) for i in temp if 'Whole' in i]
    names = ['nu_model','b_turn','r_model']
    for index,i in enumerate(res):
        line = temp[i].split()
        out[names[index]] = line[3]
        
    names = ['blue','red','black']
    res = [temp.index(i) for i in temp if 'Number' in i]
    for index,i in enumerate(res):
        line = temp[i].split()
        out[names[index]] = line[6]

    res = [temp.index(i) for i in temp if 'domain' in i]
    last = res[6:][0]
    rows = []
    for i in range(last,len(temp)):
        new_row = []
        x = temp[i].split()
        indexes = [1,2, 6,9,11]
        for index in indexes:
            if index == 2:
                new = x[2].strip(',') + ' ' + x[3].strip(',')
                new_row.append(new)
            else:    
                new_row.append(x[index].strip(','))
        rows.append(new_row)
    out['table'] = rows
    manipulate_csv(count)
    return out
