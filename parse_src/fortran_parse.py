

def process_parse(temp):

    out = {}

    out['result'] = temp[0]
    line = temp[1].split('    ')
    out['sequence_length'] = line[1]

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

    return out