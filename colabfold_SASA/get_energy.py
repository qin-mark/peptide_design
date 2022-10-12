import shutil
import os
import csv
substring = '.pdb'
substring2 = '.fxout'
def data_needed(filePath):
    file_name = list()
    for i in os.listdir(filePath):        #filename
        data_collect = ''.join(i)
        file_name.append(data_collect)
    return(file_name)

def copyfile(srcfile,dstpath):                       # copyfile
    if not os.path.isfile(srcfile):
        print ("%s not exist!"%(srcfile))
    else:
        fpath,fname=os.path.split(srcfile)
        if not os.path.exists(dstpath):
            os.makedirs(dstpath)
        shutil.copy(srcfile, dstpath + fname)

def get_complex_out_name(path,all_seq):
    name_list = os.listdir(path)
    # print(name_list)
    pdb_name = ''
    for name in name_list:
        # if name.startswith('test_9216f_unrelaxed_rank_1') and name[-4:] == '.pdb':
        if all_seq in name and name[-11:] == '_0_ST.fxout':
            output_name = name
    if output_name == '':
        return ValueError
    else:
        # print('name:'+pdb_name)
        return output_name


#submit to foldx
def caculate_energy(all_seq,result_path,foldx_path):
    pdb_path = '' + result_path + '/' + all_seq + ''
    pdb_list = data_needed(pdb_path)
    #print(pdb_list)
    # get energy score
    # print(all_seq)
    # submit to foldx
    for j in range(0, len(pdb_list)):
        if pdb_list[j].find(substring) != -1:  # find pdb
            copyfile(
                '' + result_path + '/' + all_seq + '/' + pdb_list[j] + '',
                '' + foldx_path + '/')
            command='foldx --command=Stability --pdb=' + pdb_list[j] + ''
            #command = './foldx --command=Stability --pdb=' + pdb_list[j] + ''   #linux
            #command = 'foldx -h'
            os.system(command)
    foldx_list =get_complex_out_name(''+foldx_path+'',all_seq)
    with open('' + foldx_path + '/' + foldx_list + '', 'r', encoding='utf-8') as fp2:#find out_put file
        for line in fp2:
            a = line.split()
            energy = a[1:2]  # get energy
            return energy



foldx_path = r'E:\1Code\lambo-main\foldx'
csv_path = r'E:\1Code\lambo-main\lambo\assets\fpbase\rfp_known_structures_peptide.csv'
result_path = r'E:\1Code\lambo-main\colabfold_SASA\result'
with open(csv_path,'r') as fp:
    lines = fp.readlines()
    lines = lines[1:]
    row = 0
    #print(lines)
    for line in lines:
        row = row + 1
        #count=0
        if line.split(',')[21]!='':
            continue
        #try:
            #print(row,line.split(',')[21])          #total_energy=null
        else:
            all_seq=line.split(',')[14]
            pdb_path = os.path.join(result_path, line.split(',')[14])
            #print(pdb_path)
            energy=caculate_energy(all_seq,result_path,foldx_path)
            #print(energy)
            #energy='1111111111'
            #print(all_seq,energy)
            #r = csv.reader(open(r''+csv_path+''))  # Here your csv file
            #lines = list(r)
            #l=line
            #b = ','
            #l=l.replace('\n',b)+str(energy)
            #lines[row]=[l]
            #print(lines[row])
            #writer = csv.writer(open(r'' + csv_path + '', 'w', newline=''))
            #writer.writerows(lines)

            r = csv.reader(open(
                csv_path))  # Here your csv file
            lines = list(r)
            # print(lines[2])
            lines[row][21] = energy
            print(lines[row][21])
            writer = csv.writer(open(
                csv_path,
                'w', newline=''))
            writer.writerows(lines)