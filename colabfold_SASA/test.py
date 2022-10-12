import os

from Colabfold_SASA import Colabfold_SASA
# from lambo.candidate import FoldedCandidate
from Bio.PDB import PDBParser

def get_complex_pdb_name(path):
    name_list = os.listdir(path)
    pdb_name = ''
    for name in name_list:
        # if name.startswith('test_9216f_unrelaxed_rank_1') and name[-4:] == '.pdb':
        if 'unrelaxed_rank_1' in name and name[-4:] == '.pdb':
            pdb_name = name
    if pdb_name == '':
        return ValueError
    else:
        return pdb_name

def copy_pdb():
    with open(r'/home/qin/source_code/lambo_osc/lambo-main/lambo/assets/fpbase/proxy_rfp_seed_data_peptide.csv', 'r') as f:
        lines = f.readlines()
        lines = lines[1:]
    all_seq=[]
    for line in lines:
        all_seq.append(line.split(',')[0])
    for seq in all_seq:
        print(seq)
        result_path = os.path.join(r'/home/qin/source_code/lambo-source/lambo-main/colabfold_SASA/result', seq)
        complex_pdb_name = get_complex_pdb_name(result_path)
        complex_path = os.path.join(result_path, complex_pdb_name)
        peptide_path=os.path.join(result_path, 'peptide.pdb')
        target_path='result1/'+seq
        os.system('mkdir '+target_path)
        os.system('cp '+complex_path+' '+target_path)
        os.system('cp ' + peptide_path + ' ' + target_path)

# copy_pdb()

def binding_ratio(pdbfile, mcherrystart_site=29, peptide_chain='B', AKT_chain='C', AKT_domain_start=149,
                  AKT_domain_end=409, CA_only=False):
    '''
    pdbfile: pdbfile position for AKT and peptide complex.
    in result pdb file, peptide and mchaeery are concatenate, and mcherrystart_site is the start site for mcherry.
    peptide_chain indicate the peptide chain_AKT
    AKT_chain indicate the AKT chain
    if mcherrystart_site=-1, no mcherry
    AKT_domain_start: domain start position on AKT
    AKT_domain_end: domain end position on AKT
    CA_only: whether consider all atoms or CA only, default is False same with current ppt
    '''
    parser = PDBParser()
    # read structure from file
    structure = parser.get_structure('peptide_AKT', pdbfile)

    model = structure[0]
    chain_peptide = model[peptide_chain]
    chain_AKT = model[AKT_chain]

    bind_num = 0
    cutoff = 5
    peptide_bind_sites = {}
    AKT_bind_sites = {}
    AKT_bind_sites_domain = {}
    for residue1 in chain_peptide:
        peptide_pos = residue1.id[1]
        if mcherrystart_site != -1:
            if peptide_pos >= mcherrystart_site:
                break

        for residue2 in chain_AKT:
            AKT_pos = residue2.id[1]
            if CA_only:
                # compute distance between CA atoms
                try:
                    distance = residue1['CA'] - residue2['CA']
                except KeyError:
                    ## no CA atom, e.g. for H_NAG
                    print("no ca ")
                    continue

                if distance <= cutoff:
                    peptide_bind_sites[peptide_pos] = distance
                    AKT_bind_sites[AKT_pos] = distance
                    if AKT_pos >= AKT_domain_start and AKT_pos <= AKT_domain_end:
                        AKT_bind_sites_domain[AKT_pos] = distance

            else:
                # compute distance between all atoms
                for atom1 in residue1.get_atoms():
                    for atom2 in residue2.get_atoms():
                        distance = atom1 - atom2
                        if distance <= cutoff:
                            peptide_bind_sites[peptide_pos] = distance
                            AKT_bind_sites[AKT_pos] = distance
                            if AKT_pos >= AKT_domain_start and AKT_pos <= AKT_domain_end:
                                AKT_bind_sites_domain[AKT_pos] = distance

                            break

    num_of_AKT_bind_sites = len(AKT_bind_sites)
    num_of_peptide_bind_sites = len(peptide_bind_sites)
    num_of_AKT_bind_sites_domain = len(AKT_bind_sites_domain)
    # print(num_of_AKT_bind_sites_domain)
    # print(num_of_AKT_bind_sites)
    return num_of_AKT_bind_sites_domain / (AKT_domain_end - AKT_domain_start + 1), num_of_AKT_bind_sites / 480


def copy_fa():
    with open('list.lst','r')as f:
        lines=f.readlines()
        for line in lines:
            if line[0]!='[':
                fa_path = os.path.join('fasta', line[:-1]+'.fa')
                result_path = os.path.join('fa_pdb', line[:-1]+'.fa')
                os.system('cp '+fa_path+' '+result_path)

# copy_fa()


def calculate_ratio_bingding_sites():
    all_seq=[]
    with open(r'/home/qin/source_code/lambo-main/lambo/assets/fpbase/proxy_rfp_seed_data_peptide.csv', 'r') as f:
        lines = f.readlines()
        lines = lines[1:]
    for line in lines:
        all_seq.append(line.split(',')[0])
    for seq in all_seq:
        print(seq)
        result_path=os.path.join(r'/home/qin/source_code/lambo-source/lambo-main/colabfold_SASA/result',seq)
        complex_pdb_name = get_complex_pdb_name(result_path)
        complex_path=os.path.join(result_path,complex_pdb_name)
        print(binding_ratio(complex_path))

# calculate_ratio_bingding_sites()

def copy_pdb():
    with open('list.lst','r')as f:
        lines=f.readlines()
        for line in lines:
            if line[0]!='[':
                result_path = os.path.join('result', line[:-1])
                result_list=os.listdir(result_path)
                pdb_name=''
                for name in result_list:
                    if name.startswith(line[:-1] + '_unrelaxed_rank_1') and name[-4:] == '.pdb':
                        pdb_name = name
                print()
                result_path = os.path.join('result', line[:-1],pdb_name)
                os.system('cp '+result_path+' '+os.path.join('pdb',pdb_name))

def test_colabfold():
    fa_path=os.path.join('./fasta','QNGGIHRLCQLLSKAHEMVQRRLALAPDA.fa')
    result_path=os.path.join('./result','QNGGIHRLCQLLSKAHEMVQRRLALAPDA')
    cmd='colabfold_batch --use-gpu-relax --num-recycle 3 --model-type AlphaFold2-multimer-v2 '+fa_path+' '+result_path
    os.system(cmd)

def calculate_SASA():
    all_seq=[]
    with open(r'/home/qin/source_code/lambo-main/lambo/assets/fpbase/proxy_rfp_seed_data_peptide.csv', 'r') as f:
        lines = f.readlines()
        lines = lines[1:]
    for line in lines:
        all_seq.append(line.split(',')[0])
    for seq in all_seq:
        print(seq)
        print(Colabfold_SASA(seq).run_colab_sasa())
# calculate_SASA()



def run():
    fa_list=os.listdir('fa')
    for fa in fa_list:
        fa_path=os.path.join('fa',fa)
        result_path='result/'+fa[:-3]
        cmd = 'bash -i run_complex_osc.sh ' + fa_path + ' ' + result_path
        os.system(cmd)

def getfasta():
    all_seq=[]
    with open(r'/home/qin/source_code/lambo_osc/lambo-main/lambo/assets/fpbase/proxy_rfp_seed_data_peptide.csv','r')as f:
        lines=f.readlines()
        lines=lines[1:]
    for line in lines:
        all_seq.append(line.split(',')[0])
    # with open(r'/home/qin/source_code/lambo_osc/lambo-main/lambo/assets/fpbase/rfp_known_structures_peptide.csv', 'r') as f:
    #     lines = f.readlines()
    #     lines = lines[1:]
    # for line in lines:
    #     all_seq.append(line.split(',')[0])
    akt_seq = 'MSDVAIVKEGWLHKRGEYIKTWRPRYFLLKNDGTFIGYKERPQDVDQREAPLNNFSVAQCQLMKTERPRPNTFIIRCLQWTTVIERTFHVETPEEREEWTTAIQTVADGLKKQEEEEMDFRSGSPSDNSGAEEMEVSLAKPKHRVTMNEFEYLKLLGKGTFGKVILVKEKATGRYYAMKILKKEVIVAKDEVAHTLTENRVLQNSRHPFLTALKYSFQTHDRLCFVMEYANGGELFFHLSRERVFSEDRARFYGAEIVSALDYLHSEKNVVYRDLKLENLMLDKDGHIKITDFGLCKEGIKDGATMKTFCGTPEYLAPEVLEDNDYGRAVDWWGLGVVMYEMMCGRLPFYNQDHEKLFELILMEEIRFPRTLGPEAKSLLSGLLKKDPKQRLGGGSEDAKEIMQHRFFAGIVWQHVYEKKLSPPFKPQVTSETDTRYFDEEFTAQMITITPPDQDDSMECVDSERRPHFPQFSYSASGTA'

    # ready_path=os.listdir('./')
    seqs=set(all_seq)
    for seq in seqs:
        with open(os.path.join('fa',seq+'.fa'),'w')as f:
            f.write('>'+seq+'\n'+seq+':'+akt_seq)
        # os.system('mkdir result_fa/'+seq)
        print('colabfold_batch --use-gpu-relax --num-recycle 3 --model-type AlphaFold2-multimer-v2 fa/'+seq+'.fa result_fa/'+seq)
    print()
# getfasta()


def calculate():
    all_seq=[]
    # with open(r'/home/qin/source_code/lambo-main/lambo/assets/fpbase/proxy_rfp_seed_data_peptide.csv', 'r') as f:
    #     lines = f.readlines()
    #     lines = lines[1:]
    # for line in lines:
    #     all_seq.append(line.split(',')[0])
    for seq in ['ASDADSDSDS']:
        print(seq)
        print(Colabfold_SASA(seq).run_colab_sasa())
# calculate()


def get_metrics_by_fa_path():
    file_names=os.listdir('fasta')
    names=[]
    for name in file_names:
        names.append(name[:-3])
    for name in names:
        print(name)
        print(Colabfold_SASA(name).run_colab_sasa())


get_metrics_by_fa_path()