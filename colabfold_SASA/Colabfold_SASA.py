import os
import uuid
from typing import Optional, Dict
from Bio.PDB import PDBParser, SASA
from Bio import PDB
import subprocess
from Bio.PDB import PDBParser

import our_settings
from lambo.candidate import FoldedCandidate



import torch
#colabfold_batch --templates --use-gpu-relax --num-recycle 3 --model-type AlphaFold2-multimer-v2 $fasta $result

class SurfaceArea:
    def __init__(self, probe_radius: float = 1.4, n_points: int = 100, radii_dict: Optional[Dict] = None):
        if radii_dict is None:
            radii_dict = {'X': 2.0}

        self.parser = PDBParser(QUIET=1)
        self.structure_computer = SASA.ShrakeRupley(probe_radius=probe_radius, n_points=n_points, radii_dict=radii_dict)

    def __call__(self, name, loc) -> float:
        struct = self.parser.get_structure(name, loc)
        # self.structure_computer.compute(struct, level="C")
        # return struct[0]['C'].sasa
        self.structure_computer.compute(struct, level="S")
        return struct.sasa


class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser()
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, chain_letters, overwrite=False, struct=None):

        pdb_id = self.out_dir.split('\\')[-1]
        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, self.out_dir)
        self.writer.set_structure(struct)
        self.writer.save(pdb_path, select=SelectChains(chain_letters))

        return pdb_path


class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)


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


def get_json_name(path):
    name_list = os.listdir(path)
    pdb_name = ''
    for name in name_list:
        # if name.startswith('test_9216f_unrelaxed_rank_1') and name[-4:] == '.pdb':
        if 'unrelaxed_rank_1' in name and name[-5:] == '.json':
            pdb_name = name
    if pdb_name == '':
        return ValueError
    else:
        return pdb_name


def get_sasa_difference(peptide_path, akt_pdb_path,result_path,peptide_chain,AKT_chain):
    # TODO replace the path
    path = result_path
    complex_pdb_name = get_complex_pdb_name(result_path)
    # print(complex_pdb_name)
    # peptide_pdb_name = 'peptide.pdb'
    # akt_pdb_name = 'akt.pdb'
    # splitter = ChainSplitter(os.path.join(path , complex_pdb_name))
    # splitter.make_pdb(os.path.join(path , peptide_pdb_name), peptide_chain)
    # splitter.make_pdb(os.path.join(path, akt_pdb_name), AKT_chain)
    sasa_fn = SurfaceArea()
    complex_sasa = sasa_fn(uuid.uuid4().hex,  os.path.join(path , complex_pdb_name))
    independent_sasa = sasa_fn(uuid.uuid4().hex, peptide_path)
    akt_sasa = sasa_fn(uuid.uuid4().hex, akt_pdb_path)
    # independent_sasa = sasa_fn(uuid.uuid4().hex, os.path.join(path , peptide_pdb_name))
    # akt_sasa = sasa_fn(uuid.uuid4().hex, os.path.join(path, akt_pdb_name))
    # akt_sasa = sasa_fn(uuid.uuid4().hex, path + akt_pdb_name)
    # print(akt_sasa, independent_sasa, complex_sasa)
    return akt_sasa + independent_sasa - complex_sasa


import json
import numpy as np
def get_confident_score(result_path):
    json_name = get_json_name(result_path)
    jfile = json.load(open(result_path + '/' + json_name, 'r'))
    plddt_list = jfile['plddt']
    score = np.mean(plddt_list)
    return score


import time
def run_colabfold_5(cmd):
    state=1
    for i in range(our_settings.try_num):
        if state!=0:
            try:
                state=subprocess.run(cmd,shell=True).returncode
            except:
                time.sleep(our_settings.sleep_time)
                continue
        else:
            break
    return state



from pathlib import Path
import sys
import hydra
import our_settings
class Colabfold_SASA():
    def __init__(self):
        # self.peptide_seq = peptide_seq
        # self.path=r'/home/qin/source_code/lambo_osc/lambo-main/colabfold_SASA'
        self.akt_seq = our_settings.AKT
        # self.uuid=peptide_seq
        self.project_root = hydra.utils.get_original_cwd()
        # self.project_root = Path(os.getcwd()).parents[3]
        self.path = os.path.join(self.project_root, 'colabfold_SASA')
        self.tokenizer=hydra.utils.instantiate({'_target_': 'lambo.utils.ResidueTokenizer'})
        self.work_dir=os.path.join(self.project_root,'data','experiments','test')

    def binding_ratio(self, pdbfile, anchor_before_len=0, anchor_after_len=0, peptide_chain='B', AKT_chain='C',
                      AKT_domain_start=149,
                      AKT_domain_end=409, CA_only=False):
        '''
        pdbfile: pdbfile position for AKT and peptide complex.
        in result pdb file, peptide and anchor_before,anchor_after are concatenate
        peptide_len: the length of original peptide
        anchor_before_len: the length of anchor_before
        anchor_after_len: the length of anchor_after
        peptide_chain indicate the peptide chain_AKT
        AKT_chain indicate the AKT chain
        AKT_domain_start: domain start position on AKT
        AKT_domain_end: domain end position on AKT
        CA_only: whether consider all atoms or CA only, default is False same with current ppt
        '''
        parser = PDBParser()
        structure = parser.get_structure('peptide_AKT', pdbfile)

        model = structure[0]
        chain_peptide = model[peptide_chain]
        chain_AKT = model[AKT_chain]

        AKT_domain_start=our_settings.AKT_domain_start
        AKT_domain_end=our_settings.AKT_domain_end

        peptide_len = len(chain_peptide) - anchor_before_len - anchor_after_len
        bind_num = 0
        cutoff = float(our_settings.cutoff)
        peptide_bind_sites = {}
        AKT_bind_sites = {}
        AKT_bind_sites_domain = {}
        for residue1 in chain_peptide:
            peptide_pos = residue1.id[1]
            if anchor_before_len != 0:
                if peptide_pos <= anchor_before_len:  # peptide_pos in anchor_before
                    continue

            if anchor_after_len != 0:
                if peptide_pos > anchor_before_len + peptide_len:  # peptide_pos in anchor_after
                    break

            # print(peptide_pos)
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
        # print(AKT_bind_sites)
        # print(peptide_bind_sites)
        # print(AKT_bind_sites_domain)
        return num_of_AKT_bind_sites_domain / (AKT_domain_end - AKT_domain_start + 1), num_of_AKT_bind_sites / len(
            chain_AKT), num_of_peptide_bind_sites / peptide_len

    def get_metrics(self,peptide_path, akt_pdb_path, result_path, peptide_chain, AKT_chain, single_chain, complex_path, metric):
        metrics = {}
        if metric == 'SASA'or metrics == 'confident':
            metrics[metric] = get_sasa_difference(peptide_path, akt_pdb_path,result_path, peptide_chain, AKT_chain)
            if metric == 'confident':
                metrics[metric] = get_confident_score(result_path)
        elif metric == 'energy':
            if os.path.exists(self.work_dir):
                metrics[metric] = FoldedCandidate(self.work_dir, peptide_path, [], self.tokenizer,
                                              skip_minimization=True, chain=single_chain,
                                              wild_name='test').mutant_total_energy
            else:
                os.system('mkdir '+self.work_dir)
                metrics[metric] = FoldedCandidate(self.work_dir, peptide_path, [], self.tokenizer,
                                                  skip_minimization=True, chain=single_chain,
                                                  wild_name='test').mutant_total_energy
        elif metric == 'ratio_AKT_domain':
            metrics[metric] = self.binding_ratio(pdbfile=complex_path)[0]*our_settings.ratio_weight
        elif metric == 'ratio_AKT':
            metrics[metric] = self.binding_ratio(pdbfile=complex_path)[1]*our_settings.ratio_weight
        elif metric == 'ratio_peptide':
            metrics[metric] = self.binding_ratio(pdbfile=complex_path)[2]*our_settings.ratio_weight
        return metrics

    def run_colab_sasa(self,  peptide_seq,metric):
        complex_seq = '>'+peptide_seq+'\n'+our_settings.anchor_before+peptide_seq +our_settings.anchor_after+':'+ self.akt_seq
        peptide_fa = '>' + peptide_seq + '\n' + our_settings.anchor_before + peptide_seq + our_settings.anchor_after
        akt_fa = '>akt' + '\n' + self.akt_seq

        # if our_settings.anchor_after=='' and our_settings.anchor_before=='':
        #     peptide_chain = 'B'
        #     AKT_chain = 'C'
        # elif our_settings.anchor_after!='' and our_settings.anchor_before!='':
        #     peptide_chain = 'C'
        #     AKT_chain = 'E'
        peptide_chain = 'B'
        AKT_chain = 'C'
        single_chain = 'A'

        result_path=os.path.join(self.path,'result',peptide_seq)
        peptide_result_path = os.path.join(self.path, 'pep_result', peptide_seq)
        akt_result_path = os.path.join(self.path, 'akt/')
        fa_path=os.path.join(self.path, 'fasta',peptide_seq+ '.fa')
        pep_fa_path = os.path.join(self.path, 'pep_fasta', peptide_seq + '.fa')
        akt_fa_path = os.path.join(self.path, 'akt', 'akt.fa')
        # print(fa_path)
        # peptide_path=os.path.join(self.path, 'result',peptide_seq,'peptide.pdb')
        peptide_pdb_path = os.path.join(self.path, 'pep_result', peptide_seq, 'peptide.pdb')
        akt_pdb_path = os.path.join(self.path, 'akt', 'akt.pdb')
        if os.path.exists(fa_path) and os.path.exists(result_path) and os.path.exists(peptide_pdb_path) \
                and os.path.exists(akt_pdb_path):
            complex_path=os.path.join(self.path, 'result',peptide_seq ,get_complex_pdb_name(result_path))
            metric=self.get_metrics(peptide_pdb_path,akt_pdb_path,result_path,peptide_chain,AKT_chain,single_chain,complex_path,metric)
            return metric
        else:
            # os.system('rm -rf ' + result_path)
            # os.system('mkdir ' + result_path)
            # os.system('touch '+fa_path)
            #echo "123456" | sudo -S rm -rf $result
            # os.system('echo "123456" | sudo -S rm -rf '+result_path)
            # os.system('echo "123456" | sudo -S touch ' + fa_path)
            print(complex_seq)

            with open(fa_path, 'w') as f:
                f.write(complex_seq)
            with open(pep_fa_path, 'w') as f:
                f.write(peptide_fa)

            # If don't have an environment that integrates colabfold and lambo, use a bash script to run colabfold
            # os.environ["CUDA_VISIBLE_DEVICES"] = "0"
            sh_path = os.path.join(self.path, 'run_complex_osc.sh')
            sh_single_path = os.path.join(self.path, 'run_single_osc.sh')
            cmd = 'bash -i '+sh_path+' '+fa_path+' '+result_path
            peptide_cmd = 'bash -i ' + sh_single_path + ' ' + pep_fa_path + ' ' + peptide_result_path

            if not os.path.exists(akt_pdb_path):
                with open(akt_fa_path, 'w') as a:
                    a.write(akt_fa)
                akt_cmd = 'bash -i ' + sh_single_path + ' ' + akt_fa_path + ' ' + akt_result_path
                state_akt = run_colabfold_5(akt_cmd)
                if state_akt != 0:
                    raise NameError('Colabfold have not run successfully!')
                colab_akt_name = get_complex_pdb_name(akt_result_path)
                with open(akt_pdb_path, 'w') as aw:
                    with open(os.path.join(akt_result_path, colab_akt_name), 'r') as ar:
                        aw.writelines(ar.readlines())

            #If the colabfold can run in terminal, you can run them directly with the following command
            # try:
            #     os.system('colabfold_batch', shell=True)
            # except:
            #     print('colabfold_batch have not add to environment path!')
            #     sys.exit(0)
            # os.system('mkdir ' + result_path)
            # cmd='colabfold_batch --use-gpu-relax --num-recycle 3 --model-type AlphaFold2-multimer-v2 '+fa_path+' '+result_path

            state=run_colabfold_5(cmd)
            state_pep = run_colabfold_5(peptide_cmd)
            colab_pep_name = get_complex_pdb_name(peptide_result_path)
            # colab_akt_name = get_complex_pdb_name(akt_result_path)
            with open(peptide_pdb_path, 'w') as pw:
                with open(os.path.join(peptide_result_path, colab_pep_name), 'r') as pr:
                    pw.writelines(pr.readlines())

            complex_path = os.path.join(self.path, 'result', peptide_seq, get_complex_pdb_name(result_path))
            if state == 0 and state_pep == 0:
                return self.get_metrics(peptide_pdb_path,akt_pdb_path,result_path,peptide_chain,AKT_chain,single_chain,complex_path,metric)
            else:
                raise NameError('Colabfold have not run successfully!')

