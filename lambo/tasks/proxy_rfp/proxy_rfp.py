import os

import hydra
import numpy as np
import pandas as pd
import torch
from botorch.utils.multi_objective import pareto
import math

import our_settings
from lambo.candidate import FoldedCandidate
from lambo.tasks.base_task import BaseTask
from colabfold_SASA.Colabfold_SASA import  Colabfold_SASA

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

class ProxyRFPTask(BaseTask):
    def __init__(self, tokenizer, candidate_pool, obj_dim, transform=lambda x: x,
                 num_start_examples=1024, **kwargs):
        super().__init__(tokenizer, candidate_pool, obj_dim, transform, **kwargs)
        self.op_types = ["sub"]
        # self.num_start_examples = num_start_examples
        self.num_start_examples = 3 #candidate sample numbers
        self.get_metrics=Colabfold_SASA()

    def task_setup(self, config, project_root=None, *args, **kwargs):
        project_root = hydra.utils.get_original_cwd() if project_root is None else project_root
        work_dir = f'{project_root}/{config.log_dir}/{config.job_name}/{config.timestamp}/foldx'
        rfp_known_structures = pd.read_csv(
            f'{project_root}/lambo/assets/fpbase/rfp_known_structures.csv'
        )

        all_seqs = rfp_known_structures.foldx_seq.values
        all_targets = np.stack([
            -rfp_known_structures.SASA.values,
            rfp_known_structures.foldx_total_energy.values,
        ], axis=-1)

        seed_data = pd.read_csv(
            f'{project_root}/lambo/assets/fpbase/proxy_rfp_seed_data.csv'
        )
        seed_data = seed_data.sample(self.num_start_examples)
        sample_batch_targets = np.stack([
            -seed_data.SASA.values,
            -seed_data.stability.values,
        ], axis=-1)

        all_seqs = np.concatenate((all_seqs, seed_data.foldx_seq.values))
        all_targets = np.concatenate((all_targets, sample_batch_targets))

        seq_len_mask = np.array([len(x) <= self.max_len for x in all_seqs])
        all_seqs = all_seqs[seq_len_mask]
        all_targets = all_targets[seq_len_mask]

        # filter candidate sequences by length
        foldx_seq_len = rfp_known_structures.foldx_seq.apply(lambda x: len(x))
        rfp_known_structures = rfp_known_structures[foldx_seq_len <= self.max_len]
        rfp_known_structures.reset_index(inplace=True)

        # find valid, non-dominated starting candidates
        valid_targets = np.stack([
            -rfp_known_structures.SASA.values,
            rfp_known_structures.foldx_total_energy.values,
        ], axis=-1)
        pareto_mask = pareto.is_non_dominated(-torch.tensor(valid_targets))
        base_targets = valid_targets[pareto_mask]

        base_candidates = []
        for row_idx, datum in rfp_known_structures.iterrows():
            if not pareto_mask[row_idx]:
                continue
            print(f'{datum.Name} is non-dominated, adding to start pool')
            pdb_id = datum.pdb_id.lower()
            chain_id = datum.longest_chain
            parent_pdb_path = f'{project_root}/lambo/assets/foldx/{pdb_id}_{chain_id}/wt_input_Repair.pdb'
            base_candidates.append(
                FoldedCandidate(work_dir, parent_pdb_path, [], self.tokenizer,
                                skip_minimization=True, chain=chain_id, wild_name=datum.Name)
            )
        base_candidates = np.array(base_candidates).reshape(-1)

        return base_candidates, base_targets, all_seqs, all_targets

    def task_setup_SASA_energy(self, config, project_root=None, *args, **kwargs):
        project_root = hydra.utils.get_original_cwd() if project_root is None else project_root
        work_dir = f'{project_root}/{config.log_dir}/{config.job_name}/{config.timestamp}/foldx'

        # rfp_known_structures = pd.read_csv(r'/users/PCON0022/yuy702/qwy/lambo-main/lambo/assets/fpbase/rfp_known_structures_peptide.csv')
        # rfp_known_structures = pd.read_csv(r'/home/qin/source_code/lambo-main/lambo/assets/fpbase/rfp_known_structures_peptide.csv')
        # rfp_known_structures = pd.read_csv(r'rfp_known_structures_peptide.csv')
        rfp_known_structures = pd.read_csv(f'{project_root}/lambo/assets/fpbase/rfp_known_structures_peptide.csv')

        peptide_chain = 'B'
        AKT_chain = 'C'

        base_targets = []
        all_targets = []
        all_seqs = []
        base_seqs = []
        base_candidates = []

        for row_idx, datum in rfp_known_structures.iterrows():
            metrics = {}
            for metric in our_settings.scoring_metrics:
                assert metric in ['SASA', 'energy', 'ratio_AKT', 'ratio_AKT_domain','ratio_peptide'],'metric in scoring funtion is wroung, please check the oursetting\'s metric spell'
                if math.isnan(datum[metric]):
                    metric=self.get_metrics.run_colab_sasa(datum.fpbase_seq,metric)
                    metrics.update(metric)
                else:
                    metric={metric:datum[metric]}
                    metrics.update(metric)

            peptide_path = os.path.join(f'{project_root}',
                                        f'colabfold_SASA/result/' + datum.fpbase_seq + '/peptide.pdb')
            print(datum.pdb_id)
            fold_candidate=FoldedCandidate(work_dir, peptide_path, [], self.tokenizer,
                            skip_minimization=True, chain=peptide_chain, wild_name=datum.Name)

            metrics_list=[]
            for metric_name in our_settings.scoring_metrics:
                metrics_list.append(metrics[metric_name])

            base_candidates.append(fold_candidate)
            base_targets.append(metrics_list)
            all_seqs.append(datum.fpbase_seq)
            all_targets.append(metrics_list)
            base_seqs.append(datum.fpbase_seq)

        base_candidates = np.array(base_candidates).reshape(-1)

        #load seeds data
        seed_data = pd.read_csv(f'{project_root}/lambo/assets/fpbase/proxy_rfp_seed_data_peptide.csv')

        for row_idx, datum in seed_data.iterrows():
            metrics = {}
            for metric in our_settings.scoring_metrics:
                if math.isnan(datum[metric]):
                    metric = self.get_metrics.run_colab_sasa(datum.foldx_seq, metric)
                    metrics.update(metric)
                else:
                    metric={metric:datum[metric]}
                    metrics.update(metric)

            assert len(metrics)==len(our_settings.scoring_metrics),'scoring funtion have some problem, please contact qwy'

            metrics_list = []
            for metric_name in our_settings.scoring_metrics:
                metrics_list.append(metrics[metric_name])

            all_targets.append(metrics_list)
            all_seqs.append(datum.foldx_seq)

        return base_candidates, np.array(base_targets), np.array(all_seqs), np.array(all_targets)


    def get_energy(self, config, project_root=None,bb_task=None, *args, **kwargs):
        all_seq = []
        project_root =r'/home/qin/source_code/lambo-source/lambo-main'#Temporary file storage path
        with open(os.path.join(project_root,r'lambo/assets/fpbase/proxy_rfp_seed_data_peptide.csv'), 'r') as f:
            lines = f.readlines()
            lines = lines[1:]
        for line in lines:
            all_seq.append(line.split(',')[0])
        for seq in all_seq:
            print(seq)
            result_path = r'/home/qin/source_code/lambo-source/lambo-main/colabfold_SASA/result'
            pdb_path = os.path.join(result_path, seq)
            pdb_path = os.path.join(result_path,seq,get_complex_pdb_name(pdb_path))
            energy_value=FoldedCandidate(project_root, pdb_path, [],
                            self.tokenizer,
                            skip_minimization=True, chain='B', wild_name='test').mutant_total_energy
            # bb_task.score([FoldedCandidate(project_root, pdb_path, [],
            #                                self.tokenizer,
            #                                skip_minimization=True, chain='B', wild_name='test')])

            # print(energy_value)


    def make_new_candidates(self, base_candidates, new_seqs):
        assert base_candidates.shape[0] == new_seqs.shape[0]
        new_candidates = []
        for b_cand, n_seq in zip(base_candidates, new_seqs):
            b_seq = b_cand.mutant_residue_seq
            assert len(b_seq) == len(n_seq), 'FoldX only accepts substitutions'
            mutation_ops = []
            for i, (b_char, n_char) in enumerate(zip(b_seq, n_seq)):
                if not b_char == n_char:
                    mutation_ops.append(b_cand.new_mutation(i, n_char, 'sub'))
            # mutation_ops = mutation_list(b_cand.mutant_residue_seq, n_seq)
            new_candidates.append(b_cand.new_candidate(mutation_ops))
        return np.stack(new_candidates)

    def _evaluate(self, x, out, *args, **kwargs):
        assert x.ndim == 2
        x_cands, x_seqs, f_vals = [], [], []
        for query_pt in x:
            cand_idx, mut_pos, mut_res_idx, _ = query_pt
            base_candidate = self.candidate_pool[cand_idx]
            mut_res = self.tokenizer.sampling_vocab[mut_res_idx]
            mut_list = [base_candidate.new_mutation(mut_pos, mut_res, mutation_type='sub')]
            candidate = base_candidate.new_candidate(mut_list)
            x_cands.append(candidate)
            x_seqs.append(candidate.mutant_residue_seq)
        x_seqs = np.array(x_seqs).reshape(-1)
        x_cands = np.array(x_cands).reshape(-1)

        out["X_cand"] = x_cands
        out["X_seq"] = x_seqs
        norm_scores = self.transform(self.score(x_cands))
        out["F"] = norm_scores

    def score(self, candidates):
        f_vals = []
        for cand in candidates:
            f_vals.append(np.array([
                -cand.mutant_surface_area,
                cand.mutant_total_energy,
            ]))
        return np.stack(f_vals)