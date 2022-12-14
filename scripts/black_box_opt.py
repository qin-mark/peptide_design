import hydra
import wandb
import warnings
import random
import logging
import os
from pathlib import Path
import sys
import torch

from omegaconf import OmegaConf

from upcycle.scripting import startup
from upcycle.logging.analysis import flatten_config

from gpytorch.settings import max_cholesky_size
# os.environ["CUDA_VISIBLE_DEVICES"]="-1"
import time
import our_settings


@hydra.main(config_path='../hydra_config', config_name='black_box_opt')
def main(config):
    project_path = hydra.utils.get_original_cwd()
    our_setting_conf = {
        'RESUME': our_settings.RESUME,
        'AKT': our_settings.AKT,
        'time_stamp': our_settings.time_stamp,
        'anchor_before': our_settings.anchor_before,
        'anchor_after': our_settings.anchor_after,
        'scoring_metrics': our_settings.scoring_metrics,
        'ratio_weight': our_settings.ratio_weight,
        'sleep_time': our_settings.sleep_time,
        'try_num': our_settings.try_num,
        'cutoff': our_settings.cutoff,
        'AKT_domain_start': our_settings.AKT_domain_start,
        'AKT_domain_end': our_settings.AKT_domain_end,
    }

    if our_setting_conf['time_stamp'] == '':
        our_setting_conf['time_stamp'] = str(time.time()).split('.')[0]
        os.mkdir(os.path.join(project_path, 'data', 'experiments', our_setting_conf['time_stamp']))
        # setting_conf_path = os.path.join(project_path, 'data', 'experiments', 'test',
        #                                  our_setting_conf['time_stamp'],
        #                                  'setting_config.pt')
        with open('../work_time_stamp.txt', 'a') as ti:
            ti.write('work name: ' + our_setting_conf['time_stamp'] + '\n' + 'config: ' + '\n' + '\t')
            config_all = ''
            for key, value in our_setting_conf.items():
                config_all = config_all + str(key) + ': ' + str(value) + '\n' + '\t'
            ti.write(config_all + '\n')
    #judge seetings of init lambo
    # if our_settings.RESUME == True:
    #     path_config = os.path.join(os.path.join(project_path, 'data', 'experiments', 'test', 'config.pt'))
    #     if os.path.exists(path_config):
    #         assert os.path.exists(path_config), 'you have not run once at least, please set RESUME=False in oursettings.py'
    #         config_RESUME = torch.load(path_config)
    #         for i in config_RESUME.keys():
    #             if i == 'logger' or i == 'timestamp':
    #                 continue
    #             else:
    #                 assert config_RESUME[i] == config[i] ,'Please do not update configuration before you start breakpoint optimization!'
    #
    # torch.save(config, os.path.join(project_path, 'data', 'experiments', 'test', 'config.pt'))

    # setup
    # random.seed(None)  # make sure random seed resets between multirun jobs for random job-name generation
    log_config = flatten_config(OmegaConf.to_container(config, resolve=True), sep='/')
    log_config = {'/'.join(('config', key)): val for key, val in log_config.items()}
    wandb.init(project='lambo', config=log_config, mode=config.wandb_mode,
               group=config.exp_name)
    config['job_name'] = wandb.run.name
    config, _ = startup(config)  # random seed is fixed here

    # if torch.cuda.is_available():
    #     torch.backends.cudnn.benchmark = True

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # catch really annoying BioPython warnings

        try:
            # create initial candidates, dataset
            tokenizer = hydra.utils.instantiate(config.tokenizer)
            bb_task = hydra.utils.instantiate(config.task, tokenizer=tokenizer, candidate_pool=[])

            project_root = Path(os.getcwd()).parents[2]  # changing the Hydra run dir will break this.
            # base_candidates, base_targets, all_seqs, all_targets = bb_task.task_setup(config, project_root=project_root)

            # path_checkpoint = os.path.join(os.path.join(project_root, 'data', 'experiments', 'test', 'temp_data.pt'))
            #
            # # config.surrogate.num_epochs = len(our_settings.scoring_metrics) * 128
            # if our_settings.RESUME==True and os.path.exists(path_checkpoint):                 # pool
            #     #TODO:hold on training
            #     #load base_candidates, base_targets, all_seqs, all_targets
            #
            #     # path_checkpoint = os.path.join(os.path.join(project_root, 'data', 'experiments', 'test','temp_data.pt'))
            #     assert os.path.exists(path_checkpoint),'you have not run once, please set RESUME=False in oursettings.py'
            #     checkpoint = torch.load(path_checkpoint)     #active_candidates active_targets   active_seqs
            #     base_candidates=checkpoint['base_candidate']
            #     base_targets=checkpoint['base_target']
            #     active_candidates=checkpoint['active_candidates']  #pool
            #     active_targets=checkpoint['active_targets']
            #     active_seqs=checkpoint['active_seqs']
            #     round_idx=checkpoint['round_start']
            #
            #     max_chol_sz = config.surrogate.get('max_cholesky_size', int(1e5))
            #     with max_cholesky_size(max_chol_sz):
            #         optimizer = hydra.utils.instantiate(
            #             config.optimizer,
            #             bb_task=config.task,
            #             surrogate=config.surrogate,
            #             acquisition=config.acquisition,
            #             encoder=config.encoder,
            #             tokenizer=tokenizer
            #         )
            #         metrics = optimizer.optimize(
            #             base_candidates, base_targets, active_seqs,active_targets,active_candidates,round_idx, log_prefix=config.task.log_prefix
            #         )
            # else:
            base_candidates, base_targets, all_seqs, all_targets = bb_task.task_setup_SASA_energy(config,project_root=project_root)
            # optimizer
            max_chol_sz = config.surrogate.get('max_cholesky_size', int(1e5))
            with max_cholesky_size(max_chol_sz):
                optimizer = hydra.utils.instantiate(
                    config.optimizer,
                    bb_task=config.task,
                    surrogate=config.surrogate,
                    acquisition=config.acquisition,
                    encoder=config.encoder,
                    tokenizer=tokenizer,
                    setting_conf=our_setting_conf,
                )
                metrics = optimizer.optimize(
                    base_candidates, base_targets, all_seqs, all_targets,None,1, log_prefix=config.task.log_prefix
                )

            metrics = {key.split('/')[-1]: val for key, val in metrics.items()}  # strip prefix
            ret_val = metrics['hypervol_rel']

        except Exception as err:
            logging.exception(err)
            ret_val = float('NaN')

    wandb.finish()  # necessary to log Hydra multirun output to different jobs
    return ret_val


if __name__ == "__main__":
    main()
    sys.exit()
