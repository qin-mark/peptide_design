# 1. Environment
## 1.1 Lambo environment
To install the lambo environment, please refer to the environment building section of the lambo tutorial. \
[https://github.com/samuelstanton/lambo]
## 1.2 Colabfold environment
Please refer to the colabfold installation tutorial.\
[https://github.com/sokrypton/ColabFold]
And then, Set the colabfold_batch command as a system environment variable
# 1.3 Running settings 
Before running, Please go to our_settings.py file under the project to set the running parameters.
# 2. Running command

```
python scripts/black_box_opt.py optimizer=lambo optimizer.encoder_obj=mlm task=proxy_rfp tokenizer=protein surrogate=multi_task_exact_gp acquisition=nehvi
```
# 3. Example
## 3.1 Input
For input, please refer to the lambo_document.docx document in the project path
## 3.2 Output
The project will output a 'log.csv' file, located under 'data/experiment/test'



