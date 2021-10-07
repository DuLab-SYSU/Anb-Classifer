# Anb-classifier

**This is the source code and data of "Heavy Chain Sequence-Based Classifier for the Specificity of Human Antibodies".**

To run the trained classifier, run Anb_classifier.py in the bin folder. *** Please unzip the `stacking_model.zip` file first** (in `./bin/model` folder).

```
python3 Anb_classifier.py -i input path -o output path for the predictions --olga_path The absolute path of python script: compute_pgen.py --antigen_names The folder name of IMGT HighV-QUEST results files
```

## Args:

```
-i, IMGT HighV-QUEST results file dir.
  --input dir
     --result file1(name the folder as the antigen_name)
     --result file2(name the folder as the antigen_name)
     --...
-o, Output path
--remove_redandunt, default=True,
   True or False. It is faster when running on non-redandunt data.
--format, csv or tsv. Output file format of prediction. 
--cpu, Specify how many cores to use, default=1
--olga_path, The absolute path of python script: compute_pgen.py. /xxx/xx/olga/compute_pgen.py
--pssm_path, The absolute path of PSSM result files
--mode, validate or predict, default="predict"
--antigen_names, The folder names of IMGT HighV-QUEST result files, split by ',', default="hiv,flu,pps,acpa,tt,hcv,hbv"
```

## Requirements:

```
Linux environment
python 3.7
numpy
OLGA (https://github.com/statbiophys/OLGA)
sklearn
biopython
Result files from IMGT HighV-QUEST. If you want to validate this model, please make sure each result file is named as the antigen name(hiv, flu, pps, acpa, tt, hcv, and hbv). 
PSSM files from POSSUM webserver. Please make sure all PSSM files are put into one folder.
```

