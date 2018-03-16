import os
os.chdir("P:/1_Projects/Adarsh/1_gwasQC")

import numpy as np
import pandas as pd

d = {}
#snpMISS_filename = "./outputs/QC_train_t321-0105/clean_inds_data_missing.lmiss"
snpMISS_filename = "./outputs/QC_train_t321-0106/clean_inds_data_missing.lmiss"
with open(snpMISS_filename, 'r') as f:
    next(f)
    for line in f:
        words = line.split()
        key = words[1]
        value = float(words[4])
        d[key] = value

dupSNP = []
#snpDUP_filename = "./outputs/QC_train_t321-0105/duplicatedSNPs.dupvar"
snpDUP_filename = "./outputs/QC_train_t321-0106/duplicatedSNPs.dupvar"
with open(snpDUP_filename, 'r') as f:
    for line in f:
        words = line.split()
        missing_rates = [(i, d[i]) for i in words]
        missing_rates_sorted = sorted(missing_rates, key=lambda x: x[1])
        for j in range(1, len(missing_rates_sorted)):
            dupSNP.append(missing_rates_sorted[j][0])
#with open("./outputs/QC_train_t321-0105/duplicatedSNPs.txt",'w') as file:
with open("./outputs/QC_train_t321-0106/duplicatedSNPs.txt", 'w') as file:
    for item in dupSNP:
        file.write(item + "\n")
