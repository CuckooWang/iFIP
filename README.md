# iFIP
Inference of functional interacting partners

## The description of each source code
### KnownAutophagy.py
For genes interacted with Atg9, calculate the number of interacted known Atg proteins and autophagy regulators from THANATOS (http://thanatos.biocuckoo.org/).
### GetRawScore.py
For genes interacted with Atg9, get the raw scores from mRNA, protein and PPI levels for further training.
### LGTraining.py
Training of the logistic regression (LG) model.
### GSEATraining.py
Optimize the LG model using gene set enrichment analysis (GSEA).
### demo
A small dataset to demo above codes, an example of the output of training is also provided.

## Software Requirements
### OS Requirements
Above codes have been tested on the following systems:  
Windows: Windows 7, Windos 10  
Linux: CentOS linux 7.8.2003  
### Hardware Requirements
All codes and softwares could run on a "normal" desktop computer, no non-standard hardware is needed

## Installation guide
All codes can run directly on a "normal" computer with Python 3.7.9 installed, no extra installation is required

## Additional information
Expected run time for each program is about 5 seconds.
## Contact
Dr. Yu Xue: xueyu@hust.edu.cn  
Dr. Chenwei Wang: wangchenwei@hust.edu.cn  
Dr. Di Peng: pengdi@hust.edu.cn 
