# COVID-19 pandemic re-shapes the global dispersal patterns of seasonal influenza viruses

Zhiyuan Chen<sup>1,2</sup>, Joseph L.-H. Tsui<sup>2</sup>, Bernardo Gutierrez<sup>2,3</sup>, Simon Busch Moreno<sup>2</sup>, Louis du Plessis<sup>4,5</sup>, Xiaowei Deng<sup>6</sup>, Jun Cai<sup>1</sup>, Sumali Bajaj<sup>2</sup>, Marc A. Suchard<sup>7</sup>, Oliver G. Pybus<sup>2,8,‡</sup>, Philippe Lemey<sup>9,‡,†</sup>, Moritz U. G. Kraemer<sup>2,10,‡,†</sup>, Hongjie Yu<sup>1,‡,†</sup>

1.	Department of Epidemiology, School of Public Health, Key Laboratory of Public Health Safety, Ministry of Education, Fudan University, Shanghai, China
2.	Department of Biology, University of Oxford, Oxford, UK
3.	Colegio de Ciencias Biologicas y Ambientales, Universidad San Francisco de Quito USFQ, Quito, Ecuador
4.	Department of Biosystems Science and Engineering, ETH Zürich, Basel, Switzerland
5.	Swiss Institute of Bioinformatics, Lausanne, Switzerland
6.	Department of Epidemiology, National Vaccine Innovation Platform, School of Public Health, Nanjing Medical University, Nanjing, China
7.	Departments of Biostatistics, Biomathematics and Human Genetics, University of California, Los Angeles, Los Angeles, CA, USA
8.	Department of Pathobiology and Population Sciences, Royal Veterinary College, London, UK
9.	Department of Microbiology, Immunology and Transplantation, Rega Institute, KU Leuven, Leuven, Belgium
10.	Pandemic Sciences Institute, University of Oxford, Oxford, UK

<sup>‡</sup> Contributed equally as senior authors  
<sup>†</sup> Correspondence: moritz.kraemer@biology.ox.ac.uk; philippe.lemey@kuleuven.be; yhj@fudan.edu.cn

## Before starting
1. This repository contains code and data used in the paper entitled "COVID-19 pandemic re-shapes the global dispersal patterns of seasonal influenza viruses"  
2. All code and data contained within this repository is released under the CC BY-NC-SA License. 
3. Genetic data has not been provided due to the policy restriction, and can be accessed in NCBI (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and GISAID (https://www.gisaid.org/). 
4. Map files have not been provided due to redistribution restriction, and can be accessed in GADM (https://gadm.org/data.html).  
5. Due to the large sizes of Markov jumps files used in the paper, we have reduced the file size by including few tree files as an example (complete data has been deposited in https://figshare.com/account/home#/projects/187860).

## Repository structure and usage
The structure of this repository is shown below.  
Data are deposited in "data_part" (except for large files, see above), codes and results for genomics-based analyses are deposited in "genomic_part", codes for the regression model are deposited in "model_part", and codes for generating the main figures can be found in "figure_scripts".

```
global_influenza_project/
├── data_part
│   ├── air_traffic_data
│   ├── epi_data
│   ├── map_data
├── genomic_part
│   ├── phylogenetic_analyses
│   ├── phylogeographic_analyses
│   │   ├── step1
│   │   └── step2
│   ├── post-analyses
│   │   ├── glm_log_file
│   │   ├── trunk
│   │   ├── jump_history
│   │   ├── persistence
│   │   ├── mcc_tree
│   │   ├── diversity
│   │   ├── pop_size
│   │   └── nature_selection
│   └── acknowledge_table
├── model_part
│   ├── bayesian_model.py
│   ├── regression_plots
│   └── forest_plots
├── figure_scripts
│   ├── Fig1.r
│   ├── Fig2.r
│   ├── Fig3.r
│   ├── Fig4.r
│   ├── Fig5.r
│   └── output
└── README.md
```

<h1> License </h1>
<h4>CC BY-NC-SA 4.0 </h4>

Attribution-NonCommercial-ShareAlike 4.0 International
This license requires that reusers give credit to the creator. It allows reusers to distribute, remix, adapt, and build upon the material in any medium or format, for noncommercial purposes only. If others modify or adapt the material, they must license the modified material under identical terms.

BY: Credit must be given to you, the creator.

NC: Only noncommercial use of your work is permitted.
Noncommercial means not primarily intended for or directed towards commercial advantage or monetary compensation.

SA: Adaptations must be shared under the same terms.

COVID-19 pandemic re-shapes the global dispersal patterns of seasonal influenza viruses © 2023 by Zhiyuan Chen is licensed under CC BY-NC-SA 4.0 

Copyright (c) 2023 Zhiyuan Chen, Fudan University & University of Oxford
