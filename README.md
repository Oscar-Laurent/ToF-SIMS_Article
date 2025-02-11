# Supplementary Data for "High Entropy Alloys: Assessing Atomic-Scale Mixing and Surface Passivation with Time-of-Flight Secondary Ion Mass Spectrometry"

![image](https://github.com/user-attachments/assets/d884f7ee-d3bc-4f6d-a611-c4a3be7295eb)

## Repository Structure

The repository is organized as follows:
```
ToF-SIMS_Article/
│
├── README.md                        
├── data/                             # Contains all data files of the article               
│   ├── ToF-SIMS_data/              
│   ├── XRD_raw_spectra/              
│   ├── figures/                     
│   └── supplementary_material/       
│
├── scripts/                          # Scripts for data analysis and figure generation
│   ├── XRD_figure.py                 
│   ├── data_processing.R             
│   ├── figures_generation.py         
│   ├── peaks_decomposition.py        
│   └── utilities/                    
│
└── LICENSE   
```

## Usage
To process the raw data and generate figures, use the scripts provided in the `scripts/` folder. The workflow typically involves:

1.	Loading the data from the data/ directory.
2.	Running the respective scripts for preprocessing (data_processing.R) and analysis (e.g., peaks_decomposition.py).
3.	Generating figures using the figures_generation.py script.

## Contact 
For any questions or collaborations, please feel free to reach out to the authors:

* Oscar Laurent (UCLouvain) – oscar.laurent@uclouvain.be
* Arnaud Delcorte (UCLouvain) – arnaud.delcorte@uclouvain.be
* Damien P. Debecker (UCLouvain) – damien.debecker@uclouvain.be
