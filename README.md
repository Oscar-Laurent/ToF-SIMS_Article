# ToF-SIMS Article: High Entropy Alloys

This repository contains the data, scripts, and supplementary materials for the article titled **“High entropy alloys: assessing atomic-scale mixing and surface passivation with time-of-flight secondary ion mass spectrometry.”**
![image](https://github.com/user-attachments/assets/d884f7ee-d3bc-4f6d-a611-c4a3be7295eb)

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
