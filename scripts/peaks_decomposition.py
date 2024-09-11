"""
File: peaks_decomposition.py
Author: Oscar Laurent
Date: 2024-09-11
Description: A Python script to perform the peaks decomposition of raw ToF-SIMS spectra into a cvs file
"""

# ===================================================================
#  MODULE IMPORTATION
# ===================================================================
from utilities.molecule import Molecule, Element
from utilities.periodic_table import get_periodic_table
from utilities.SIMS_Spectra_class import SIMS_Spectra
import numpy as np
import pandas as pd
import os
import shutil
import re
from scipy.ndimage import gaussian_filter1d
import plotly.graph_objects as go
from pybaselines import Baseline
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, ScalarFormatter)


# ===================================================================
# - DATA RETRIEVAL
# ===================================================================
# Scan over the data/ToF-SIMS_data/raw_spectra folder to get a list of path of all the spectra

positive_filename_list = []
negative_filename_list = []
directory_name = "data/ToF-SIMS_data/raw_spectra"
for filename in os.listdir(directory_name):
    if re.match("P_", filename):
        positive_filename_list.append(os.path.join(directory_name, filename))
    elif re.match("N_", filename):
        negative_filename_list.append(os.path.join(directory_name, filename))
positive_filename_list.sort(reverse=True)
negative_filename_list.sort(reverse=True)


# ===================================================================
# - POSITIVE DATAFRAME CREATION 
# ===================================================================
single_element = ["H", "K", "Na", "Cu", "Zn", "Al", "CH3", "C2H5"]
other_elements = []

# Creation of a dict holding all the isotopes lists and isotope dfs to avoid to compute them each time 
isotope_df_dict = {"Ru":None,
                   "Ru-Pt":None,
                   "Ru-Pd":None,
                   "Ru-Ir":None,
                   "Ru-Rh":None,
                   "Ru-Pt-Pd":None,
                   "Ru-Pt-Pd-Ir":None,
                   "Ru-Pt-Pd-Ir-Rh": None,
                   "HEA19": None,
                   "HEA 1.2": None}

df_all = pd.DataFrame()
for path in positive_filename_list:
    sample_name_index = np.where(["_"+x+"_" in path for x in isotope_df_dict.keys()])[0].item()
    sample_name = list(isotope_df_dict.keys())[sample_name_index]
    
    if sample_name == "HEA19" or sample_name == "HEA 1.2":
        sample_name = "Ru-Pt-Pd-Ir-Rh"
    sample_components = sample_name.split("-")

    if isotope_df_dict.get(sample_name) is None :
        spectra = SIMS_Spectra(path,
                                components = sample_components, 
                                single_other_component= single_element,
                                other_elements=other_elements,
                                indice_max=5,
                                max_molecules=4,
                                highest_peak_mz=1800)
        df = spectra.get_output_df()
        isotope_df_dict[sample_name] = (spectra.isotopes, spectra.isotopes_df)
    
    else : 
        sample_isotopes = isotope_df_dict.get(sample_name)[0]
        sample_isotope_df = isotope_df_dict.get(sample_name)[1]
        spectra = SIMS_Spectra(path,
                               single_other_component= single_element,
                               other_elements= other_elements,
                               indice_max=5,
                               max_molecules=4,
                               isotopes_df=sample_isotope_df,
                               isotopes=sample_isotopes,
                               highest_peak_mz=1800)
        df = spectra.get_output_df()
    
    df_all = pd.concat([df_all, df], ignore_index=True)
df_all.fillna(0, inplace=True)

# Remove Columns where the quantification value is 0 for all the samples
df_all_cleaned = df_all.loc[:,(df_all !=0).any().values]
df_all_cleaned.to_csv("data/ToF-SIMS_data/processed_data/AllSample_PositiveMode_ProcessedData.csv", index=False)


# ===================================================================
# - NEGATIVE DATAFRAME CREATION 
# ===================================================================

single_element = ["H", "O", "OH", "Cl", "Cu", "Zn", "Al", "CH3", "C2H"]
other_elements = ["O"]

# Creation of a dict holding all the isotopes lists and isotope dfs to avoid to compute them each time 
isotope_df_dict = {"Ru":None,
                   "Ru-Pt":None,
                   "Ru-Pd":None,
                   "Ru-Ir":None,
                   "Ru-Rh":None,
                   "Ru-Pt-Pd":None,
                   "Ru-Pt-Pd-Ir":None,
                   "Ru-Pt-Pd-Ir-Rh": None,
                   "HEA19": None,
                   "HEA 1.2": None}
df_all = pd.DataFrame()

for path in negative_filename_list:
    sample_name_index = np.where(["_"+x+"_" in path for x in isotope_df_dict.keys()])[0].item()
    sample_name = list(isotope_df_dict.keys())[sample_name_index]

    if isotope_df_dict.get(sample_name) is None :
        if sample_name == "HEA19" or sample_name == "HEA 1.2":
            sample_components = ["Ru", "Rh", "Pd", "Pt", "Ir"]
        else : 
            sample_components = sample_name.split("-")

        spectra = SIMS_Spectra(path,
                                components = sample_components, 
                                single_other_component= single_element,
                                other_elements=other_elements,
                                indice_max=4,
                                max_molecules=2)
        df = spectra.get_output_df()
        isotope_df_dict[sample_name] = (spectra.isotopes, spectra.isotopes_df)
    
    else : 
        sample_isotopes = isotope_df_dict.get(sample_name)[0]
        sample_isotope_df = isotope_df_dict.get(sample_name)[1]
        spectra = SIMS_Spectra(path,
                               single_other_component= single_element,
                               other_elements= other_elements,
                               indice_max=4,
                               max_molecules=2,
                               isotopes_df=sample_isotope_df,
                               isotopes=sample_isotopes)
        df = spectra.get_output_df()
    
    df_all = pd.concat([df_all, df], ignore_index=True)
df_all.fillna(0, inplace=True)

# Remove Columns where the quantification value is 0 for all the samples
df_all_cleaned = df_all.loc[:,(df_all !=0).any().values]
#df_all_cleaned.to_csv("data/ToF-SIMS_data/processed_data/AllSample_NegativeMode_ProcessedData.csv", index=False)
