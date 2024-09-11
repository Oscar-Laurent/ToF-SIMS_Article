"""
File: XRD_figure.py
Author: Oscar Laurent
Date: 2024-09-11
Description: A Python script to generate the XRD figures in the article and in the supplementary information
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import math
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


# =======================================================
#  - COMBINED XRD PLOT
# =======================================================

# Theoretical phases creation
FCC_Rh = pd.DataFrame({"2theta":np.linspace(6,90, 2500), "Intensity":np.zeros(2500)}).set_index("2theta")
FCC_Rh.loc[np.logical_and(FCC_Rh.index < 41.2, FCC_Rh.index > 41.1), "Intensity"] = 1
FCC_Rh.loc[np.logical_and(FCC_Rh.index < 47.9, FCC_Rh.index > 47.8), "Intensity"] = 0.2
FCC_Rh.loc[np.logical_and(FCC_Rh.index < 70.1, FCC_Rh.index > 70.0), "Intensity"] = 0.2
FCC_Rh.loc[np.logical_and(FCC_Rh.index < 84.6, FCC_Rh.index > 84.5), "Intensity"] = 0.1
FCC_Rh["Intensity"] += 0.25

FCC_Ru = pd.DataFrame({"2theta":np.linspace(6,90, 2500), "Intensity":np.zeros(2500)}).set_index("2theta")
FCC_Ru.loc[np.logical_and(FCC_Ru.index < 40.8, FCC_Ru.index > 40.7), "Intensity"] = 1
FCC_Ru.loc[np.logical_and(FCC_Ru.index < 47.5, FCC_Ru.index > 47.4), "Intensity"] = 0.52
FCC_Ru.loc[np.logical_and(FCC_Ru.index < 69.4, FCC_Ru.index > 69.3), "Intensity"] = 0.35
FCC_Ru.loc[np.logical_and(FCC_Ru.index < 83.75, FCC_Ru.index > 83.6), "Intensity"] = 0.4
FCC_Ru.loc[np.logical_and(FCC_Ru.index < 88.4, FCC_Ru.index > 88.3), "Intensity"] = 0.12
FCC_Ru["Intensity"] += 0.25

HCP_Ru = pd.read_excel("Maria XRD/theoritical_phases.xls", sheet_name = "Ru", index_col= 0 , header = None, usecols=[0,1], names = ["2theta","Intensity"])

FCC_Pt = pd.DataFrame({"2theta":np.linspace(6,90, 2500), "Intensity":np.zeros(2500)}).set_index("2theta")
FCC_Pt.loc[np.logical_and(FCC_Pt.index < 39.8, FCC_Pt.index > 39.7), "Intensity"] = 1
FCC_Pt.loc[np.logical_and(FCC_Pt.index < 46.3, FCC_Pt.index > 46.2), "Intensity"] = 0.53
FCC_Pt.loc[np.logical_and(FCC_Pt.index < 67.5, FCC_Pt.index > 67.4), "Intensity"] = 0.4
FCC_Pt.loc[np.logical_and(FCC_Pt.index < 85.7, FCC_Pt.index > 85.6), "Intensity"] = 0.15
FCC_Pt["Intensity"] += 0.25

FCC_Ir = pd.DataFrame({"2theta":np.linspace(6,90, 2500), "Intensity":np.zeros(2500)}).set_index("2theta")
FCC_Ir.loc[np.logical_and(FCC_Ir.index < 40.7, FCC_Ir.index > 40.6), "Intensity"] = 1
FCC_Ir.loc[np.logical_and(FCC_Ir.index < 47.4, FCC_Ir.index > 47.3), "Intensity"] = 0.52
FCC_Ir.loc[np.logical_and(FCC_Ir.index < 69.2, FCC_Ir.index > 69.1), "Intensity"] = 0.35
FCC_Ir.loc[np.logical_and(FCC_Ir.index < 83.5, FCC_Ir.index > 83.4), "Intensity"] = 0.43
FCC_Ir.loc[np.logical_and(FCC_Ir.index < 88.1, FCC_Ir.index > 88.0), "Intensity"] = 0.12
FCC_Ir["Intensity"] += 0.25

FCC_Pd = pd.DataFrame({"2theta":np.linspace(6,90, 2500), "Intensity":np.zeros(2500)}).set_index("2theta")
FCC_Pd.loc[np.logical_and(FCC_Pd.index < 40.2, FCC_Pd.index > 40.1), "Intensity"] = 1
FCC_Pd.loc[np.logical_and(FCC_Pd.index < 46.7, FCC_Pd.index > 46.6), "Intensity"] = 0.52
FCC_Pd.loc[np.logical_and(FCC_Pd.index < 68.2, FCC_Pd.index > 68.1), "Intensity"] = 0.35
FCC_Pd.loc[np.logical_and(FCC_Pd.index < 82.2, FCC_Pd.index > 82.1), "Intensity"] = 0.4
FCC_Pd.loc[np.logical_and(FCC_Pd.index < 86.7, FCC_Pd.index > 86.6), "Intensity"] = 0.12
FCC_Pd["Intensity"] += 0.25

# Data Retrieval
RuRu = pd.read_csv("Own XRD/maria Ru 2 before co_exported.xy",sep =" ", names=["Intensity"], index_col=0, skiprows= 1 )
RuRh = pd.read_csv("Maria XRD/RuRh 15-09-22.csv", names=["Intensity"], header=None)
RuPt = pd.read_csv("Maria XRD/RuPt 80-20.csv", names=["Intensity"], header=None)
RuPd = pd.read_csv("Maria XRD/RuPd.csv", names=["Intensity"], header=None)
RuIr = pd.read_csv("Maria XRD/RuIr.csv", names=["Intensity"], header=None)
RuPdPt = pd.read_csv("Maria XRD/RuPtPd 15-09-22.csv", names=["Intensity"], header=None)
RuPdPtIr = pd.read_csv("Maria XRD/RuPtPdIr 15-09-22.csv", names=["t", "Intensity"], header=None)
RuPdPtIrRh = pd.read_csv("Maria XRD/sample 21_exported.xy", names=["t", "Intensity"], header=2, sep = " ", index_col=0)

start_value = 6
end_value = HCP_Ru.index[0]  # The current starting value of 2theta in HCP_Ru
# Create a new dataframe with the required range of 2theta values
new_index = pd.Index(np.arange(start_value, end_value, 0.02), name='2theta')
new_data = pd.DataFrame(0, index=new_index, columns=['Intensity'])
HCP_Ru = pd.concat([new_data, HCP_Ru])
HCP_Ru = HCP_Ru[HCP_Ru.index <=90]

# Set style and figure size
plt.style.use("default")
fig = plt.figure(figsize=(7, 10))
ax = fig.add_subplot(111)
plt.xticks(np.arange(0, 80, 10))
# plot the HCP Ru Phase
plt.plot(HCP_Ru["Intensity"], color= "#e63b48")
plt.text(77, 0.45, "HCP Ru", va='baseline', color = "#e63b48")

# Plotting each sample with its color and adding text labels
for i, sample in enumerate(samples): 
    if isinstance(sample, pd.DataFrame):
        plt.plot(sample["Intensity"] + i + 1, color=sample_colors[i])
        plt.text(77, i + intensity_at_81[i] + 1, sample_names[i], va='baseline', color = sample_colors[i])

# Plot FCC phase of Ru, Pt, Pd, Ir, Rh
plt.plot(FCC_Ru["Intensity"] + num_samples + 1, color = "#e63b48", zorder = 3)
plt.text(77, num_samples + 1 + 0.7, "FCC Ru", va='baseline', color = "#e63b48" )

plt.plot(FCC_Pt["Intensity"] + num_samples + 1, color = "#ffd880")
plt.text(77, num_samples + 1 + 0.9 , "FCC Pt", va='baseline', color = "#ffd880")

plt.plot(FCC_Pd["Intensity"] + num_samples + 1, color = "#39deb3")
plt.text(77, num_samples + 1 + 1.1, "FCC Pd", va='baseline', color = "#39deb3")

plt.plot(FCC_Rh["Intensity"] + num_samples + 1, color = "#1a8eb5")
plt.text(77, num_samples + 1 + 1.3, "FCC Rh", va='baseline', color = "#1a8eb5")

plt.plot(FCC_Ir["Intensity"] + num_samples + 1, color = "#204f5e")
plt.text(77, num_samples + 1 + 1.5, "FCC Ir", va='baseline', color = "#204f5e")

# Customizing axis ticks and labels
ax.set_yticks([])
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.xaxis.label.set_fontsize(15)
ax.yaxis.label.set_fontsize(15)

plt.xlabel(r"2$\theta$ (°)")
plt.ylabel("Normalized Intensity")
plt.savefig("data/figures/All_XRDs", dpi=700, bbox_inches='tight')
