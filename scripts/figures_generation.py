"""
File: figure_generation.py
Author: Oscar Laurent
Date: 2024-09-11
Description: A Python script to generate the figures on the peaks decomposition process used in the article and in the supplementary information
"""

from utilities.molecule import Molecule, Element
from utilities.periodic_table import get_periodic_table
from utilities.SIMS_Spectra_class import SIMS_Spectra
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, ScalarFormatter)



# =======================================================
#  - COMBINED SPECTRA PLOT
# =======================================================

isotope_df = pd.DataFrame({"bin_mz": np.arange(95, 110.5, 0.5)})

df_list = [
    pd.DataFrame({"mz" : Element("Pd", count = 1).isotopic_weight, "Pd" : Element("Pd", count = 1).isotopic_ratios}),
    pd.DataFrame({"mz" : Element("Ru", count = 1).isotopic_weight, "Ru" : Element("Ru", count = 1).isotopic_ratios}),
    pd.DataFrame({"mz" : Element("Rh", count = 1).isotopic_weight, "Rh" : Element("Rh", count = 1).isotopic_ratios}),
           ]

for df in df_list:
    df.iloc[:, 1] = df.iloc[:, 1]/df.iloc[:, 1].max()
    # Binning the isotopic distribution
    df["bin_mz"] = np.round(df.mz * 0.5**(-1)) * 0.5
    bin_distribution_df = df.groupby("bin_mz").max().reset_index().drop("mz", axis=1)
    isotope_df = pd.merge(isotope_df, bin_distribution_df, how="left", on = "bin_mz")

isotope_df = isotope_df.fillna(0)

# Create a figure with a custom grid layout
fig = plt.figure(figsize=(20, 10))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 2], wspace=0.4)

# Create the first three subplots stacked vertically in the first column
gs_left = gs[0, 0].subgridspec(3, 1, hspace=0.05)
ax0 = fig.add_subplot(gs_left[0])
ax1 = fig.add_subplot(gs_left[1])
ax2 = fig.add_subplot(gs_left[2])

# Create the fourth subplot in the second column
ax3 = fig.add_subplot(gs[0, 1])

# Plot the data in the first three subplots
ax0.bar(isotope_df.bin_mz, isotope_df.Ru, width=0.4, alpha=1, label="Ru", color = "#247e89")
ax1.bar(isotope_df.bin_mz, isotope_df.Pd, width=0.4, alpha=1, label="Pd", color = "#D3811D")
ax2.bar(isotope_df.bin_mz, isotope_df.Rh, width=0.4, alpha=1, label="Rh", color = "#DE1631")

# Add multicolored text to the top left corner of the right subplot
ax3.text(0.15 , 0.95, 'Ru', transform=ax3.transAxes, fontsize=22, verticalalignment='top', fontweight = "bold", color='#247e89')
ax3.text(0.2 , 0.95, '+',  transform=ax3.transAxes, fontsize=22, verticalalignment='top', fontweight = "bold", color='#14213d')
ax3.text(0.23, 0.95, 'Pd', transform=ax3.transAxes, fontsize=22, verticalalignment='top', fontweight = "bold", color='#D3811D')
ax3.text(0.28, 0.95, '+',  transform=ax3.transAxes, fontsize=22, verticalalignment='top', fontweight = "bold", color='#14213d')
ax3.text(0.31, 0.95, 'Rh', transform=ax3.transAxes, fontsize=22, verticalalignment='top', fontweight = "bold", color='#DE1631')

# Add text to the top left corner of each subplot in the left column
ax0.text(0.01, 0.95, 'Ru', transform=ax0.transAxes, fontsize=22, verticalalignment='top', fontweight = "bold", color = "#247e89")
ax1.text(0.01, 0.95, 'Pd', transform=ax1.transAxes, fontsize=22, verticalalignment='top', fontweight = "bold", color = "#D3811D")
ax2.text(0.01, 0.95, 'Rh', transform=ax2.transAxes, fontsize=22, verticalalignment='top', fontweight = "bold", color = "#DE1631")

# Remove y-axis labels for the first three subplots
ax0.get_yaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)
ax3.get_yaxis().set_visible(False)

# Remove x-axis labels for the first two subplots
ax0.get_xaxis().set_visible(False)
ax1.get_xaxis().set_visible(False)

# Customize x-axis for the third subplot
ax2.xaxis.set_major_locator(MultipleLocator(2))
#ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
ax2.tick_params(axis='x', which='major', labelsize=15, length=6, width=2)  # Major ticks
ax2.tick_params(axis='x', which='minor', labelsize=10, length=5, width=1)  # Minor ticks

ax3.xaxis.set_major_locator(MultipleLocator(2))
#ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
ax3.tick_params(axis='x', which='major', labelsize=15, length=6, width=2)  # Major ticks
ax3.tick_params(axis='x', which='minor', labelsize=10, length=5, width=1)  # Minor ticks

# Plot something in the fourth subplot (dummy data used for example)
ax3.bar(isotope_df.bin_mz, isotope_df.Ru, width=0.4, alpha=1, label='Ru', color='#247e89')
ax3.bar(isotope_df.bin_mz, isotope_df.Pd, width=0.4, alpha=1, label='Pd', color='#D3811D', bottom=isotope_df.Ru)
ax3.bar(isotope_df.bin_mz, isotope_df.Rh, width=0.4, alpha=1, label='Rh', color='#DE1631', bottom=isotope_df.Ru + isotope_df.Pd)

# Remove the spines (box) around the right subplot
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['left'].set_visible(False)

ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

for ax in [ax0, ax1, ax2, ax3]:
    for spine in ax.spines.values():
        spine.set_linewidth(2)  # Set the thickness
        spine.set_edgecolor('#14213d')  # Set the color
       
    ax.tick_params(axis='both', which='both', width=2, color='#6D7471')  # Set tick parameters
    ax.tick_params(axis='x', colors='#4F5452')  # Set x-tick labels to grey

# Add arrow annotation
arrow = FancyArrowPatch((0.38, 0.5), (0.48, 0.5), transform=fig.transFigure,
                         mutation_scale=100, lw=2, color='#14213d')
fig.patches.append(arrow)

plt.subplots_adjust(hspace=0.0)
plt.savefig("data/figures/combined_spectra.png", dpi=700, bbox_inches='tight')



# =======================================================
#  - QUANTIFICATION PLOT
# =======================================================

single_element = ["H", "K", "Na", "Cu", "Zn", "Al", "CH3", "C2H5"]
other_elements = []
spectra = SIMS_Spectra("spectrum txt cleaned/P_HEA 1.2_Bi5_1.txt",
                       single_other_component = single_element,
                       indice_max=3,
                       max_molecules=4,
                       other_elements= other_elements)

xlim = (95, 111)
index_range = np.where(np.logical_and(spectra.bin_mz >= xlim[0], spectra.bin_mz <= xlim[1]))
bin_mz = spectra.bin_mz[index_range]
bin_spectra = spectra.bin_intensities[index_range]
computed_spectra = spectra.computed_intensities[index_range]

Ru = spectra.individual_molecule_spectra_df["Ru"].values[index_range]
Rh = spectra.individual_molecule_spectra_df["Rh"].values[index_range]
Pd = spectra.individual_molecule_spectra_df["Pd"].values[index_range]

# Create subplots
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(10, 10))

# First subplot: Compare observed spectra and computed spectra
ax1.plot(bin_mz, bin_spectra,linewidth = 2, color ="#494B70", label=r"Observed Spectrum  ($\tilde{\mathbf{Y}}_{obs}$)")
ax1.plot(bin_mz, computed_spectra, color='red', label=r"Computed Spectrum ($\tilde{\mathbf{Y}}_{calc}$)", linewidth=2, linestyle ="--", dashes=(2, 2))
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.tick_params(axis='y', labelsize=14)
ax1.yaxis.get_offset_text().set_fontsize(14)
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.get_xaxis().set_visible(False)
ax1.set_ylabel("Counts", fontsize=16)
ax1.legend(fontsize=15, frameon=False)

# Add annotation "A" to the first subplot
ax1.text(0.025, 0.94, 'A', transform=ax1.transAxes, fontsize=18, va='center', ha='center',
         bbox=dict(facecolor='#DBDBDB', edgecolor='grey', boxstyle='square,pad=0.3'))

# Second subplot: Decomposition of observed spectra into individual elements
ax2.plot(bin_mz, bin_spectra, linewidth=2,  alpha=1, color = "#494B70", label=r"Observed Spectrum ($\tilde{\mathbf{Y}}_{obs}$)")
ax2.bar(bin_mz, Ru, width=0.4, alpha=1, color='#247e89')
ax2.bar(bin_mz, Rh, width=0.4, alpha=1, color='#DE1631', bottom=Ru)
ax2.bar(bin_mz, Pd, width=0.4, alpha=1, color='#D3811D', bottom=Ru+Rh)
ax2.xaxis.set_major_locator(MultipleLocator(2))
ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
ax2.tick_params(axis='x', which='major', labelsize=12)
ax2.legend(fontsize=15, frameon=False)
ax2.tick_params(axis='y', labelsize=14)
ax2.yaxis.get_offset_text().set_fontsize(14)
ax2.set_ylabel("Counts", fontsize=16)
ax2.set_xlabel("m/z", fontsize=16)

# Add annotations with arrows for individual elements in the lower plot
Ru_quantification = spectra.quantification_values[np.where(np.array(spectra.isotopes_str) =="Ru")].item()
Rh_quantification = spectra.quantification_values[np.where(np.array(spectra.isotopes_str) =="Rh")].item()
Pd_quantification = spectra.quantification_values[np.where(np.array(spectra.isotopes_str) =="Pd")].item()

ax2.annotate(r"$\mathbf{\tilde{x}_{Ru}} = $" + f"{Ru_quantification:.0f}" , xy=(bin_mz[14] - 0.05, Ru[14] +1000), xytext=(bin_mz[14] - 3, Ru[14] + 8000),
             arrowprops=dict(facecolor='#247e89', edgecolor='#247e89', shrink=0.05, width = 3), fontsize=15, fontweight = "bold",  color='#247e89')

ax2.annotate(r"$\mathbf{\tilde{x}_{Rh}} = $" +  f"{Rh_quantification:.0f}", xy=(bin_mz[16] + 0.09, Ru[16] + Rh[16] ), xytext=(bin_mz[16] +3 , Ru[16] + Rh[16] - 10000),
             arrowprops=dict(facecolor='#DE1631', edgecolor='#DE1631', shrink=0.05, width = 3), fontsize=15, fontweight = "bold", color='#DE1631')

ax2.annotate(r"$\mathbf{\tilde{x}_{Pd}} = $" + f"{Pd_quantification:.0f}", xy=(bin_mz[22], Pd[22]), xytext=(bin_mz[22] + 2, Pd[22] + 8000),
             arrowprops=dict(facecolor='#D3811D', edgecolor='#D3811D', shrink=0.05, width = 3), fontsize=15, fontweight = "bold", color='#D3811D')

# Add annotation "B" to the second subplot
ax2.text(0.025, 0.94, 'B', transform=ax2.transAxes, fontsize=18, va='center', ha='center',
         bbox=dict(facecolor='#DBDBDB', edgecolor='grey', boxstyle='square,pad=0.3'))

# Adjust layout
plt.ylim(0, 40000)
plt.xlim(95, 111)
plt.tight_layout()
plt.savefig("quantifications.png", dpi=700, bbox_inches='tight')



# =======================================================
#  - QUANTIFICATION SUBPLOTS
# =======================================================

single_element = ["H", "K", "Na", "Cu", "Zn", "Al", "CH3", "C2H5"]
other_elements = []
spectra = SIMS_Spectra("../data/ToF-SIMS_data/raw_spectra/P_HEA19_Bi5_1.txt",
                       single_other_component = single_element,
                       indice_max=3,
                       #components= ["Ru", "Pt", "Pd"],
                       max_molecules=4,
                       other_elements= other_elements)

xlim = (650, 1540)
index_range = np.where(np.logical_and(spectra.bin_mz >= xlim[0], spectra.bin_mz <= xlim[1]))
bin_mz = spectra.bin_mz[index_range]
bin_spectra = spectra.bin_intensities[index_range]
computed_spectra = spectra.computed_intensities[index_range]
df = spectra.individual_molecule_spectra_df.iloc[index_range]
df = df.iloc[:, (df !=0).any().values]

# Create subplots
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=True, figsize=(10, 10))

# First subplot: Compare observed spectra and computed spectra
ax1.plot(bin_mz, bin_spectra,linewidth = 0.5, color ="#494B70", label=r"Observed Spectrum ($\tilde{\mathbf{Y}}_{obs}$)")
ax1.plot(bin_mz, computed_spectra, color='red', label=r"Computed Spectrum ($\tilde{\mathbf{Y}}_{calc}$)", alpha = 0.5, linewidth = 0.5)
ax1.set_xlim(xlim)
ax1.xaxis.set_major_locator(MultipleLocator(100))
ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
ax1.tick_params(axis='x', which='major', labelsize=12)
ax1.tick_params(axis='y', labelsize=14)
ax1.yaxis.get_offset_text().set_fontsize(14)
ax1.get_xaxis().set_visible(True)
ax1.set_ylabel("Counts", fontsize=16)
ax1.legend(fontsize=15, frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.15),  ncol =5)
inset_xlim = (950, 1050)

# Filter the data for the inset plot
inset_index_range = np.where(np.logical_and(bin_mz >= inset_xlim[0], bin_mz <= inset_xlim[1]))
inset_bin_mz = bin_mz[inset_index_range]
inset_bin_spectra = bin_spectra[inset_index_range]
inset_computed_spectra = computed_spectra[inset_index_range]

# Create the inset plot for the first subplot
axins1 = ax1.inset_axes([0.5, 0.5, 0.47, 0.47], xlim = inset_xlim, ylim = (0.75, 75))
axins1.plot(inset_bin_mz, inset_bin_spectra, linewidth=0.5, color="#494B70", label=r"Observed Spectrum ($\tilde{\mathbf{Y}}_{obs}$)")
axins1.plot(inset_bin_mz, inset_computed_spectra, color='red', label=r"Computed Spectrum ($\tilde{\mathbf{Y}}_{calc}$)", alpha = 0.75, linewidth = 0.75)
axins1.yaxis.set_visible(False)
axins1.xaxis.set_visible(False)
ax1.indicate_inset_zoom(axins1, edgecolor="black")

# Add annotation "A" to the first subplot
ax1.text(0.03, 0.94, 'A', transform=ax1.transAxes, fontsize=18, va='center', ha='center',
         bbox=dict(facecolor='#DBDBDB', edgecolor='grey', boxstyle='square,pad=0.3'))

# Second subplot: Decomposition of observed spectra into individual elements
ax2.plot(bin_mz, bin_spectra, linewidth=0.5,  alpha=1, color = "#494B70", label=r"Observed Spectrum ($\tilde{\mathbf{Y}}_{obs}$)")
for column in df.columns[1:]:
    ax2.plot(bin_mz, df[column], label=column, alpha = 0.75)

ax2.xaxis.set_major_locator(MultipleLocator(100))
ax2.xaxis.set_minor_locator(AutoMinorLocator(10))
ax2.tick_params(axis='x', which='major', labelsize=12)
#ax2.legend(fontsize=15, frameon=False)
ax2.tick_params(axis='y', labelsize=14)
ax2.yaxis.get_offset_text().set_fontsize(14)
ax2.set_ylabel("Counts", fontsize=16)
ax2.set_xlabel("m/z", fontsize=16)
inset_xlim = (950, 1050)

# Add annotation "B" to the second subplot
ax2.text(0.03, 0.94, 'B', transform=ax2.transAxes, fontsize=18, va='center', ha='center',
         bbox=dict(facecolor='#DBDBDB', edgecolor='grey', boxstyle='square,pad=0.3'))

# Filter the data for the inset plot
inset_index_range = np.where(np.logical_and(bin_mz >= inset_xlim[0], bin_mz <= inset_xlim[1]))
inset_bin_mz = bin_mz[inset_index_range]
inset_bin_spectra = bin_spectra[inset_index_range]
inset_df = df.iloc[inset_index_range]
bottom_values = np.zeros_like(inset_bin_mz)

# Create the inset plot
axins = ax2.inset_axes([0.5, 0.5, 0.47, 0.47], xlim = inset_xlim, ylim = (0.75, 75), yticklabels=[], yticks = [])
axins.yaxis.set_visible(False)
axins.xaxis.set_visible(False)
axins.plot(inset_bin_mz, inset_bin_spectra, linewidth=0.5, alpha=1, color="#494B70", label=r"Observed Spectrum ($\mathbf{Y}$)")
for column in inset_df.columns[1:]:
    axins.plot(inset_bin_mz, inset_df[column], label=column, alpha=0.75)
    #axins.bar(inset_bin_mz, inset_df[column], label=column, alpha=0.75, bottom=bottom_values)
    #bottom_values += inset_df[column]

# Adjust layout
ax2.indicate_inset_zoom(axins, edgecolor="black")
plt.xlim(xlim[0], xlim[1])
plt.tight_layout()
plt.savefig("data/figures/quantifications_subplots.png", dpi=700, bbox_inches='tight')



# =======================================================
#  SUPPLEMENTARY MATERIAL - RAW SPECTRUM PLOT 
# =======================================================

spectra = pd.read_csv("data/ToF-SIMS_data/raw_spectra/P_HEA 1.2_Bi5_3.txt", sep='\t', skiprows=3, usecols=[1, 2], names=['mz', 'intensity'])

figure, axis = plt.subplots(2, 1, figsize=(11, 6)) 

axis[0].plot(spectra["mz"], spectra["intensity"], linewidth = 0.5)
axis[0].ticklabel_format(axis = "y", style = "sci")
axis[0].set_yscale("linear")
axis[0].set_xticklabels([])
axis[0].get_xaxis().set_visible(False)
axis[0].tick_params(axis='y', labelsize=12)

axis[1].plot(spectra["mz"], spectra["intensity"], linewidth = 0.5)
axis[1].set_yscale("log")

plt.xlabel("m/z", size = 14)
plt.xticks(size = 12)
plt.yticks(size = 12)

plt.subplots_adjust(hspace=0.0)
plt.savefig("data/supplementary_material/raw_spectrum", dpi=700, bbox_inches='tight')



# =======================================================
#  SUPPLEMENTARY MATERIAL - BINNED SPECTRA PLOT 
# =======================================================

raw = pd.DataFrame({"mz": spectra.mz, "intensity": spectra.intensities})
binned = pd.DataFrame({"mz": spectra.bin_mz, "intensity": spectra.bin_intensities})

fig, ax = plt.figure(figsize = (11,6)) , plt.gca()

plt.plot(raw.mz, raw.intensity, label = "Raw Spectrum", color = "royalblue", linewidth = 2)
plt.bar(binned.mz, binned.intensity, label = "Binned Spectrum", alpha = 0.6, width=0.45, color = "#fb8500")
plt.xlim(199,221)
plt.ylim(0,1600)
plt.legend(prop={'size': 20})
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))

ax.tick_params(axis='both', which='major', labelsize=14, length=7, width=2)  # Major ticks
ax.tick_params(axis='both', which='minor', labelsize=14, length=6, width=1)  # Minor ticks

# Adding labels and customizing their size
ax.set_xlabel("m/z", fontsize=16)
ax.set_ylabel("Counts", fontsize=16)

plt.savefig("data/supplementary_material/binned_spectrum", dpi=500, bbox_inches='tight')


# =======================================================
#  - REVIEWERS COMMENTS : INTERMEDIATE M/Z RANGE PLOT
# =======================================================

single_element = ["H", "K", "Na", "Cu", "Zn", "Al", "CH3", "C2H5"]
other_elements = []
spectra = SIMS_Spectra("../data/ToF-SIMS_data/raw_spectra/P_HEA19_Bi5_1.txt",
                       single_other_component = single_element,
                       indice_max=3,
                       #components= ["Ru", "Pt", "Pd"],
                       max_molecules=4,
                       other_elements= other_elements)


xlim = (150, 640)
index_range = np.where(np.logical_and(spectra.bin_mz >= xlim[0], spectra.bin_mz <= xlim[1]))
bin_mz = spectra.bin_mz[index_range]
bin_spectra = spectra.bin_intensities[index_range]
computed_spectra = spectra.computed_intensities[index_range]
df = spectra.individual_molecule_spectra_df.iloc[index_range]
df = df.iloc[:, (df !=0).any().values]

# Create subplots
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=True, figsize=(10, 10))

# First subplot: Compare observed spectra and computed spectra
ax1.plot(bin_mz, bin_spectra,linewidth = 0.5, color ="#494B70", label=r"Observed Spectrum ($\tilde{\mathbf{Y}}_{obs}$)")
ax1.plot(bin_mz, computed_spectra, color='red', label=r"Computed Spectrum ($\tilde{\mathbf{Y}}_{calc}$)", alpha = 0.5, linewidth = 0.5)
ax1.set_xlim(xlim)
ax1.xaxis.set_major_locator(MultipleLocator(100))
ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
ax1.tick_params(axis='x', which='major', labelsize=12)
ax1.tick_params(axis='y', labelsize=14)
ax1.yaxis.get_offset_text().set_fontsize(14)
ax1.get_xaxis().set_visible(True)
ax1.set_ylabel("Counts", fontsize=16)
ax1.legend(fontsize=15, frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.15),  ncol =5)

# Filter the data for the inset plot
inset_xlim = (480, 540)
inset_index_range = np.where(np.logical_and(bin_mz >= inset_xlim[0], bin_mz <= inset_xlim[1]))
inset_bin_mz = bin_mz[inset_index_range]
inset_bin_spectra = bin_spectra[inset_index_range]
inset_computed_spectra = computed_spectra[inset_index_range]

# Create the inset plot for the first subplot
axins1 = ax1.inset_axes([0.5, 0.5, 0.47, 0.47], xlim = inset_xlim, ylim = (0.75, 340))
axins1.plot(inset_bin_mz, inset_bin_spectra, linewidth=0.5, color="#494B70", label=r"Observed Spectrum ($\tilde{\mathbf{Y}}_{obs}$)")
axins1.plot(inset_bin_mz, inset_computed_spectra, color='red', label=r"Computed Spectrum ($\tilde{\mathbf{Y}}_{calc}$)", alpha = 0.75, linewidth = 0.75)
axins1.yaxis.set_visible(False)
axins1.xaxis.set_visible(False)
ax1.indicate_inset_zoom(axins1, edgecolor="black")

# Add annotation "A" to the first subplot
ax1.text(0.03, 0.94, 'A', transform=ax1.transAxes, fontsize=18, va='center', ha='center',
         bbox=dict(facecolor='#DBDBDB', edgecolor='grey', boxstyle='square,pad=0.3'))

# Second subplot: Decomposition of observed spectra into individual elements
ax2.plot(bin_mz, bin_spectra, linewidth=0.5,  alpha=1, color = "#494B70", label=r"Observed Spectrum ($\tilde{\mathbf{Y}}_{obs}$)")
for column in df.columns[1:]:
    ax2.plot(bin_mz, df[column], label=column, alpha = 0.75)

ax2.xaxis.set_major_locator(MultipleLocator(100))
ax2.xaxis.set_minor_locator(AutoMinorLocator(10))
ax2.tick_params(axis='x', which='major', labelsize=12)
#ax2.legend(fontsize=15, frameon=False)
ax2.tick_params(axis='y', labelsize=14)
ax2.yaxis.get_offset_text().set_fontsize(14)
ax2.set_ylabel("Counts", fontsize=16)
ax2.set_xlabel("m/z", fontsize=16)


# Add annotation "B" to the second subplot
ax2.text(0.03, 0.94, 'B', transform=ax2.transAxes, fontsize=18, va='center', ha='center',
         bbox=dict(facecolor='#DBDBDB', edgecolor='grey', boxstyle='square,pad=0.3'))

# Filter the data for the inset plot
inset_index_range = np.where(np.logical_and(bin_mz >= inset_xlim[0], bin_mz <= inset_xlim[1]))
inset_bin_mz = bin_mz[inset_index_range]
inset_bin_spectra = bin_spectra[inset_index_range]
inset_df = df.iloc[inset_index_range]
bottom_values = np.zeros_like(inset_bin_mz)

# Create the inset plot
axins = ax2.inset_axes([0.5, 0.5, 0.47, 0.47], xlim = inset_xlim, ylim = (0.75, 340), yticklabels=[], yticks = [])
axins.yaxis.set_visible(False)
axins.xaxis.set_visible(False)
axins.plot(inset_bin_mz, inset_bin_spectra, linewidth=0.5, alpha=1, color="#494B70", label=r"Observed Spectrum ($\mathbf{Y}$)")
for column in inset_df.columns[1:]:
    axins.plot(inset_bin_mz, inset_df[column], label=column, alpha=0.75)
    #axins.bar(inset_bin_mz, inset_df[column], label=column, alpha=0.75, bottom=bottom_values)
    #bottom_values += inset_df[column]

# Adjust layout
ax2.indicate_inset_zoom(axins, edgecolor="black")
plt.xlim(xlim[0], xlim[1])
plt.tight_layout()