# =======================================================
#  - COMBINED SPECTRA PLOT
# =======================================================

from matplotlib.patches import FancyArrowPatch
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


