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
