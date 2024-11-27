import pandas as pd 
import numpy as np 
from dataclasses import dataclass, field
from .molecule import Molecule, Element
from .periodic_table import get_periodic_table
import re
import os 
import scipy.signal as sp
from scipy import stats 
from scipy.ndimage import gaussian_filter1d
from scipy import optimize
import itertools
import plotly.graph_objects as go
from tqdm import tqdm
import re
from scipy.optimize import nnls
import sys 

@dataclass
class SIMS_Spectra: 
    # MetaData
    sample_name: str = field(init=False)
    state: str = field(init=False)
    polarity: str = field(init=False)
    primary_ion: str = field(init=False)
    replicat: int = field(init=False)

    # Hidden Propreties
    filepath: int = field(repr=False)

    mz: np.ndarray = field(init=False, repr=False)
    intensities: np.ndarray = field(init=False, repr=False)
    bin_intensities: float = field(init=False, repr=False)
    bin_mz: np.ndarray = field(init=False, repr=False)
    bin_intensities: np.ndarray = field(init=False, repr=False)
    computed_intensities: np.ndarray = field(init=False, repr=False)

    components: list[str] = field(repr=False, default_factory= lambda: ["Ru", "Rh", "Pd", "Pt", "Ir"])
    single_other_component: list[str] = field(repr=False, default_factory= lambda: [])
    other_elements: list[str] = field(repr=False, default_factory= lambda: [])

    isotopes: list[Molecule] = field(default_factory=lambda: None,  repr=False)
    isotopes_str: list[str] = field(init=False,  repr=False)
    isotopes_df: pd.DataFrame = field(default_factory=lambda: None,  repr=False)
    individual_molecule_spectra: pd.DataFrame = field(init=False, repr=False)

    quantified_molecules_df: pd.DataFrame = field(init=False, repr=False)
    quantification_values: np.ndarray = field(init=False, repr=False)
    residuals: float = field(init=False, repr=False)

    indice_max: int = field(repr=False, default_factory=lambda: 5)
    max_molecules: int = field(repr=False, default_factory=lambda: 4)
    mode: str = field(repr=False, default_factory=lambda: "nnls")
    highest_peak_mz: float = field(repr=False, default_factory=lambda: None)
    number_of_peak_found: int = field(init=False, repr=False)

    
    def __post_init__(self):
        self.extract_metadata_from_filename()
        self.extract_spectra_from_file()
        self.compute_bin_spectra()
        print(f"\nSpectra loaded from '{self.filepath}'\n")

        # Change the highest mz
        if self.highest_peak_mz is None:
            self.highest_peak_mz = self.bin_mz.max()
        else: 
            self.highest_peak_mz = min(self.bin_mz.max(), self.highest_peak_mz)
            
        if self.isotopes_df is None or self.isotopes is None: 
            print(f"Creation of the Isotope DataFrame with {self.components} up to {self.indice_max} indices for the mz range of {self.mz.min():.0f} to {self.highest_peak_mz:.0f}")
            self.compute_isotope_combinaison(indice_max=self.indice_max,
                                             max_molecules=self.max_molecules,
                                             highest_peak_mass=self.highest_peak_mz,
                                             other_elements=self.other_elements)
            print(f"\nIsotope DataFrame created, {len(self.isotopes)} isotopes considered for quantification")

        # Compute the list of isotope names
        self.isotopes_str = list(map(lambda isotope : isotope.molecular_formula if isinstance(isotope, Molecule) else 
                                isotope.symbol, self.isotopes))
        
        self.isotope_quantification(mode = self.mode)
        s = f"""Quantification done --> {len(self.quantified_molecules_df.columns)-1} molecules out of
        {len(self.isotopes)} where quantified in the spectra"""
        print(re.sub(r"\s{2,}", " ", s))
        print(f"Residual of the fit: {self.residuals:.0f}")
        print("\n----------------------------------------------------------------------------------------------")

    
    def get_metadata(self) -> dict: 
        metadata_dict = {"sample_name": self.sample_name,
                         "state": self.state,
                         "polarity": self.polarity, 
                         "primary_ion": self.primary_ion,
                         "replicat": self.replicat}
        return metadata_dict
    

    def get_output_df(self, norm=False) -> pd.DataFrame:
        """
        Return a pd DataFrame with the metatdata sample_name | state | polarity | primary_ion | replicat | Molecules...
        """
        df_metadata = pd.DataFrame(data=self.get_metadata(), index = [0])
        if norm : 
            df_isotope = pd.DataFrame([self.quantification_values/self.intensities.sum()], columns=self.isotopes_str)
        else : 
            df_isotope = pd.DataFrame([self.quantification_values], columns=self.isotopes_str)
        output_df = pd.concat([df_metadata, df_isotope], axis = 1)
        return output_df
    

    def extract_metadata_from_filename(self) -> None:
        """
        Extract metadata from filepath using regex 
        The filename shoud have the following synthax  : 
        N_HEA 1.2 spent_Bi3_1.txt or P_HEA long shot_Bi5_3.txt 
        """
        
        filename = os.path.basename(self.filepath)
        try:
            self.polarity = re.search("[NP]" ,filename).group(0)
            filename = re.sub("[NP]_", "", filename)
            self.primary_ion = re.search("Bi\\d\\+*", filename).group(0)
            filename = re.sub("Bi\\d\\+*_", "", filename)
            self.replicat = re.search("(\\d).txt", filename).group(1)
            filename = re.sub("_(\\d).txt", "", filename)

            if re.search("spent", filename):
                self.state = "spent"
                filename = re.sub("spent", "", filename)
            else : 
                self.state = "new"

            self.sample_name = filename.strip().replace(" ", "_")

        except Exception as e:
            raise ValueError(f"Failed to read metadata from filename {filename}: {e}")

    
    def extract_spectra_from_file(self) -> None:
        # Read the mz and intensity array from the .txt file 
        try:
            data = pd.read_csv(self.filepath, sep='\t', skiprows=3, usecols=[1, 2], names=['mz', 'intensity'])
            self.mz = data['mz'].values
            self.intensities = data['intensity'].values
        except Exception as e:
            raise ValueError(f"Failed to read spectra data from filename {self.filepath}: {e}")
        

    def compute_isotope_combinaison(self, indice_max, max_molecules, highest_peak_mass, other_elements) -> None:
        """
        metals : a liste of melals ["Ru", "Pt"]
        other_elements : a list of the additional elements to be considerete bounded with the 5 metals ["H", "O"] 
        other_single_elements : a list of the other elements thta stay alone ["Cl", "K"]
        indice_max : the maximum indice to consider 
        max_molecules : the maximum indice to consider 
        highest_peak_mass : the highest mass to consider for the molecule
        """

        if max_molecules < 1 or max_molecules > 4:
            raise ValueError(f"max_molecules = {max_molecules}, please chose a value between 1 and 4")

        combinaison_list = []

        # Removing the element that are in other_elements,ts and other_single_elements
        single_other_component = [component for component in self.single_other_component if component not in other_elements]

        # Addition of the single element 
        if max_molecules >= 1 :

            # Addition of the single other element such as Cu, Ar, CH3, ...
            for other_element_str in single_other_component:
                try:
                    other_element = Element(other_element_str, 1)
                except:
                    other_element = Molecule(other_element_str)
                combinaison_list.append(other_element)
            
            # Creation of a list of single componennt molecules ["Ru", "Ru2",..., "Ir5", "H", "H2",...]
            Ax_list = [re.sub("(?<=\D{2})1$", "", f"{component}{indice}") for component, indice in 
                       itertools.product(self.components + other_elements, range(1,indice_max + 1))]

            # Addition of single element 
            for Ax_str in Ax_list:
                if bool(re.search("\d", Ax_str)) or len(re.findall(r'[A-Z]',Ax_str)) >= 2 : 
                    Ax = Molecule(Ax_str)
                    if Ax.molecular_weight < highest_peak_mass:
                        combinaison_list.append(Ax)
                else:
                    Ax = Element(Ax_str, 1)
                    combinaison_list.append(Ax)
                    

        # Addition of combinaison of 2 elements clusters like RuRh, RuPd, ..., Pt5Ir5 
        # But also the combinaison with other element such as RuO, RuH, ..., IrH
        if max_molecules >= 2 :
            for combi_tuple in itertools.combinations(self.components + other_elements, 2):
                    A, B  = combi_tuple
                    if A in other_elements and B in other_elements :
                        pass
                    else : 
                        list_A = [x for x in Ax_list if A in x]
                        list_B = [x for x in Ax_list if B in x]
                        for combi_tuple in itertools.product(list_A, list_B):
                            A_xB_x = Molecule("".join(combi_tuple))
                            if A_xB_x.molecular_weight < highest_peak_mass:
                                combinaison_list.append(A_xB_x)

        # Addition of combinaison of 3 element clusters 
        if max_molecules >= 3:
            for combi_tuple in itertools.combinations(self.components + other_elements, 3):
                A, B ,C = combi_tuple
                if [x in other_elements for x in combi_tuple].count(True) >= 2:
                    pass
                else : 
                    list_A = [x for x in Ax_list if A in x]
                    list_B = [x for x in Ax_list if B in x]
                    list_C = [x for x in Ax_list if C in x]
                    for combi_tuple in itertools.product(list_A, list_B, list_C):
                        A_xB_xC_x = Molecule("".join(combi_tuple))
                        if A_xB_xC_x.molecular_weight < highest_peak_mass:
                            combinaison_list.append(A_xB_xC_x)
    
        # Addition of combinaison of 4 element clusters 
        if max_molecules >= 4:
            for combi_tuple in itertools.combinations(self.components + other_elements, 4):
                A, B ,C, D = combi_tuple
                if [x in other_elements for x in combi_tuple].count(True) >= 3:
                    pass
                else : 
                    list_A = [x for x in Ax_list if A in x]
                    list_B = [x for x in Ax_list if B in x]
                    list_C = [x for x in Ax_list if C in x]
                    list_D = [x for x in Ax_list if D in x]
                    for combi_tuple in itertools.product(list_A, list_B, list_C, list_D):
                        A_xB_xC_xD_x = Molecule("".join(combi_tuple))
                        if A_xB_xC_xD_x.molecular_weight < highest_peak_mass:
                            combinaison_list.append(A_xB_xC_xD_x)

        self.isotopes = combinaison_list

        # Fetching of the istopic distribution of the Molecules object and grouping them into a DataFrame
        isotope_df = pd.DataFrame(self.bin_mz[self.bin_mz <= self.highest_peak_mz], columns=["bin_mz"])
        for molecule in tqdm(self.isotopes): 
            if isinstance(molecule, Element): 
                distribution_tuple = (np.array(molecule.isotopic_weight), np.array(molecule.isotopic_ratios))
                molecule_str = molecule.symbol
            else : 
                distribution_tuple = molecule.isotopic_distribution()
                molecule_str = molecule.molecular_formula

            distribution_df = pd.DataFrame({"mz" : distribution_tuple[0], molecule_str: distribution_tuple[1]})

            # Scaling of the intensities values between 0 and 1 
            distribution_df[molecule_str] = distribution_df[molecule_str]/distribution_df[molecule_str].max()

            # Binning of the isotopic distributions 
            distribution_df["bin_mz"] = np.round(distribution_df.mz * self.bin_size**(-1)) * self.bin_size
            bin_distribution_df = distribution_df.groupby("bin_mz").max().reset_index().drop("mz", axis=1)

            # Mapping the proportion to the hole spectra 
            isotope_df = pd.merge(isotope_df, bin_distribution_df, how="left", on = "bin_mz")
            isotope_df.fillna(0, inplace=True)

        self.isotopes_df = isotope_df
        

    def compute_peak_detection(self, dis=1, h=1, p=1, display=False, on = "raw") -> go.Figure:
        """
        Perform the detection of the peaks.

        Parameters:
        on (str): either "raw" or "binned", the intensity array to perform the detection

        Returns:
        np.ndarray: The smoothed intensities.
        """

        if on == "raw":
            peaks, _ = sp.find_peaks(self.intensities, distance=dis, height=h, prominence=p)
            max_mz = self.mz[peaks[-1]]
            self.highest_peak_mz = round(max_mz, 2)
            self.number_of_peak_found = len(peaks)
            print(self.number_of_peak_found, "Peaks found")
            print("Highest Peak mz:", self.highest_peak_mz)
            if display:
                fig = go.Figure()

                fig.add_trace(go.Scatter(x=self.mz,
                                        y=(self.intensities),
                                        name='Raw Spectra',
                                        opacity = 1,
                                        line=dict(color='royalblue', width=1)))

                fig.add_trace(go.Scatter(x=self.mz[peaks],
                                         y=self.intensities[peaks],
                                         name='Peaks Found',
                                         mode="markers",
                                         marker=dict(color='firebrick', size = 5, symbol="x")))
                
                fig.update_layout(title= f"{self.sample_name}, {self.primary_ion}, {self.replicat}",
                                  height = 600,
                                  template = "none",
                                  margin=dict(t=50, b=50, l = 50, r = 50))

                fig.show()

            
        elif on == "binned":
            peaks, _ = sp.find_peaks(self.bin_intensities, distance=dis, height=h, prominence=p)
            max_mz = self.mz[peaks[-1]]
            self.highest_peak_mz = round(max_mz, 2)
            self.number_of_peak_found = len(peaks)
            print(self.number_of_peak_found, "Peaks found")
            print("Highest Peak mz:", self.highest_peak_mz)
            if display:
                fig = go.Figure()

                fig.add_trace(go.Bar(x=self.bin_mz,
                                     y=self.bin_intensities,
                                     name=' Binned Spectra',
                                     opacity = 1))

                fig.add_trace(go.Scatter(x=self.bin_mz[peaks],
                                         y=self.bin_intensities[peaks],
                                         name='Peaks Found',
                                         mode="markers",
                                         marker=dict(color='firebrick', size = 5, symbol="x")))
                
                fig.update_layout(title= f"{self.sample_name}, {self.primary_ion}, {self.replicat}",
                                  height = 600,
                                  template = "none",
                                  margin=dict(t=50, b=50, l = 50, r = 50))

                fig.show()

        else : 
            raise ValueError(f"'{on}' is not a valide argument chose between 'raw' and 'binned'")
      

    def compute_gaussian_approx(self, sigma: int = 2) -> go.Figure:
        """
        Smooth the spectra using a Gaussian filter.

        Parameters:
        sigma (int): The standard deviation of the Gaussian filter.

        Returns:
        np.ndarray: The smoothed intensities.
        """
        smoothed_spectra = gaussian_filter1d((self.intensities), sigma=sigma)

        fig = go.Figure()

        fig.add_trace(go.Scatter(x=self.mz, y=(self.intensities), name='Raw Spectra',
                         line=dict(color='royalblue', width=1)))
        fig.add_trace(go.Scatter(x=self.mz, y=smoothed_spectra, name='Smoothed',
                         line=dict(color='firebrick', width=1)))
        
        fig.update_layout(title= f"{self.sample_name}, {self.primary_ion}, {self.replicat}",
                      height = 600,
                      template = "none",
                      margin=dict(t=50, b=50, l = 50, r = 50))

        fig.show()


    def compute_bin_spectra(self, bin_size: float = 0.5, plot: bool = False) -> None :
        """
        Compute the binned spectra with each bin being bin_size m/z.

        Parameters:
        bin_size (float): The size of each bin in m/z.

        Returns:
        pd.DataFrame: A DataFrame with binned m/z and corresponding intensities.
        """

        raw_df = pd.DataFrame({'mz': self.mz, 'intensity': self.intensities})
        raw_df["bin_mz"] = np.round(raw_df.mz * bin_size**(-1)) * bin_size
        bin_df = raw_df.groupby("bin_mz").max().reset_index().drop("mz", axis=1)  # Get the max of the intensity in each bin

        if plot : 
            fig = go.Figure()

            fig.add_trace(go.Scatter(
                x=self.mz, 
                y=(self.intensities), 
                name='Raw Spectra', 
                opacity= 0.6,
                line=dict(color='royalblue', width=1)))
            
            fig.add_trace(go.Bar(
                x=bin_df.bin_mz,
                y=bin_df.intensity,
                name='Binned Spectra',
                opacity=0.6))

            fig.update_layout(title= f"{self.sample_name}, {self.primary_ion}, {self.replicat}",
                          height = 600,
                          template = "none",
                          margin=dict(t=50, b=50, l = 50, r = 50))
            fig.show()
        
        self.bin_mz = bin_df.bin_mz.values
        self.bin_intensities = bin_df.intensity.values
        self.bin_size = bin_size
        

    def isotope_quantification(self, quantification_threshold: float = 0.1, mode = "nnls") -> None:
        if not isinstance(self.isotopes_df, pd.DataFrame) : 
            raise ValueError(f"'{self.isotopes_df}' is not a DataFrame. Please pass a valid isotope DataFrame")
        
        if mode != "nnls" and mode != "nnls_admm":
            raise ValueError(f"'{mode}' is not a valid mode, please use 'nnls' or 'nnls_admm'")

        nrow_isotope_df = len(self.isotopes_df)
        nrow_bin_mz = len(self.bin_mz)
        if nrow_isotope_df != nrow_bin_mz:
            print(f"isotope df of shape {self.isotopes_df.shape} does not match with bin_mz array of shape {self.bin_mz.shape}...")
            if nrow_bin_mz > nrow_isotope_df :
                A = self.isotopes_df.drop("bin_mz", axis=1).to_numpy() 
                Y = self.bin_intensities[:nrow_isotope_df]
                print(f"bin_mz array reshaped to {Y.shape} ")
            else : 
                A = self.isotopes_df[:min(nrow_bin_mz, nrow_isotope_df)].drop("bin_mz", axis=1).to_numpy()
                Y = self.bin_intensities
                print(f"Isotope df reshaped to {A.shape} ")
        else :
            A = self.isotopes_df.drop("bin_mz", axis=1).to_numpy()
            Y = self.bin_intensities
        
        if mode == "nnls":
            X, residuals = nnls(A, Y)
        elif mode == "nnls_admm":
            X, residuals = nnls_admm(A, Y, Y, tol= 1e-10, max_iter = 2000)
        
        self.quantification_values = X
        self.computed_intensities = np.dot(A, X)
        self.residuals = residuals

        # Taking only molecules quantified above the threhold value
        molecule_indices = np.where(X > quantification_threshold)[0]
        indices = np.insert(molecule_indices+1, 0, 0) # addition the the indice of the bin_mz col
        quantified_molecules_df = self.isotopes_df.iloc[:, indices]
        self.quantified_molecules_df = quantified_molecules_df  

        # Compute the spectra for each molecules separetly 
        individual_molecule_spectra = self.quantified_molecules_df.drop("bin_mz", axis =1)*X[X>quantification_threshold]
        self.individual_molecule_spectra_df = individual_molecule_spectra
        
        return


    def plot_results(self, yaxis_type = "linear", plot_computed_pectra = False) -> go.Figure:
        # Plotting the Compounds into the actual spectra
        fig = go.Figure()

        fig.add_trace(go.Scatter(x=self.bin_mz, y=self.bin_intensities, name='True spectra', 
                                 line=dict(color='royalblue', width=1, dash = "solid")))
        
        if plot_computed_pectra:
            fig.add_trace(go.Scatter(x=self.bin_mz, y=self.computed_intensities, name='Computed spectra', 
                                     line=dict(color='firebrick', width=1, dash = "dot")))

        # Plot the different compound 
        for molecule_name, molecule_spectra in self.individual_molecule_spectra_df.items(): 
            fig.add_trace(go.Scatter(x = self.bin_mz, 
                                     y = molecule_spectra,
                                     name = molecule_name,
                                     mode='lines',
                                     opacity= 0.6, 
                                     line = dict(width=1, dash = "solid")))

        fig.update_layout(title= f"{self.sample_name}, {self.primary_ion}, {self.replicat}",
                          height = 600,
                          template = "none",
                          showlegend=True,
                          yaxis=dict(type= yaxis_type),
                          margin=dict(t=50, b=50, l = 50, r = 50))



        fig.show()

