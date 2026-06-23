import numpy as np
from dataclasses import dataclass

@dataclass
class DataReader():
    """class to read in data from experimental files"""

    ### define initialization arguments ###
    files : list = None
    x_cuts : list = None
    y_cuts : list = None

    ### normalization arguments ###
    norm : list = None
    # whether to normalize at a specific x-value
    norm_at : list = None
    devide : list = None

    ### moving average ###
    ma : list = None
    ma_npoints : list = None

    ### formating options ###
    delimiter: list = None
    skiprows: list = None
    usecols : list = None

    ### baseline correction ###
    baseline : list = None
    baseline_at : list = None

    # whether to use TDM representation 
    TDM : bool = False

    # if you want to compute absorptance (for comparison to excitation spectra)
    absorptance : bool = False

    # whether to extract wavenumbers or wavelengths 
    wn : bool = True

    # whether data is IR or vis
    IR : bool = False

    # if emission data 
    em : bool = False
    corr: bool = True

    ### transient absorption data reading arguments ### 
    two_dim : bool = False
    experiment : str = 'femto'
    delay : list = None
    scatter : list = None
    wavelength : list = None
    slicing : list = None

    # temporary list to store normalization values for dA 
    movnorm : list = None

    def read_data(self) -> None:
        """method to read/cut/normalize data"""

        ### initialize list of data arrays for files ###
        self.x = np.empty(len(self.files), dtype=object)
        self.y = np.empty(len(self.files), dtype=object)
        if self.two_dim:
            self.z = np.empty(len(self.files), dtype=object)
        if self.ma:
            self.x_ma = np.empty(len(self.files), dtype=object)
            self.y_ma = np.empty(len(self.files), dtype=object)

        ### define standard values if not set ###
        # set_standard value method is implemented in plotter class
        self.delimiter = self.set_standard_value(self.delimiter, default=',')
        self.skiprows = self.set_standard_value(self.skiprows, default=2)
        self.usecols = self.set_standard_value(self.usecols, default=(0,1))
        self.x_cuts = self.set_standard_value(self.x_cuts, default=False)
        self.y_cuts = self.set_standard_value(self.y_cuts, default=False)
        self.norm = self.set_standard_value(self.norm, default=False)
        self.norm_at = self.set_standard_value(self.norm_at, default=False)
        self.devide = self.set_standard_value(self.devide, default=False)
        self.baseline = self.set_standard_value(self.baseline, default=False)
        self.baseline_at = self.set_standard_value(self.baseline_at, default=False)
        self.ma = self.set_standard_value(self.ma, default=False)
        self.ma_npoints = self.set_standard_value(self.ma_npoints, default=5)
        self.scatter = self.set_standard_value(self.scatter, default=False)
        self.slicing = self.set_standard_value(self.slicing, default=False)

        ### loop through files ###
        for i in range(len(self.files)):
            ### 2D Data = TA data ###
            if self.two_dim:
                if '.npy' in self.files[i]:
                    t, wl, dA = np.load(self.files[i], allow_pickle=True)
                if '.pdat' in self.files[i]:
                    from pyrene.standard.misc import load_pdat
                    t, wl, dA = load_pdat(self.files[i])
                # remove scatter 
                if self.scatter[i]:
                    dA[:, (wl >= self.scatter[i][0]) & (wl <= self.scatter[i][1])] = np.nan  
                # extract spectrum at delay
                if self.delay:
                    from pyrene.standard.misc import find_index
                    self.y[i] = dA[find_index(t, self.delay[i]), :]
                    self.x[i] = wl
                # extract kinetics at wavelength
                elif self.wavelength:
                    from pyrene.standard.misc import find_index
                    self.y[i] = dA[:, find_index(wl, self.wavelength[i])]
                    self.x[i] = t                    
                else:
                    self.x[i] = wl
                    self.y[i] = t
                self.z[i] = dA
            ### 1D Data = standard spectra ###
            else:
                if not '.txt' in self.files[i]:
                    # normal csv file
                    if not self.IR:
                        data = np.loadtxt(self.files[i], skiprows=self.skiprows[i], delimiter=self.delimiter[i], usecols=self.usecols[i])
                        self.x[i] = data[:,0]
                        self.y[i] = data[:,1]   
                        # correct baseline
                        if self.baseline[i]:
                            base = np.loadtxt(self.baseline[i], skiprows=self.skiprows[i], delimiter=self.delimiter[i], usecols=self.usecols[i]) 
                            self.y[i] -= base[:, 1] 
                    else:
                        # load FTIR file
                        from brukeropus import read_opus
                        opus_file = read_opus(self.files[i])
                        self.x[i] = np.array(opus_file.a.x)
                        self.y[i] = np.array(opus_file.a.y) 
                else:
                    # exported file
                    data = np.loadtxt(self.files[i], skiprows=1, delimiter=',')
                    if not self.IR:
                        self.x[i] = data[:,0]
                        self.y[i] = data[:,2]                       
                    else:
                        self.x[i] = data[:,1]
                        self.y[i] = data[:,-1]      

            # slice 2D countour along t if wanted
            if self.slicing[i]:
                self.z[i] = self.z[i][0::self.slicing[i],:]
                self.y[i] = self.y[i][0::self.slicing[i]]                                     

            # cut data if wanted
            if self.x_cuts[i]:
                mask = (self.x[i]>self.x_cuts[i][0])&(self.x[i]<self.x_cuts[i][1])
                self.x[i] = self.x[i][mask]
                if self.contour:
                    self.z[i] = self.z[i][:, mask]
                else:
                    self.y[i] = self.y[i][mask]
            if self.y_cuts[i]:
                mask = (self.y[i]>self.y_cuts[i][0])&(self.y[i]<self.y_cuts[i][1])
                self.y[i] = self.y[i][mask]
                if self.contour:
                    self.z[i] = self.z[i][mask, :]
                else:
                    self.x[i] = self.x[i][mask]

            # correct baseline using a single value
            if self.baseline_at[i]:
                from pyrene.standard.misc import find_index
                self.y[i] -= self.y[i][find_index(self.x[i], self.baseline_at[i])] 

            # photometric correction for emission spectra 
            if self.em and self.corr and not '.txt' in self.files[i]:
                import os
                corr_file = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'FMax_lamp_20151217.txt'))
                c = np.interp(self.x[i], corr_file[:,0], corr_file[:,1])
                self.y[i] *= c

            # convert intensity for emission I(nu) = I(lambda) * lambda^2
            if self.em and self.wn:
                self.y[i] *= self.x[i]**2

            # apply transient dipole moment representation if wanted 
            if self.TDM:
                wn = 1/self.x[i]
                if self.em: 
                    self.y[i] /= wn**3
                else:
                    self.y[i] /= wn

            # convert absorbance to absorptance 
            if self.absorptance: 
                self.y[i] = 1 - 10**(-self.y[i])

            # normalize if needed
            if self.norm[i]:
                if self.norm_at[i]:
                    from pyrene.standard.misc import find_index
                    if not self.contour:
                        self.y[i] /= np.abs(self.y[i][find_index(self.x[i], self.norm_at[i])])
                    else:
                        if self.movnorm:
                            self.movnorm[i] = np.abs(self.z[i][find_index(self.y[i], self.norm_at[i][1]), find_index(self.x[i], self.norm_at[i][0])]) 
                        self.z[i] /= np.abs(self.z[i][find_index(self.y[i], self.norm_at[i][1]), find_index(self.x[i], self.norm_at[i][0])])
                else:
                    if not self.contour:
                        self.y[i] /= np.abs(np.nanmax(self.y[i]))
                    else:
                        if self.movnorm:
                            self.movnorm[i] = np.abs(np.nanmax(self.z[i])) 
                        self.z[i] /= np.abs(np.nanmax(self.z[i]))

            if self.devide[i]:
                if not self.contour:
                    self.y[i] /= self.devide[i]
                else:
                    self.z[i] /= self.devide[i]

            # calculate moving average if wanted 
            if self.ma[i]:
                from pyrene.standard.misc import moving_average
                self.x_ma[i] = moving_average(self.x[i], self.ma_npoints[i])
                self.y_ma[i] = moving_average(self.y[i], self.ma_npoints[i])

            # convert to wavenumber if needed
            if self.wn and not self.IR:
                self.x[i] = 1e4/self.x[i]

    def __str__(self):
        return "class to read in files that contain data"