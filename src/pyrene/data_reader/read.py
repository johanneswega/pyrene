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
    min_norm : list = None
    devide : list = None

    ### moving average ###
    ma : list = None
    ma_npoints : list = None

    ### formating options ###
    delimiter: str = ','
    skiprows: int = 0
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
        if self.x_cuts:
            if len(self.x_cuts)==1:
                self.x_cuts = [self.x_cuts[0] for _ in self.files]
        if self.y_cuts:
            if len(self.y_cuts)==1:
                self.y_cuts = [self.y_cuts[0] for _ in self.files]
        if self.norm:
            if len(self.norm)==1:
                self.norm = [self.norm[0] for _ in self.files]    
        if self.norm_at:
            if len(self.norm_at)==1:
                self.norm_at = [self.norm_at[0] for _ in self.files]         
        if not self.norm: 
            self.norm = [False for _ in self.files]
        if not self.min_norm:
            self.min_norm = [False for _ in self.files]
        if not self.devide:
            self.devide = [False for _ in self.files]
        if not self.ma: 
            self.ma = [False for _ in self.files]
        else:
            if len(self.ma)==1:
                self.ma = [True for _ in self.files]
        if self.ma and not self.ma_npoints:
            self.ma_npoints = [5 for _ in self.files]
        if not self.scatter:
            self.scatter = [False for _ in self.files]
        else:
            if len(self.scatter)==1:
                self.scatter = [self.scatter[0] for _ in self.files]

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
                else:
                    self.x[i] = wl
                    self.y[i] = t
                self.z[i] = dA
            ### 1D Data = standard spectra ###
            else:
                if not '.txt' in self.files[i]:
                    # normal csv file
                    if not self.IR:
                        data = np.loadtxt(self.files[i], skiprows=self.skiprows, delimiter=self.delimiter, usecols=self.usecols)
                        self.x[i] = data[:,0]
                        self.y[i] = data[:,1]   
                        # correct baseline
                        if self.baseline:
                            if self.baseline[i]:
                                base = np.loadtxt(self.baseline[i], skiprows=self.skiprows, delimiter=self.delimiter, usecols=self.usecols) 
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

            # cut data if wanted
            if self.x_cuts:
                mask = (self.x[i]>self.x_cuts[i][0])&(self.x[i]<self.x_cuts[i][1])
                self.x[i] = self.x[i][mask]
                if self.contour:
                    self.z[i] = self.z[i][:, mask]
                else:
                    self.y[i] = self.y[i][mask]
            if self.y_cuts:
                mask = (self.y[i]>self.y_cuts[i][0])&(self.y[i]<self.y_cuts[i][1])
                self.y[i] = self.y[i][mask]
                if self.contour:
                    self.z[i] = self.z[i][mask, :]
                else:
                    self.x[i] = self.x[i][mask]

            # correct baseline using a single value
            if self.baseline_at:
                from pyrene.standard.misc import find_index
                if self.baseline_at[i]:
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
                if self.norm_at:
                    from pyrene.standard.misc import find_index
                    if not self.contour:
                        self.y[i] /= self.y[i][find_index(self.x[i], self.norm_at[i])]
                    else:
                        self.z[i] /= self.z[i][find_index(self.y[i], self.norm_at[i][0]), find_index(self.x[i], self.norm_at[i][1])]     
                else:
                    if not self.contour:
                        self.y[i] /= np.max(self.y[i]) 
                    else:
                        self.z[i] /= np.max(self.z[i])

            if self.min_norm[i]:
                self.y[i] *= -1 

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