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
    norm : bool = False
    # whether to normalize at a specific x-value
    norm_at : float = None

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

    # whether 1 D or 2 D data 
    two_dim : bool = False

    # whether to extract wavenumbers or wavelengths 
    wn : bool = True

    # whether data is IR or vis
    IR : bool = False

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

        ### standard values of not set ###
        if self.x_cuts:
            if len(self.x_cuts)==1:
                self.x_cuts = [self.x_cuts[0] for _ in self.files]
        if self.y_cuts:
            if len(self.y_cuts)==1:
                self.y_cuts = [self.y_cuts[0] for _ in self.files]
        if self.norm:
            if len(self.norm)==1:
                self.norm = [self.norm[0] for _ in self.files]            
        if not self.norm: 
            self.norm = [False for _ in self.files]
        if not self.ma: 
            self.ma = [False for _ in self.files]
        if self.ma and not self.ma_npoints:
            self.ma_npoints = [5 for _ in self.files]

        ### loop through files ###
        for i in range(len(self.files)):

            # read data
            if self.two_dim:
                ...
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
                        self.y[i] = data[:,-1]    
                    else:
                        self.x[i] = data[:,1]
                        self.y[i] = data[:,-1]                             

            # cut data if wanted
            if self.x_cuts:
                mask = (self.x[i]>self.x_cuts[i][0])&(self.x[i]<self.x_cuts[i][1])
                self.y[i] = self.y[i][mask]
                self.x[i] = self.x[i][mask]
            if self.y_cuts:
                mask = (self.y[i]>self.y_cuts[i][0])&(self.y[i]<self.y_cuts[i][1])
                self.x[i] = self.x[i][mask]
                self.y[i] = self.y[i][mask]

            # correct baseline using a single value
            if self.baseline_at:
                from pyrene.standard.misc import find_index
                if self.baseline_at[i]:
                    self.y[i] -= self.y[i][find_index(self.x[i], self.baseline_at[i])] 

            # normalize if needed
            if self.norm[i]:
                if self.norm_at:
                    from pyrene.standard.misc import find_index
                    self.y[i] /= self.y[i][find_index(self.x[i], self.norm_at[i])]     
                else:
                    self.y[i] /= np.max(self.y[i])  

            # calculate moving average if wanted 
            if self.ma[i]:
                from pyrene.standard.misc import moving_average
                self.x_ma[i] = moving_average(self.x[i], self.ma_npoints[i])
                self.y_ma[i] = moving_average(self.y[i], self.ma_npoints[i])

            # invert if needed
            if self.wn and not self.IR:
                self.x[i] = 1e4/self.x[i]

    def __str__(self):
        return "class to read in files that contain data"