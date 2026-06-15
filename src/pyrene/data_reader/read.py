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

        ### loop through files ###
        for i in range(len(self.files)):

            # read data
            if self.two_dim:
                ...
            else:
                data = np.loadtxt(self.files[i], skiprows=self.skiprows, delimiter=self.delimiter, usecols=self.usecols)
                self.x[i] = data[:,0]
                self.y[i] = data[:,1]   
                if self.baseline:
                    if self.baseline[i]:
                        base = np.loadtxt(self.baseline[i], skiprows=self.skiprows, delimiter=self.delimiter, usecols=self.usecols) 
                        self.y[i] -= base[:, 1]         

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
            if self.norm:
                if self.norm_at:
                    from pyrene.standard.misc import find_index
                    self.y[i] /= self.y[i][find_index(self.x[i], self.norm_at[i])]     
                else:
                    self.y[i] /= np.max(self.y[i])  

            # calculate moving average if wanted 
            if self.ma:
                if self.ma[i]:
                    from pyrene.standard.misc import moving_average
                    self.x_ma[i] = moving_average(self.x[i], self.ma_npoints[i])
                    self.y_ma[i] = moving_average(self.y[i], self.ma_npoints[i])
            else:
                self.ma = [None for _ in self.files]

    def __str__(self):
        return "class to read in files that contain data"