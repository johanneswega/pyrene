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

    ### formating options ###
    delimiter: str = ','
    skiprows: int = 0
    usecols : list = None

    # whether 1 D or 2 D data 
    two_dim : bool = False

    def __post_init__(self):
        """post initialization to fill in standard arguments"""
        self.ma_npoints = [5 for _ in self.files]

    def read_data(self) -> None:
        """method to read/cut/normalize data"""

        ### initialize list of data arrays for files ###
        self.x = np.empty(len(self.files), dtype=object)
        self.x_inv = np.empty(len(self.files), dtype=object)
        self.y = np.empty(len(self.files), dtype=object)
        if self.two_dim:
            self.z = np.empty(len(self.files), dtype=object)

        ### loop through files ###
        for i in range(len(self.files)):

            # read data
            if self.two_dim:
                ...
            else:
                data = np.loadtxt(self.files[i], skiprows=self.skiprows, delimiter=self.delimiter, usecols=self.usecols)
                self.x[i] = data[:,0]
                self.y[i] = data[:,1]             

            # cut data if wanted
            if self.x_cuts:
                mask = (self.x[i]>self.x_cuts[i][0])&(self.x[i]<self.x_cuts[i][1])
                self.y[i] = self.y[i][mask]
                self.x[i] = self.x[i][mask]
            if self.y_cuts:
                mask = (self.y[i]>self.y_cuts[i][0])&(self.y[i]<self.y_cuts[i][1])
                self.x[i] = self.x[i][mask]
                self.y[i] = self.y[i][mask]

            # normalize if needed
            if self.norm:
                if self.norm_at:
                    from pyrene.standard.misc import find_index
                    self.y[i] /= self.y[i][find_index(self.x[i], self.norm_at[i])]     
                else:
                    self.y[i] /= np.max(self.y[i])  

            # calculate inverse of x-axis
            self.x_inv[i] = ((1/self.x[i])*1e4)     

    def __str__(self):
        return "class to read in files that contain data"