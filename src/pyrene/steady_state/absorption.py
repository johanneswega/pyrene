from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class Absorption(DataReader, Plotter):

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):
        self.skiprows = 2
        self.usecols = [0, 1]
        self.read_data()
        self.plot_data()

# test implementation
if __name__ == "__main__":
    files = ['../../../examples/steady_state/absorption/abs_file1.csv',
             '../../../examples/steady_state/absorption/abs_file2.csv']
    a = Absorption(files=files, x_cuts=[(300, 500), (300, 500)], 
                   labels=['1', '2'], 
                   colors=['b', 'r'], marker=['-', '--'], 
                   norm=[True, True], norm_at=[409, 350])
    plt.show()