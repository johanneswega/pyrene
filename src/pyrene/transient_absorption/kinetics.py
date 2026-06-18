from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class Kinetics(DataReader, Plotter, DataExporter):
    """class to plot, analyze and compare spectral slices at different delays"""

    # wheter to do rainbow plot (overview)
    overview : bool = False

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):

        # get axis label
        if self.norm or self.devide:
            self.ylabel = r'norm. $\Delta A$'
        else:
            self.ylabel = r'$\Delta A / 10^{-3}$'
        if self.experiment=='femto':
            self.xlabel = r'$\Delta t / \text{ps}$'
        else:
            self.xlabel = r'$\Delta t / \text{ns}$'

        ### read data ###
        self.wn = False
        self.two_dim = True
        self.read_data()
        self.plot_data()

    def show(self):
        self.show_plot(self.ax)
        self.save_fig()
        plt.show()