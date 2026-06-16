from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class Emission(DataReader, Plotter, DataExporter):
    """class to plot and analyze steady-state emission spectroscopy data"""

    # whether to do plot
    do_plot : bool = True

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):
        ### read data ###
        self.em = True
        self.skiprows = 2
        self.usecols = [0, 1]
        self.read_data()

        ### plot data ###
        if not self.norm:
            self.ylabel = 'intensity / counts'
        if self.TDM:
            self.ylabel = r'$F(\tilde{\nu}) / \tilde{\nu}^3$'
        else:
            self.ylabel = 'norm. intensity'

        if not self.slave:
            self.plot_data()

    def plot_absorption(self, AbsClass, ExClass=None):
        """method to plot absorption spectrum into the same figure"""
        AbsClass.plot_data(master_ax=self.ax)
        if ExClass:
            ExClass.plot_data(master_ax=self.ax)
        # make pretty axis labels
        self.ylabel = None
        self.yticks = False
        ax2 = self.ax.twinx()
        ax2.set_yticks([])
        if not self.TDM:
            self.ax.set_ylabel('norm. absorbance') 
            if ExClass:
                self.ax.set_ylabel('norm. absorbance / excitation') 
            ax2.set_ylabel('norm. intensity')
        else:
            self.ax.set_ylabel(r'$\varepsilon(\tilde{\nu}) / \tilde{\nu}$') 
            ax2.set_ylabel(r'$F(\tilde{\nu}) / \tilde{\nu}^3$')

    def show(self):
        self.show_plot()
        plt.show()