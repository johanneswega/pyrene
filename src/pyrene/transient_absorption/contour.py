from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class Contour(DataReader, Plotter, DataExporter):
    """class to plot, analyze and compare TA contours"""

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):
        ### read data ###
        self.two_dim = True
        self.contour = True
        self.read_data()
        self.set_contour_plot_settings()

        if self.norm or self.devide:
            self.zlabel = r'norm. $\Delta A$'
        else:
            self.zlabel = r'$\Delta A / 10^{-3}$'

        if self.experiment=='femto':
            self.ylabel = r'$\Delta t / \text{ps}$'
        else:
            self.ylabel = r'$\Delta t / \text{ns}$'

        if len(self.files)>1:
            self.fig, self.ax = plt.subplots(1, len(self.files), figsize=self.figsize, sharey=True)    

        for i in range(len(self.files)):
            if len(self.files)==1:
                self.plot_data()
            else:
                self.plot_data(master_ax=self.ax[i], contour_index=i)

    def show(self):
        if len(self.files)==1:
            self.show_plot(self.ax)
        else:
            for i in range(len(self.files)):
                if self.titles:
                    self.show_plot(self.ax[i], title=self.titles[i])
                else:
                    self.show_plot(self.ax[i])
        self.save_fig()
        plt.show()