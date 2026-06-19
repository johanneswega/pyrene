from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class Spectra(DataReader, Plotter, DataExporter):
    """class to plot, analyze and compare spectral slices at different delays"""

    # wheter to do rainbow plot (overview)
    overview : bool = False
    # for different steady-state spectra in the different axes in subplots
    # simple put the normal steady state spectra in yet another list 
    # e.g. steady_state_ax = [[s1, s2], [s3]]
    # s1 and s2 will be plotted in ax[0] while s3 will be plotted in ax[1]
    steady_state_ax : list = None
    export_overview : bool = False

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):

        if self.norm or self.devide:
            self.ylabel = r'norm. $\Delta A$'
        else:
            self.ylabel = r'$\Delta A / 10^{-3}$'
    
        ### read data ###
        self.two_dim = True
        if not self.overview:
            self.read_data()
            self.plot_data()
        else:
            ### for overview plot ###
            self.nfiles = self.files
            devide = self.devide
            if len(self.nfiles)>1:
                self.fig, self.ax = plt.subplots(1, len(self.nfiles), figsize=self.figsize, sharey=True)
            for i in range(len(self.nfiles)):
                if not self.delay:
                    if self.experiment=='femto':
                        self.delay = [0.2, 0.5, 1, 2, 5, 10, 50, 100, 250, 500, 1000, 1200, 1800]
                    else:
                        self.delay = [2, 5, 10, 20, 50, 80, 100, 250, 200, 500, 1e+3, 5e+3, 20e+3, 50e+3, 100e+3, 200e+3, 500e+3]
                if devide:
                    self.devide = [devide[i] for _ in self.delay]
                if self.steady_state_ax:
                    self.steady_state = self.steady_state_ax[i]
                self.files = [self.nfiles[i] for _ in self.delay]
                self.labels = self.get_delay_labels(self.delay)
                from pyrene.standard.misc import rainbow
                self.colors = rainbow(self.delay, r=True)
                self.read_data()
                if len(self.nfiles)==1:
                    self.plot_data()
                else:
                    self.plot_data(master_ax=self.ax[i])
                if self.export_overview:
                    self.export()

    def show(self):
        if self.overview:
            if len(self.nfiles)==1:
                self.show_plot(self.ax)
            else:
                for i in range(len(self.nfiles)):
                    if self.titles:
                        self.show_plot(self.ax[i], title=self.titles[i])
                    else:
                        self.show_plot(self.ax[i])
        else:
            self.show_plot(self.ax)
        self.save_fig()
        plt.show()