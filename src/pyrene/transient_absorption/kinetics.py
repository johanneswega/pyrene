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

    def fit(self, file_index=0, model=None, p0=None):
        """method to fit kinetic trace"""
        x = self.x[file_index]
        y = self.y[file_index]
        xfine = np.linspace(np.min(x), np.max(x), 1000)
        p, pcov = curve_fit(model, x, y, p0=p0)
        self.ax.plot(xfine, model(xfine, *p), '-', color=self.colors[file_index], label=self.labels[file_index] + ' fit')

        # plot residuals
        fig, ax = plt.subplots(1,1, figsize=(7, 3))
        ax.plot(x, model(x, *p) - y, '-', color=self.colors[file_index])
        ax.set_ylabel('residuals')
        ax.set_xscale(self.xscale)
        ax.set_yscale(self.yscale)
        ax.axhline(y=0, color='k')
        ax.set_xlabel(self.xlabel)
        fig.tight_layout()

        # print results
        print("")
        print(f"### fitted file {self.files[file_index]} with {model.__name__} ###")
        print("")
        if "exp" in model.__name__:
            from pyrene.standard.misc import print_results_of_exp_fit
            print_results_of_exp_fit(self, model, p, pcov)

    def show(self):
        self.show_plot(self.ax)
        self.save_fig()
        plt.show()