from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class Emission(DataReader, Plotter, DataExporter):
    """class to plot and analyze steady-state emission spectroscopy data"""

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):
        ### read data ###
        self.em = True
        self.skiprows = 2
        self.usecols = [0, 1]
        self.read_data()

        ### plot data ###
        if not self.norm[0]:
            self.ylabel = 'intensity / counts'
        if self.devide[0]:
            self.ylabel = 'norm. intensity'
        if self.TDM:
            self.ylabel = r'$F(\tilde{\nu}) / \tilde{\nu}^3$'
        else:
            self.ylabel = 'norm. intensity'

        if not self.slave:
            self.plot_data()

        # for export 
        self.em_spec = True

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

    def solv_chrom(self, AbsClass=None, lim=None, save=False):
        """method to plot emission solvatochromism"""
        from pyrene.standard.figures import solvchrom_figures          
        
        # get wavenumber at maxima
        wns_em = np.array([self.x[i][self.y[i]==np.max(self.y[i])] for i in range(len(self.files))])

        # plot and get solvent parameters
        n, er, fe, fn, df = solvchrom_figures(wns_em, self.labels, self.colors, label='em', 
                                              stokes=False, lim=lim, save=save)
        
        # if absorption file provided perform abs solvchrom 
        if AbsClass: 
            wns_abs, n, er, fe, fn, df = AbsClass.solvchrom(lim=lim, save=save) 
            # calculate Stokes shift in cm-1 
            stokes_shift = (wns_abs - wns_em)*1000
            solvchrom_figures(stokes_shift, self.labels, self.colors, label='stokes', 
                                                stokes=True, lim=lim, save=save)
            return wns_abs, wns_em, stokes_shift, n, er, fe, fn, df            
        else:
            return wns_em, n, er, fe, fn, df  

    def show(self):
        self.show_plot(self.ax)
        self.save_fig()
        plt.show()