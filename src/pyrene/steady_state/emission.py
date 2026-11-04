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
        wns_em = np.array([self.x[i][self.y[i]==np.max(self.y[i])][0] for i in range(len(self.files))])

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

    def calc_FluQY(self, file_ind_S, file_ind_R, A_S, A_R, phi_R, solvents):
        """method to calculate fluorescence quantum yield based on the relative method
        A_R, A_S : absorbances of reference and sample at their respective excitation wavelengths"""

        # get refractive indices
        from pyrene.standard.solvents import solvent_dic
        n_S, n_R = solvent_dic[solvents[0]][0], solvent_dic[solvents[1]][0]
        print(f'\n refractive indices at 25°C of \n nS ({solvents[0]}) = {n_S} \n nR ({solvents[1]}) = {n_R} \n')

        # integrate fluorescence
        I_S = np.trapezoid(self.y[file_ind_S], x=self.x[file_ind_S])
        I_R = np.trapezoid(self.y[file_ind_R], x=self.x[file_ind_R])
        print(I_S, I_R)

        # calculate fluoresence quantum yield
        fS = 1 - 10**(-A_S)
        fR = 1 - 10**(-A_R)
        phiQflu = phi_R * (fR/fS) * (I_S / I_R) * (n_S/n_R)**2

        ### error estimation ### 
        # take error in Aex as 0.002
        errAex = 0.002
        # take error in intensities as relative error of 2 % 
        errint_S = 0.02*I_S
        errint_R = 0.02*I_R
        # take error in refractive indicies as 0.0002
        errn = 0.0002
        # take error in reference quantum yield as 2 % 
        errphiR = 0.02

        # partial derivatives
        dphiR = 10**(-A_R + A_S)*I_S*n_S**2*(10**A_R - 1)/(I_R*n_R**2*(10**A_S - 1))
        dI_S = 10**(-A_R + A_S)*n_S**2*phi_R*(10**A_R - 1)/(I_R*n_R**2*(10**A_S - 1))
        dI_R = -10**(-A_R + A_S)*I_S*n_S**2*phi_R*(10**A_R - 1)/(I_R**2*n_R**2*(10**A_S - 1))
        dA_R = 10**(-A_R + A_S)*I_S*n_S**2*phi_R*np.log(10)/(I_R*n_R**2*(10**A_S - 1))
        dA_S = -10**(-A_R + A_S)*I_S*n_S**2*phi_R*(10**A_R - 1)*np.log(10)/(I_R*n_R**2*(10**A_S - 1)**2)
        dn_S = 2*10**(-A_R + A_S)*I_S*n_S*phi_R*(10**A_R - 1)/(I_R*n_R**2*(10**A_S - 1))
        dn_R = -2*10**(-A_R + A_S)*I_S*n_S**2*phi_R*(10**A_R - 1)/(I_R*n_R**3*(10**A_S - 1))

        # calculate error 
        errQflu = np.sqrt( (dphiR * errphiR)**2 + (dI_S * errint_S)**2 + (dI_R * errint_R)**2 + (dA_R * errAex)**2 + (dA_S * errAex)**2 + (dn_S * errn)**2 + (dn_R * errn)**2 )
        print(f'\n Fluorescence Quantum Yield = ({phiQflu*100 : .3f} +- {errQflu*100 : .3f}) percent \n ({100*errQflu/phiQflu : .3f} perc. rel error)')

    def show(self):
        self.show_plot(self.ax)
        self.save_fig()
        plt.show()