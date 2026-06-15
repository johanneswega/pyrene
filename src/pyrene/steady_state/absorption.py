from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class Absorption(DataReader, Plotter):
    """class to plot and analyze steady-state absorption spectroscopy data"""

    # concentration list in M
    c : list = None
    # pathlength list in cm 
    l : list = None
    # provide extinction coefficient at certain wavelength [(eps_at_x, x)]
    eps : list = None

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):
        ### read data ###
        self.skiprows = 2
        self.usecols = [0, 1]
        self.read_data()

        ### plot data ###
        if not self.c:
            if not self.norm:
                self.ylabel = 'absorbance'
            else:
                self.ylabel = 'norm. absorbance'
        else:
            # calculate extinction coefficient A = eps c l --> eps = A / c l
            for i in range(len(self.files)):
                self.y[i] /= self.c[i] * self.l[i]
            self.ylabel = r'$\varepsilon / \text{M}^{-1}\,\text{cm}^{-1}$'
        
        self.plot_data()

    def calc_oscillator_strength(self, limits, n, nu0, file_index):
        """ 
        method to calculate oscillator strength and krad
        need to provide: limits in nm over which to calculate the band integral
        refractive index n, main transition freq. in kK
        """
        # plot integral region
        x = self.x[file_index][ (self.x[file_index] > limits[0]) & (self.x[file_index] < limits[1]) ]
        y = self.y[file_index][ (self.x[file_index] > limits[0]) & (self.x[file_index] < limits[1]) ] 
        self.ax.fill_between(x, 0, y, color='r', alpha=0.1)
        
        # flip arrays if not in acending order 
        if x[-1]<x[0]:
            x = np.flip(x)
            y = np.flip(y)

        # integrate 
        # convert to cm-1
        x = x*1000
        f = 4.3e-9 * np.trapezoid(y, x=x)
        # convert nu0 to from kK to s-1
        nu0 = sc.c*100 * nu0*1000
        krad = 7.4e-22 * ((n**2 + 2)/3) * nu0**2 * n * f
        tau_rad = (1/krad)*10**9
        # calculate TDM
        fL = 3*n**2/(2*n**2 + 1)
        TDM = 9.584e-2 * (n/fL**2)**(0.5) * (np.trapezoid(y/x, x=x))**(0.5)
        print("")
        print('The oscillator strenth calculated over the region is: f = %.3g'%(f))
        print(r'The TDM is: $\mu_{\text{TDM}} = %.3g\,\text{D}$'%(TDM))
        print(r'Radiative rate constant $k_{\text{rad}} = \SI{%.3g}{\per \second}$'%(krad))
        print(r'Radiative lifetime $\tau_\text{rad} = %.3g\,\text{ns}'%(tau_rad))
        print("")

    def find(self, file_index, where=None):
        """
        method to get absorbance at specified x-value (where) or at maximum (max==True)  (and extinction coefficient)
        """
        from pyrene.standard.misc import find_index
        if not where:
            # find wavelength at maximum absorbance
            Afound = np.max(self.y[file_index])
            xfound = self.x[file_index][self.y[file_index] == np.max(self.y[file_index])][0]
        else:
            Afound = self.y[file_index][find_index(self.x[file_index], where)]
            xfound = self.x[file_index][find_index(self.x[file_index], where)]
        if xfound>100:
            wlfound = xfound
            wnfound = 1e4/wlfound
        else:
            wnfound = xfound
            wlfound = 1e4/xfound
        if self.c:
            eps = Afound 
            print("")
            print(r"%s --> \varepsilon_{\text{%.3g nm}} = \SI{%.3g}{\text{M}^{-1}\,\text{cm}^{-1}} at wl = %.4g nm / wn = %.4g kK"%(self.files[file_index], wlfound, eps, wlfound, wnfound))
            print("")
        else:
            print("")
            print("%s --> A = %.3g at wl = %.4g nm / wn = %.4g kK"%(self.files[file_index], Afound, wlfound, wnfound))
            print("")
    
    def get_concentration(self, file_index) -> float:
        """ 
        to calculate concentration of solution
        need to procvide eps = [ (eps_at_x_look, x_look) ] at intialization
        """
        from pyrene.standard.misc import find_index
        x_look = self.eps[file_index][1]
        conc = self.y[file_index][find_index(self.x[file_index], x_look)] / (self.eps[file_index][0] * self.l[file_index])
        if 10**(-3) <= conc < 1:
            print("%s -> c = %.3g mM"%(self.files[file_index], conc*10**3))
        elif 10**(-6) <= conc < 10**(-3):
            print("%s -> c = %.3g µM"%(self.files[file_index], conc*10**6))
        else:
            print("%s -> c = %.3g M"%(self.files[file_index], conc))
        return conc

    def plot_diff(self, file1, file2, scale=False):
        """ 
        calculate and plot the difference spectrum between two files
        """
        x = self.x[file1]
        A1 = self.y[file1]
        A2 = self.y[file2]

        # interploate if needed
        if len(A1)!=len(A2):
            import numpy as np 
            if len(A1)>len(A2):
                x = self.x[file2]
                A1 = np.interp(x, self.x[file1], A1)
            else:
                x = self.x[file1]
                A2 = np.interp(x, self.x[file2], A2)  

        if scale!=False:
            dA = A1*scale[0] - A2*scale[1]
            self.ax.plot(x, dA, '--k', label=r'Difference $\times$ %i'%(round(scale[1])))
        else:
            dA = A1 - A2
            self.ax.plot(x, dA, '--k', label=r'Difference')
            
    def show(self):
        self.show_plot()
        plt.show()

# test implementation
if __name__ == "__main__":
    files = ['../../../examples/steady_state/absorption/abs_file1.csv',
             '../../../examples/steady_state/absorption/abs_file2.csv']
    a = Absorption(files=files, x_cuts=[(300, 700), (300, 700)], 
                   labels=['1', '2'], yticks=False, baseline_at=[700, 700], 
                   colors=['b', 'r'])
    plt.show()