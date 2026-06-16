from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class Absorption(DataReader, Plotter, DataExporter):
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
            if self.TDM:
                self.ylabel = r'$\varepsilon(\tilde{\nu}) / \tilde{\nu}$'
            else:
                self.ylabel = 'norm. absorbance'
        else:
            # calculate extinction coefficient A = eps c l --> eps = A / c l
            for i in range(len(self.files)):
                self.y[i] /= self.c[i] * self.l[i]
            self.ylabel = r'$\varepsilon / \text{M}^{-1}\,\text{cm}^{-1}$'
        
        if not self.slave:
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
            
    def plot_calculated_spectrum(self, logfile, cuts, fwhm, color, label, scale=None, shift=None, norm=True):
        import cclib
        # load data
        data = cclib.io.ccread(logfile)

        # IR = Vibrational frequencies (cm^-1) / vis = wavenumber (kK)
        if self.IR:
            xp = np.array(data.vibfreqs)
        else:
            xp = np.array(data.etenergies) / 1000

        # IR = intensities (km/mol) / vis = oscillator strengths
        if self.IR:
            yp = np.array(data.vibirs)
        else:
            yp = np.array(data.etoscs)

        # scale and shift
        label += '\n'
        if scale:
            xp = scale*xp
            label += f'scale = {scale}'
            label += '\n'
        if shift:
            xp += shift
            if self.IR:
                label += r'shift = $%.4g\,\text{cm}^{-1}$'%(shift)
                label += '\n'
            else:
                label += r'shift = $%.4g\,\text{cm}^{-1}$'%(shift*1000)
                label += '\n'

        if self.IR:
            label += r'FWHM = $%.4g\,\text{cm}^{-1}$'%(fwhm)
            label += '\n'
        else:
            label += r'FWHM = $%.4g\,\text{cm}^{-1}$'%(fwhm*1000)
            label += '\n'            

        # cut data
        yp = yp[(xp>cuts[0]) & (xp<cuts[1])]
        xp = xp[(xp>cuts[0]) & (xp<cuts[1])]

        # construct spectrum
        x = np.linspace(cuts[0], cuts[1], 5000)
        y = np.zeros(len(x))
        gamma = fwhm/2
        sigma = fwhm/(2*(2*np.log(2))**(0.5))
        for i in range(len(yp)):
            if self.IR:
                # compute intensity in (M-1 cm-1)
                y += ((100/np.log(10) * yp[i])/np.pi)*(gamma/(4*(x - xp[i])**2 + gamma**2))
            else:
                yp[i] = 1.3062974e8*yp[i]/(sigma*1000)
                y += yp[i]*np.exp(-((x - xp[i])**2) / (2*sigma**2))
        if norm:
            yp /= np.max(y)
            y /= np.max(y)

        # plot results
        self.ax.plot(x, y, '-', color=color, label=label)
        for i in range(len(yp)):
            self.ax.plot([xp[i], xp[i]], [0, yp[i] - 0.2*yp[i]], '-', color=color, alpha=0.5)
        
        return x, y    

    def solvchrom(self, lim=None, save=False):
        """method to plot absorption solvatochromism"""
        from pyrene.standard.solvents import solvent_dic
        fig1, ax1 = plt.subplots(1,1,figsize=(5, 3.5))
        fig2, ax2 = plt.subplots(1,1,figsize=(5, 3.5))
        fig3, ax3 = plt.subplots(1,1,figsize=(5, 3.5))
        for i in range(len(self.files)):
            n = solvent_dic[self.labels[i]][0]
            er = solvent_dic[self.labels[i]][1]
            fe = (2*(er - 1)/(2*er+1))
            fn = (2*(n**2 - 1)/(2*n**2+1))
            df = fe - fn
            ax1.plot(df, self.x[i][self.y[i]==np.max(self.y[i])], 'o', color=self.colors[i])
            ax2.plot(fe, self.x[i][self.y[i]==np.max(self.y[i])], 'o', color=self.colors[i])
            ax3.plot(fn, self.x[i][self.y[i]==np.max(self.y[i])], 'o', color=self.colors[i])
        ax1.set_ylabel(r'$\Tilde{\nu}_{\text{max, abs}} / 10^3$ cm$^{-1}$')
        ax2.set_ylabel(r'$\Tilde{\nu}_{\text{max, abs}} / 10^3$ cm$^{-1}$')
        ax3.set_ylabel(r'$\Tilde{\nu}_{\text{max, abs}} / 10^3$ cm$^{-1}$')
        ax1.set_xlabel(r'$\Delta f = f(\varepsilon_r) - f(n^2)$')
        ax2.set_xlabel(r'$f(\varepsilon_r)$')
        ax3.set_xlabel(r'$f(n^2)$')
        if lim!=None:
            ax1.set_ylim(lim)
            ax2.set_ylim(lim)
            ax3.set_ylim(lim)
        fig1.tight_layout()
        fig2.tight_layout()
        fig3.tight_layout()
        if save!=False:
            fig1.savefig('df.svg', transparent=True)
            fig2.savefig('f_epsr.svg', transparent=True)
            fig3.savefig('f_n2.svg', transparent=True)

    def epsilon_regression(self, M, m, Vstock, Vcuv, Vadd, l, wl):
        """
        Method to perform a linear regression from absorption spectra measured at different concentration

        M - molar mass in g/mol
        m - weighted masses in mg 
        Vstock - volume used to dissolve weighted masses in mL 
        Vcuv - volume [mL] of solution in which Vadd [uL] of stock is added
        l - pathlength [cm]
        wl - wavelength at which to analyze absorbance
        
        """

        from pyrene.standard.misc import find_index
        from scipy.optimize import curve_fit

        # get mass in g 
        m = m/1000
        # get volumia in L 
        Vstock = Vstock/1000
        Vadd = Vadd/10**6
        Vcuv = Vcuv/1000

        # calculate stock concentrations 
        cstock = m / (Vstock * M)
        # error for mass = 0.01 mg
        err_m = 0.01e-3
        # error of Vstock = 0.05 mL
        err_Vstock = 0.05e-3
        dcdm = 1 / (Vstock * M)
        dcdV = -m / (Vstock**2 * M) 
        cstockerr = np.sqrt( (dcdm * err_m)**2 + (dcdV * err_Vstock)**2)

        # print stock concentrations with error
        print('\n Stock concentrations: ')
        for i in range(len(cstock)):
            print(f'\n c = ({cstock[i]*1000 : .3f} +- {cstockerr[i]*1000 : .3f}) mM')

        # calculate diluted concentration 
        Vfull = Vcuv + Vadd
        err_Vcuv = (Vcuv/Vstock)*err_Vstock
        err_Vadd = 1e-6
        err_Vfull = np.sqrt(err_Vcuv**2 + err_Vadd**2)
        c = (cstock * Vadd) / Vfull 
        dcdc = Vadd/Vfull 
        dcdVa = cstock/Vfull 
        dcdVc = -(cstock * Vadd) / (Vfull **2)
        cerr = np.sqrt( (dcdc * cstockerr)**2 + (dcdVa * err_Vadd)**2 + (dcdVc * err_Vfull)**2 )

        # print concentration in cuvette with error
        print('\n Cuvette concentrations: ')
        for i in range(len(cstock)):
            print(f'\n c = ({c[i]*10**6 : .3f} +- {cerr[i]*10**6 : .3f}) µM')

        # do fit 
        fig, ax = plt.subplots(1,1, figsize=(5, 3.5))
        A = []
        for i in range(len(self.files)):
            A.append(self.y[i][find_index(self.x[i], wl)])
        ax.errorbar(x=A, y=c*10**6, yerr=cerr*10**6, ecolor='k', capsize=3, fmt='or')
        # fit to lamber beer law 
        p, pcov = curve_fit(lambda A, eps: A/(eps*l), A, c, sigma=cerr, absolute_sigma=True)
        eps = p[0]
        err_eps = pcov[0,0]**(0.5)
        Afine = np.linspace(A[0]-0.1, A[-1]+0.1, 100)
        ax.plot(Afine, Afine/(eps*10**-6 * l), '--k')
        s = f"{eps:.1e}"  
        power = int(s.split('e')[1])
        ax.set_title(r'$\varepsilon_{%i \text{nm}} = (%.2g \pm %.2g) \times 10^%i\,\text{M}^{-1}\,\text{cm}^{-1}$'%(wl, eps*10**-power, err_eps*10**-power, power))
        ax.set_xlabel(r'$A_{%i \text{nm}}$'%(int(round(wl))))
        self.ax.axvline(x = wl, color='k', linestyle='--', alpha=0.3)
        ax.set_ylabel(r'$c$ / µM')
        fig.tight_layout()
        fig.savefig('eps_regression.svg', transparent=True)
        print("")
        print(r"$\varepsilon_{%i \text{nm}} = (%.2g \pm %.2g) \times 10^%i\,\text{M}^{-1}\,\text{cm}^{-1}$"%(wl, eps*10**-power, err_eps*10**-power, power))
        print(f"relative error = {err_eps*100/eps : .2f} percent")
        print("")

    def show(self):
        self.show_plot()
        plt.show()

# test implementation
if __name__ == "__main__":
    files = ['../../../examples/steady_state/absorption/Data/abs_file1.csv',
             '../../../examples/steady_state/absorption/Data/abs_file2.csv']
    a = Absorption(files=files, x_cuts=[(300, 700), (300, 700)], 
                   labels=['1', '2'], yticks=False, baseline_at=[700, 700], 
                   colors=['b', 'r'])
    plt.show()