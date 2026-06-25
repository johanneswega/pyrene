from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass
from scipy.signal import convolve as conv

@dataclass
class TCSPC(DataReader, Plotter, DataExporter):
    """class to plot and analyze steady-state absorption spectroscopy data"""

    # compare different traces in the same plot
    compare : bool = False
    # wether to plot also the fitted traces on top in the comparison plot
    plot_fit : bool = False

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):
        self.experiment = 'nano'
        self.wn = False
        self.tcspc = True
        self.delimiter = [None]
        self.skiprows = [10]
        self.skipfooter = [2]
        # temporary list used for fit normalization
        self.norm_value = []

        ### read data ###
        if not self.compare:
            self.read_data()
            self.ylabel = r'counts'
        else:
            self.norm = [True]
            self.read_data()
            self.ylabel = r'norm. intensity'
            if self.plot_fit:
                self.alphas = [0.2]

        ### plot data ###
        self.xlabel = r'$\Delta t / \text{ns}$'
        self.plot_data()

        ### plot fit if wanted ###
        if self.plot_fit:
            for i in range(len(self.files)):
                if not 'IRF' in self.files[i]:
                    fit_file = self.files[i][:self.files[i].find('.')] + '_fit.txt'
                    data = np.loadtxt(fit_file, skiprows=1, delimiter=',')
                    x = data[:,0]
                    y = data[:,1]
                    y-=y[0]
                    y/=self.norm_value[i]
                    self.ax.plot(x, y, '-', color=self.colors[i], label=self.labels[i], linewidth=2)

    def get_fit_function(self, t, *params):
        """
        generates fit function used in the fit method
        params: A1, tau1, A2, tau2, ..., An, taun, shift
        """

        if not self.tail:
            shift = params[-1]
            params = params[:-1]
            # interpolate IRF
            IRF_int = np.interp(t, t + shift, self.irf)

        fit = np.zeros_like(t, dtype=float)
        # pairs: amplitude, lifetime
        for A, tau in zip(params[::2], params[1::2]):
            decay = np.exp(-t / tau)
            if not self.tail:
                conv_decay = conv(decay, IRF_int, mode="full")[:len(t)]
                # normalize each component like your original code
                conv_decay /= np.max(np.abs(conv_decay))
                fit += A * conv_decay
            else:
                fit += A * decay

        return fit + self.bg

    def fit(self, file_index=0, model=None, p0=None, tail=None, IRF=False, IRF_cleanup=0):
        """method to fit kinetic trace"""

        fig, ax = plt.subplots(3,1, figsize=(9, 7), sharex=True, gridspec_kw={'height_ratios':[1,1,3]}) 

        # for tail fit
        x = self.x[file_index]
        y = self.y[file_index]

        # only take counts larger than zero
        mask = y>0
        x = x[mask]
        y = y[mask]

        # calculate background 
        self.bg = np.nanmean(y[x<-3.0])
        ax[2].axhline(y=self.bg, linestyle='--', color='k')

        # plot raw data
        ax[2].plot(x, y, '-b', alpha=0.8, linewidth=1, label=r'\text{data} (%s)'%(self.labels[file_index]))

        if tail:
            self.tail = True
            y = y[x>tail]
            x = x[x>tail]
            ax[2].axvline(x=tail, linestyle='--', color='k')
        if IRF: 
            self.tail = False
            IRF_data = self.y[-1]
            IRF_data = IRF_data[mask]
            IRF_data /= np.max(IRF_data)
            IRF_data *= np.max(y)
            IRF_data[IRF_data<IRF_cleanup] = 0
            ax[2].fill_between(x, 0, IRF_data, color='k', alpha=0.1, label=r'IRF')       

        # get uncertainty
        sigma = np.sqrt(y)

        # fit trace
        try:
            if IRF:
                self.irf = IRF_data
                model = self.get_fit_function
                p, pcov = curve_fit(model, x, y, p0=p0, sigma=sigma, absolute_sigma=True)
            else:
                model = self.get_fit_function
                p, pcov = curve_fit(model, x, y, p0=p0, sigma=sigma, absolute_sigma=True)
        except Exception as e:
            print("")
            print(f"### fit failed for file {self.files[file_index]}###")
            print(f"Error: {e}")
            print("")
            return [], [], 1e3

        # compute red. chi2 and plot fir
        chi2 = (1/(len(y) - len(p)))*np.sum(((y - model(x, *p))/sigma)**2)
        ax[2].plot(x, model(x, *p), '-r', linewidth=2, label='fit')

        # save fit as .txt
        np.savetxt('%s_fit.txt'%self.files[file_index][:self.files[file_index].find('.')], 
                   np.column_stack([x, model(x, *p)]),
                   header='time delay / ns, fit / counts', delimiter=',')

        # plot residuals
        res = (model(x, *p) - y)/sigma
        ax[1].plot(x, res, '-r', label=r'$\chi^2_\nu = %.3g$'%(chi2), linewidth=0.5)

        # plot autocorrelation
        from pyrene.standard.misc import acf
        acf_res = acf(res)
        dt = np.mean(np.diff(x))
        lags = np.arange(len(acf_res)) * dt
        ax[0].plot(lags, acf_res, '-g', linewidth=0.5)
        ax[0].axhline(y=0, color='k')
        ax[0].set_ylim([-0.15, 0.15])
        ax[0].set_ylabel(r'autocorr')

        # stylistic stuff
        ax[1].set_ylabel(r'res / $\sigma$')
        ax[1].axhline(y=0, color='k')
        ax[2].set_ylim(self.ylim)
        ax[2].set_xlabel(self.xlabel)
        ax[2].set_ylabel(self.ylabel)
        ax[2].set_yscale(self.yscale)
        ax[2].set_xscale(self.xscale)
        ax[1].legend()
        ax[2].legend()

        # print results
        print("")
        print(f"### fit results for {self.files[file_index]} ###")
        print("")

        from pyrene.standard.misc import print_results_of_exp_fit
        if IRF:
            p = p[:-1]
            pcov = pcov[:-1,:-1]
        label, tau, tau_errs, amp_frac, area_frac = print_results_of_exp_fit(self, model, p, pcov, nomodel=True)
        ax[2].set_title(label[:-1], size=14)
        fig.tight_layout()
        fig.savefig('%s_fit.svg'%self.files[file_index][:self.files[file_index].find('.')], transparent=True)

        # save fit parameters as .txt
        np.savetxt('%s_fit_parameters.txt'%self.files[file_index][:self.files[file_index].find('.')], 
                   np.column_stack([tau, tau_errs, amp_frac, area_frac]),
                   header='component lifetime / ns, error in lifetime / ns, amplitude fraction, area fraction', delimiter=',')
    
    def fit_IRF(self):
        """method to fit IIRF"""
        fig, ax = plt.subplots(1,1, figsize=(5, 3.5))
        from pyrene.fitting.fit_functions import gaussian
        p, pcov = curve_fit(gaussian, self.x[-1], self.y[-1], p0=[1e3, 0, 0.5])
        FWHM = p[-1]*1000
        errFWHM = (pcov[-1, -1])**(0.5)*1000
        ax.plot(self.x[-1], self.y[-1], 'ok', alpha=0.3)
        tfine = np.linspace(-1.2, 1.2, 1000)
        ax.plot(tfine, gaussian(tfine, *p), '-r')
        ax.set_title(r'FWHM = $(%.3g \pm %.3g)\,\text{ps}$'%(FWHM, errFWHM))
        ax.set_ylabel('counts')
        ax.set_xlabel(r'$\Delta t$ / ns')
        ax.set_xlim([-1.2, 1.2])
        fig.tight_layout()

    def show(self):
        self.show_plot(self.ax)
        self.save_fig()
        plt.show()
