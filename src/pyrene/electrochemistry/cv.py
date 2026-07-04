from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class CV(DataReader, Plotter, DataExporter):
    """class to plot and analyze cyclic voltammetry data"""

    # list of E0(Fc+/Fc) for x-axis calibaration
    Fc : list = None
    # shifts Fc+/Fc calibrated CVs such that x-axis corresponds to vs. SCE
    SCE : bool = False
    # for US current convention 
    US: bool = False

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):
        self.CV = True
        self.wn = False
        self.zeroline = False

        if self.norm or self.devide:
            self.ylabel = 'normalized current'
        else:
            self.ylabel = 'current / µA' 
               
        self.read_data()
        
        # calibrate axis
        if self.Fc:
            for i in range(len(self.files)):
                if not self.SCE:
                    self.x[i] -= self.Fc[i]
                    self.xlabel = r'potential / V vs. Fc$^+$/Fc'
                else:
                    self.x[i] = self.x[i] - self.Fc[i] + 0.38
                    self.xlabel = r'potential / V vs. SCE'      
        else:
            self.xlabel = r'potential / V vs. Ag$^+$/Ag'            

        self.plot_data()

    def get_half_wave_potential(self, file_index, E_range, plot=False):
        """function to find peak and potentially also E0"""
        # file index and range
        E = self.x[file_index]
        I = self.y[file_index]
        # cut in desried voltage range
        I = I[(E>=E_range[0])&(E<=E_range[1])]
        E = E[(E>=E_range[0])&(E<=E_range[1])]
        ipa = np.max(I)
        ipc = np.min(I)
        print('ipc = %.5g µA'%(ipc))
        print('ipa = %.5g µA'%(ipa))
        Epa = E[np.argmax(I)]
        Epc = E[np.argmin(I)]
        E0 = (np.array(Epa) + np.array(Epc))/2
        print('')
        print(f'E0 = {E0 : .4f} V')
        if plot==True:
            self.ax.plot(E, I, '-g')
            self.ax.plot(Epa, ipa, 'ok')
            self.ax.plot(Epc, ipc, 'ok')
            self.ax.axvline(x=E0, color='k', linestyle='--')
            self.ax.set_title(r'$E_{1/2} = %.3g\,\text{V}$'%E0)
        return E0, ipc, ipa
    
    def peak_separation(self, scan_rates, E_range):
        """method to plot peak separation as function of the scan rate"""
        fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
        for i in range(len(scan_rates)):
            E = self.x[i]
            I = self.y[i]
            # cut in desried voltage range
            I = I[(E>=E_range[0])&(E<=E_range[1])]
            E = E[(E>=E_range[0])&(E<=E_range[1])]      
            Epa = E[np.argmax(I)]
            Epc = E[np.argmin(I)]  
            ax.plot(scan_rates[i], np.abs(Epa - Epc)*1e3, 'o', markerfacecolor=self.colors[i], markeredgecolor='k')
        ax.set_ylabel(r'$\Delta E_\text{p}$ / mV')
        ax.set_xlabel(r'$v / \text{V} \, \text{s}^{-1}$')
        ax.set_ylim([35, 155])
        fig.tight_layout()
        fig.savefig('peak_sep.svg', transparent=True)

    @staticmethod
    def cottrell(t, k, tprime, m, b):
        """used for cottrell fit below"""
        return k/((t - tprime)**0.5)  + m*t + b

    def cottrell_fit(self, scan_rates, E_range, plot=True):
        """method to perform a cotrell-fit to obtain the peak current ration |ipf/ipr| as a function of the scan rate"""
        ratio = []
        # make big overview figure
        if plot==True:
            figbig, axbig = plt.subplots(len(scan_rates), 1, figsize=(5, len(scan_rates)*3))
        # go through all scan rates
        for i in range(len(scan_rates)):
            v = scan_rates[i]
            E = self.x[i]
            I = self.y[i]
            # cut in desried voltage range
            I = I[(E>=E_range[0])&(E<=E_range[1])]
            E = E[(E>=E_range[0])&(E<=E_range[1])]    
            # get scan direction
            Estart = E[0]
            if E[1]<E[5]:
                scan_dir = 'ox'
                # get switching potential
                Elam = np.max(E)
            # reductive scan
            else:
                scan_dir = 'red'
                # get switching potential
                Elam = np.min(E) 
            # calculate switching time
            tlam = np.abs(Elam-Estart)/v
            # potentialstep
            dE = np.abs(E[1]-E[0])
            # calculate time axis
            t = np.array([i*dE/v for i in range(len(E))])
            # find forward peak
            if scan_dir=='ox':
                tpf = t[I==np.max(I)]
                Ipf0 = I[I==np.max(I)]
            else:
                tpf = t[I==np.min(I)]
                Ipf0 = I[I==np.min(I)]    
            # fit linear baseline on the first 20 points
            m, b = np.polyfit(t[:20], I[:20], 1)
            # make fine axis to extrapolate baseline
            lin_for = np.linspace(0, tpf, 1000)
            # calculate baseline
            I_for = m*lin_for + b
            # find time at peak potential on baseline
            linpf = I_for[np.argmin(np.abs(lin_for - tpf))]
            # calculate actual peak current
            Ipf = Ipf0 - linpf

            # specify fitting range for Cottrell-fit
            # as +0.05 V above peak and -0.05 V before switching potential
            lim = 0.05/v
            tfit = t[(t>tpf+lim)&(t<tlam-lim)]
            Ifit = I[(t>tpf+lim)&(t<tlam-lim)]

            # find peak values for back peak
            if scan_dir=='ox':
                tpb = t[I==np.min(I)]
                Ipb0 = I[I==np.min(I)]
            else:
                tpb = t[I==np.max(I)]
                Ipb0 = I[I==np.max(I)]    

            # do curve fit to obtain cottrell fit
            p, pcov = curve_fit(lambda t, k, tprime: self.cottrell(t, k, tprime, m, b), tfit, Ifit, p0=[Ipf0, tpf])
            # make axis to extrapolate cottrell fit
            t_fit_long = np.linspace(tfit[0], tpb, 1000)
            k = p[0]
            tprime = p[1]
            # calculate fit current
            Ifit_long = []
            for n in range(len(t_fit_long)):
                if t_fit_long[n]<=tlam:
                    Ifit_long.append(k/(np.sqrt(t_fit_long[n] - tprime)) + m*t_fit_long[n] + b)
                else:
                    Ifit_long.append(k/(np.sqrt(t_fit_long[n] - tprime)) - m*(t_fit_long[n] - tlam) - b)
            Ifit_long = np.array(Ifit_long)
            # calculate back current
            cotpb = Ifit_long[np.argmin(np.abs(t_fit_long - tpb))]
            Ipb = Ipb0 - cotpb
            Ipf = Ipf[0]
            Ipb = Ipb[0]
            ratio.append(np.abs(Ipb/Ipf))

            # plot if wanted
            if plot==True:
                axbig[i].plot(t, I, '-b')

                axbig[i].plot(lin_for, m*lin_for + b, '--k')
                axbig[i].plot([tpf, tpf], [linpf, Ipf0], 'o--g', label=r'$i_{pf}$ = %.4g µA'%(Ipf))
                axbig[i].plot(tfit, Ifit, '.k')

                axbig[i].plot([tpb, tpb], [cotpb, Ipb0], 'o--r', label=r'$i_{pr}$ = %.4g µA'%(Ipb))
                axbig[i].plot(t_fit_long, Ifit_long, '--k')
                axbig[i].axvline(x=tlam, color='k', linestyle='--', alpha=0.2, label=r'$t_{\lambda}$')

                if v>=1:
                    axbig[i].set_title(r'$v = %.3g$ V/s   $|i_{pr}/i_{pf}|$ = %.3g'%(v, np.abs(Ipb/Ipf)))
                else:
                    axbig[i].set_title(r'$v = %.3g$ mV/s   $|i_{pr}/i_{pf}|$ = %.3g'%(v*1000, np.abs(Ipb/Ipf)))
                axbig[i].plot(lin_for, m*lin_for + b, '--k')
                axbig[i].set_xlabel(r'$t$ / s')
                axbig[i].set_ylabel(r'$i$ / µA')
                axbig[i].legend(loc='upper left', fontsize=10)
        if plot==True:
            figbig.tight_layout()
            figbig.savefig('fits.pdf', transparent=True)

        # make figure for peak current ratio 
        fig, ax = plt.subplots(1,1,figsize=(5,3.5))
        ax.axhline(y=1, color='k', linestyle='--', alpha=0.2)
        for i in range(len(scan_rates)):
            ax.plot(scan_rates[i], ratio[i], 'o', markerfacecolor=self.colors[i], markeredgecolor='k')
        ax.set_ylabel(r'$|i_{pr}/i_{pf}|$')
        ax.set_xlabel(r'$v$ / V$\cdot$s$^{-1}$')
        ax.set_ylim([-0.1, 2.1])
        fig.tight_layout()
        fig.savefig('peak_current_ratio.svg', transparent=True)

    def randles_sevcik(self, scan_rates, E_range, conc=None, A=None):
        """methof to perform a Randles-Ševčík to obtain the diffusion coefficient"""
        Ip = []
        # file index and range
        for i in range(len(scan_rates)):
            E = self.x[i]
            I = self.y[i]
            if E[1]<E[5]:
                scan_direction = 'ox'
            # reductive scan
            else:
                scan_direction = 'red'
            # cut in desried voltage range
            I = I[(E>=E_range[0])&(E<=E_range[1])]
            E = E[(E>=E_range[0])&(E<=E_range[1])]
            if scan_direction=='ox':
                Ip.append(np.max(I))
            else:
                Ip.append(-1*np.min(I))
        fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
        p, pcov = np.polyfit(np.sqrt(scan_rates), Ip, 1, cov=True)
        ax.plot(np.sqrt(scan_rates), p[0]*np.sqrt(scan_rates) + p[1], '--k')
        for i in range(len(scan_rates)):
            ax.plot(np.sqrt(scan_rates[i]), Ip[i], 'o', markerfacecolor=self.colors[i], markeredgecolor='k')
        if scan_direction=='red':
            ax.set_ylabel(r'$|i_{pc}|$ / µA')
        else:
            ax.set_ylabel(r'$|i_{pa}|$ / µA')
        # conc in mM
        if conc:
            # convert slope to A s^{-1/2}
            slope = p[0]*10**-6
            err_slope = (pcov[0,0]**(0.5))*10**-6
            # number of e– transfered
            n = 1
            # constant in the equation C mol^-1 V^(-0.5)
            cont = 2.69e5
            # we need to convert the concentration from mM to mol/cm^3
            conc = conc*10**(-6)
            # calculate diffusion coefficient
            D = (slope/(cont * n**(1.5) * A * conc))**2
            print(D)
            Derr = 2*slope*err_slope/(cont * n**(1.5) * A * conc)**2
            print(r'D = (%.3g \pm %.3g) \times \, 10^{-5} \, \text{cm}^2 \cdot \text{s}^{-1}'%(D*10**5, Derr*10**5))
            ax.set_title(r'$D = (%.2g \pm %.2g) \times \, 10^{-5} \, \text{cm}^2 \cdot \text{s}^{-1}$'%(D*10**5, Derr*10**5))
        ax.set_xlabel(r'$v^{1/2} / \text{V}^{1/2} \, \text{s}^{-1/2}$')
        fig.tight_layout()
        fig.savefig('randles_sevcik.svg', transparent=True)

    def EC_pot_shift(self, scan_rates, E_range):
        """method to plot peak shift of irreversible CV and obtain slope"""
        Ep = []
        # file index and range
        for i in range(len(scan_rates)):
            E = self.x[i]
            I = self.y[i]
            if E[1]<E[5]:
                scan_direction = 'ox'
            # reductive scan
            else:
                scan_direction = 'red'
            # cut in desried voltage range
            I = I[(E>=E_range[0])&(E<=E_range[1])]
            E = E[(E>=E_range[0])&(E<=E_range[1])]
            if scan_direction=='ox':
                Ep.append(E[I==np.max(I)]*1000)
            else:
                Ep.append(E[I==np.min(I)]*1000)
        Ep = np.array(Ep)
        # make figure
        fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
        for i in range(len(scan_rates)):
            ax.plot(np.log10(scan_rates[i]), Ep[i], 'o', markerfacecolor=self.colors[i], markeredgecolor='k')
        p, cov = np.polyfit(np.log10(scan_rates), Ep, 1, cov=True)
        ax.plot(np.log10(scan_rates), p[0]*np.log10(scan_rates) + p[1], '--k')
        err = cov[0,0][0]**(0.5)
        slope = p[0][0]
        ax.set_title(r'$\partial E_\text{p} / \partial \log v = (%.3g \pm %.3g)\,\text{mV}$'%(slope, err))
        ax.set_xlabel(r'$\log(v / \text{V}^{-1} \cdot \text{s})$')
        ax.set_ylabel(r'$E_\text{p}$ / mV')
        fig.tight_layout()
        fig.savefig('EC_peak_shift.svg', transparent=True)
    
    def show(self):
        self.show_plot(self.ax)
        if self.US:
            self.ax.invert_xaxis()
        self.save_fig()
        plt.show()