from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from pyrene.transient_absorption import Contour, Spectra, Kinetics
from pyrene.standard.misc import save
from dataclasses import dataclass
import shutil

@dataclass
class GlobalAnalysis(DataReader, Plotter, DataExporter):
    """class to preform a global kinetic analyis of a TA data set"""

    p0 : list = None
    IRF : list = None
    error : list = None
    model : str = 'sequential'
    # list of fit parameters to fix, i.e. if you want to fix tau1 for file0 then fix = [(False, True, False, ...)]
    fix : list = None
    C0 : list = None
    K : list = None

    def fit(self):
        """class to perform global fit"""
        self.error = self.set_standard_value(self.error, False)
        self.IRF = self.set_standard_value(self.IRF, False)
        self.p0 = self.set_standard_value(self.p0, False)
        self.fix = self.set_standard_value(self.fix, False)
        self.C0 = self.set_standard_value(self.C0, False)
        self.K = self.set_standard_value(self.K, False)
        self.scale =  self.set_standard_value(self.scale, default=[-20, 20])
        self.scatter =  self.set_standard_value(self.scatter, default=None)
        self.wavelength =  self.set_standard_value(self.wavelength, default=np.arange(350, 750, 20))

        for i in range(len(self.files)):
            self.current_file = i
            # show data matrix 
            c = Contour(files=[self.files[i]], x_cuts=self.x_cuts, y_cuts=self.y_cuts, titles=['experimental data'],  experiment=self.experiment,
                        scale=self.scale, scatter=self.scatter, extend=['both'], yscale=self.yscale, lines=[True], zeroline=False, IR=self.IR)
            self.t = c.y[0]
            if not self.IR:
                self.wl = 1e4/c.x[0]
            else:
                self.wl = c.x[0]
            self.dA = c.z[0]
            # remove nans 
            mask = ~np.isnan(self.dA).any(axis=0)
            self.wl = self.wl[mask]
            self.dA = self.dA[:, mask]
            c.show()

            # get error matrix if wanted
            if self.error[i]:
                c = Contour(files=[self.error[i]], x_cuts=self.x_cuts, y_cuts=self.y_cuts, titles=['error matrix'], experiment=self.experiment,
                            scale=[(-0.5, 0.5)], scatter=self.scatter, extend=['both'], yscale=self.yscale, lines=[True], zeroline=False, IR=self.IR)
                self.sigma = c.z[0]
                # remove nans 
                mask = ~np.isnan(self.sigma).any(axis=0)
                self.sigma = self.sigma[:, mask]
                c.show()                

            # perform SVD
            S = np.linalg.svd(self.dA)
            sigma = S[1][:15]
            N = np.linspace(1,15,15)
            fig, ax = plt.subplots(1,1,figsize=(3.5, 3.5))
            ax.plot(N, np.log10(sigma), '.b')
            ax.set_ylabel(r'$\log{(\sigma_i)}$')
            ax.set_xlabel(r'$i$')
            ax.set_title("singular values")
            fig.tight_layout()
            # make results folder 
            self.name = self.files[i][:self.files[i].find('.')]
            if os.path.exists(self.name + '_global_analysis'):
                shutil.rmtree(self.name + '_global_analysis')
            os.mkdir(self.name + '_global_analysis')
            fig.savefig(self.name + '_global_analysis/SVD.png', transparent=True)
            plt.show()

            # perform fit
            b = self.get_bounds(self.p0[i], self.IRF[i])
            self.fit_result = least_squares(self.least_squares_fit, self.p0[i], bounds=b)
            self.p = self.fit_result.x

            # save fitted array and residuals
            self.C = self.calculate_C(self.p)
            self.Sim = np.dot(self.C, np.dot(np.linalg.pinv(self.C), self.dA))
            if not self.error[self.current_file]:
                self.Res = (self.Sim - self.dA)
                chi2 = np.sum((self.Res)**2)
                lab = r' ($\chi^2 = %.3g$)'%chi2
            else:
                self.Res = ((self.Sim - self.dA)/self.sigma)
                chi2 = np.sum((self.Res)**2) / (self.Res.size - len(self.p0[i]))
                lab = r' ($\chi_\nu^2 = %.3g$)'%chi2
            save(self.name + '_global_analysis/' + 'Sim.npy', self.t, self.wl, self.Sim)
            save(self.name + '_global_analysis/' + 'Res.npy', self.t, self.wl, self.Res)

            # print errors and get lables
            title = self.print_errors()

            # plot fit (self.scale[i][0]/10, self.scale[i][1]/10)
            c = Contour(files=[self.files[i], self.name + '_global_analysis/' + 'Sim.npy', self.name + '_global_analysis/' + 'Res.npy'],
                        figsize=[17, 4.5], x_cuts=self.x_cuts, y_cuts=self.y_cuts, 
                        titles=['experimental data', 'global fit', 'residuals' + lab],
                        scale=[self.scale[i], self.scale[i], self.scale[i]], 
                        cmap=['RdBu_r', 'RdBu_r', 'seismic'], experiment=self.experiment,
                        scatter=[self.scatter[i]], extend=['both'], yscale=self.yscale, IR=self.IR,
                        lines=[True, True, False], zeroline=False, savefig=self.name + '_global_analysis/' + '/fit.svg')
            if self.scatter[self.current_file]!=None:
                import matplotlib.patches as patches
                rect = [(1/self.scatter[self.current_file][1])*10**4, (1/self.scatter[self.current_file][0])*10**4]
                c.ax[2].add_patch(patches.Rectangle((rect[0], c.ax[2].get_ylim()[0]), rect[1]-rect[0], c.ax[2].get_ylim()[1]-c.ax[2].get_ylim()[0], facecolor='white', zorder=2))
                c.ax[1].add_patch(patches.Rectangle((rect[0], c.ax[1].get_ylim()[0]), rect[1]-rect[0], c.ax[1].get_ylim()[1]-c.ax[1].get_ylim()[0], facecolor='white', zorder=2))
            c.show()

            # plot fitted traces at selected wavelengths
            self.plot_kinetics()
            # plot species spectra
            self.plot_EADS(title=title)

    def get_bounds(self, p0, IRF):
        """methods to estimate bounds of the fit"""
        b = []
        if IRF: 
            # lower bound for tau
            b0 = [0 for _ in range(len(p0)-2)]
            # upper bound for tau
            b1 = [np.inf for _ in range(len(p0)-2)]
            # add bounds for t0 and fwhm
            b0.append(-10)
            b0.append(-10)
            b1.append(10)
            b1.append(10)
        else:
            b0 = [0 for _ in range(len(p0))]
            b1 = [np.inf for _ in range(len(p0))]
        b.append(b0)
        b.append(b1)
        return b
    
    def least_squares_fit(self, p):
        """method to perform least squares fit"""
        C = self.calculate_C(p)
        self.Sim = np.dot(C, np.dot(np.linalg.pinv(C), self.dA))
        if not self.error[self.current_file]:
            return (self.Sim - self.dA).flatten()
        else:
            return ((self.Sim - self.dA)/self.sigma).flatten()

    def calculate_C(self, p):
        """method to construct C-matrix"""
        # get parameters
        if not self.IRF[self.current_file]:
            k = 1/p
        else:
            k = 1/p[:-2]
            t0 = p[-2]
            fwhm = p[-1]
            sigma = fwhm/(2*(2*np.log(2))**0.5)

        # add custom constraints between rate constants 
        # k = np.array([k[0], (0.625/(1 - 0.625))*k[0], k[1]]) 

        # fix parameters if desired
        if self.fix[self.current_file]:
            for i in range(len(self.fix[self.current_file])):
                if self.fix[self.current_file][i]:
                    k[i] = 1/self.p0[self.current_file][i]

        # get K-matrix for parallel model (DADS)
        if self.model=='parallel':
            K = np.diag(-k)
            # initial values
            C0 = [1 for _ in range(len(K[0]))]

        # get K-matrix for sequential model (EADS)
        if self.model=='sequential':
            K = np.zeros((len(k), len(k)))
            K[0,0] = -1*k[0]
            for i in range(1, len(K)):
                K[i,i-1] = k[i-1]
                K[i,i] = -1*k[i]
            # initial values
            C0 = [1 if i == 0 else 0 for i in range(len(K[0]))]

        if self.model=='target':
            # build the target matrix from the string matrix 
            K = self.build_target_K(k)

        if self.C0[self.current_file]:
            C0 = np.array(self.C0[self.current_file])

        # Calculate eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eig(K)

        # Create the diagonal matrix of exponential eigenvalues
        # expm_Lambda = np.diag(np.exp(eigenvalues))

        # Matrices of eigenvectors
        U = eigenvectors

        # Calculate the inverse of U
        U_inv = np.linalg.inv(U)

        # Calculate matrix exponential
        # expm_K = np.dot(U, np.dot(expm_Lambda, U_inv))

        # calculate alpha
        alpha = np.dot(U_inv, C0)

        # calculate theta
        if not self.IRF[self.current_file]:
            theta = np.zeros((len(eigenvalues), len(self.t)))
            for i in range(len(eigenvalues)):
                theta[i,:] = alpha[i]*np.exp(eigenvalues[i]*self.t)
        else:
            from scipy.special import erf
            psi = np.zeros((len(eigenvalues), len(self.t)))
            theta = np.zeros((len(eigenvalues), len(self.t)))
            for i in range(len(eigenvalues)):
                psi[i,:] = 0.5*np.exp(eigenvalues[i]*(self.t - t0 + eigenvalues[i]*sigma**2/2))*(1 + erf((self.t-t0+eigenvalues[i]*sigma**2)/(sigma*2**0.5)))
            for i in range(len(eigenvalues)):
                theta[i,:] = alpha[i]*psi[i,:]

        # calculate C
        C = np.transpose(np.dot(U, theta))
        return C

    # method to print errors
    def print_errors(self):
        print("")
        from pyrene.standard.misc import round_to_significant_digit
        dic = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I'}
        # Calculate the covariance matrix
        covariance_matrix = np.linalg.pinv(np.dot(self.fit_result.jac.T, self.fit_result.jac))
        # Calculate the parameter errors
        p_err = np.sqrt(np.diagonal(covariance_matrix))

        if self.IRF[self.current_file]: 
            tau_labels = self.get_delay_labels(self.p[:-2])
            err_tau_labels = self.get_delay_labels(p_err[:-2])
        else:
            tau_labels = self.get_delay_labels(self.p)
            err_tau_labels = self.get_delay_labels(p_err)
        
        label = ''
        for i in range(len(tau_labels)):
            print(f"tau_{dic[i]} = ", tau_labels[i], " +- ", err_tau_labels[i])
            pre = round_to_significant_digit(float(tau_labels[i][:tau_labels[i].find(' ')]))
            suff = tau_labels[i][tau_labels[i].find(' '):]
            if pre>=10:
                pre = str(pre)
                pre = pre[:pre.find('.')]
            else:
                pre = str(pre)
            if self.model=='sequential':
                if p_err[i]<1e6:
                    label += r"%s $\xrightarrow{\text{%s}}$ "%(dic[i], pre + suff)
                else:
                    label += r"%s $\xrightarrow{\infty}$ "%(dic[i])
            else:
                if p_err[i]<1e3:
                    label += r" $\tau_{\text{%s}} = $ "%dic[i] + pre + suff + ','
                else:
                    label += r" $\tau_{\text{%s}} = $ "%dic[i] + r'$\infty$'

        label = label[:-1]
        # save lifetimes + error
        if self.experiment=='femto':
            np.savetxt(self.name + '_global_analysis/' + 'lifetimes.txt', np.column_stack((self.p, p_err)), delimiter=',', header='species lifetime / ps, error /ps')
        else:
            np.savetxt(self.name + '_global_analysis/' + 'lifetimes.txt', np.column_stack((self.p, p_err)), delimiter=',', header='species lifetime / ps, error /ns')
        print("")
        return label

    def plot_kinetics(self):
        """method to plot kinetic traces and residuals of the fit at a single wavelength level"""
        from pyrene.standard.misc import find_index, rainbow
        self.wavelength[self.current_file] = np.array(self.wavelength[self.current_file])
        # make rainbow colormap for all chosen kinectics
        colors = rainbow(self.wavelength[self.current_file])
        colors_darker = colors.copy()
        colors_darker = [(0.4*r, 0.4*g, 0.4*b, a) for r, g, b, a in colors]
        # generate figure
        fig, ax = plt.subplots(2, 1, figsize=(9, 5), gridspec_kw={'height_ratios':[1,3]}, sharex=True)
        # go through chosen wavenumbers and plot raw self.dAta, fit as well as the self.Residuals
        for i in range(len(self.wavelength[self.current_file])):
            # find wavenumber index closest to the wavenumber you picked 
            index = find_index(self.wl, self.wavelength[self.current_file][i])
            # plot experimental data
            if self.IR==False:
                ax[1].plot(self.t, self.dA[:, index], 'o', alpha=0.3, markersize=5, color=colors[i], label=r'%.3g nm'%(self.wavelength[self.current_file][i]))
            else:
                ax[1].plot(self.t, self.dA[:, index], 'o', alpha=0.3, markersize=5, color=colors[i], label=r'%i cm$^{-1}$'%(round(self.wavelength[self.current_file][i])))
            # plot fit
            ax[1].plot(self.t, self.Sim[:, index], '-', color=colors_darker[i], linewidth=2.5)
            # plot residuals
            ax[0].plot(self.t, self.Res[:, index], '-', color=colors[i])
        # set x-lables
        if self.experiment=='nano':
            ax[1].set_xlabel(r'$\Delta t / \text{ns}$')
            ax[1].set_xscale('log')
        else:
            ax[1].set_xlabel(r'$\Delta t / \text{ps}$')
        
        # stylistic stuff
        ax[1].set_xscale(self.yscale)
        ax[0].axhline(y=0, color='k')
        if not self.error[self.current_file]:
            ax[0].set_ylabel(r'Res.')
        else:
            ax[0].set_ylabel(r'Res. / $\sigma$')
        ymin, ymax = ax[0].get_ylim()
        if abs(ymin)>ymax:
            m = abs(ymin)
        else:
            m = ymax
        ax[0].set_ylim([-m, m])
        ax[1].set_ylabel(r'$\Delta{A} / 10^{-3}$')
        ax[1].axhline(y=0, color='k')
        ax[1].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)
        fig.tight_layout()  
        fig.savefig(self.name + '_global_analysis/' + 'kinetics.svg', transparent=True)
        plt.show()

    def plot_EADS(self, title):
        """method to plot S and C matrices"""
        fig,ax = plt.subplots(ncols=2,nrows=1, figsize=(10, 4)) 
        dic = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I'}
        col = ['r', 'b', 'g', 'orange', 'purple', 'k']
        EAS = np.dot(np.linalg.pinv(self.C), self.dA)
        if not self.IR:
            l = 1e4/self.wl
        else:
            l = self.wl

        # time evolution
        for i in range(len(self.C[0,:])):
            ax[0].plot(self.t, self.C[:,i], '-', color=col[i], label='%s'%(dic[i]))
        ax[0].legend()
        ax[0].set_ylabel('rel. concentration')
        if self.experiment=='nano':
            ax[0].set_xlabel(r'$\Delta t / \text{ns}$')
            if np.min(self.t)>0:
                ax[0].set_xscale('log')
            else:
                ax[0].set_xscale('symlog')
        else:
            ax[0].set_xlabel(r'$\Delta t / \text{ps}$')
            if np.min(self.t)>0:
                ax[0].set_xscale('log')
            else:
                ax[0].set_xscale('symlog')
        if self.model=='sequential':
            ax[1].set_ylabel('EADS / mOD')
        if self.model=='parallel':
            ax[1].set_ylabel('DADS / mOD')
        if self.model=='target':
            ax[1].set_ylabel('SADS / mOD')
        ax[0].set_title(title)
        ax[0].set_ylim([-0.1,1.1])

        # spectra
        for i in range(len(self.C[0,:])):
            ax[1].plot(l, EAS[i,:], '-', color=col[i], zorder=1)
            if not self.IR:
                np.savetxt(self.name + '_global_analysis/' + '%s.txt'%(dic[i]), np.column_stack([l, (1/l)*10**4, EAS[i,:]]), delimiter=',',
                            header='wavelength / nm, wavenumber / 10^3 cm-1, EADS / mOD')
            else:
                np.savetxt(self.name + '_global_analysis/' + '%s.txt'%(dic[i]), np.column_stack([(1/l)*10**4, l, EAS[i,:]]), delimiter=',',
                            header='wavelength / µm, wavenumber / cm-1, EADS / mOD')                
        ax[1].axhline(y=0, color='k')
        if self.IR==False:
            ax[1].invert_xaxis()
            ax[1].set_xlabel(r'$\tilde{\nu} / 10^{3} \, \text{cm}^{-1}$')
            ax2 = ax[1].secondary_xaxis("top", functions=(lambda x: (1/x)*10**+4,lambda x: (1/x)*10**-4))
            ax2.set_xlabel(r'$\lambda / \text{nm}$')
        else:
            ax[1].set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
        fig.tight_layout()
        if self.scatter[self.current_file]!=None:
            import matplotlib.patches as patches
            rect = [(1/self.scatter[self.current_file][1])*10**4, (1/self.scatter[self.current_file][0])*10**4]
            ax[1].add_patch(patches.Rectangle((rect[0], ax[1].get_ylim()[0]), rect[1]-rect[0], ax[1].get_ylim()[1]-ax[1].get_ylim()[0], facecolor='white'))
        fig.savefig(self.name + '_global_analysis/' + 'spectra.svg', transparent=True)
        plt.show()

    def build_target_K(self, k):
        """helper function to build target matrix from string"""
        import re
        Kt = np.zeros(self.K[self.current_file].shape, dtype=float)

        for i in range(self.K[self.current_file].shape[0]):
            for j in range(self.K[self.current_file].shape[1]):

                entry = self.K[self.current_file][i, j].replace(" ", "")

                if entry == "0":
                    continue

                tokens = re.findall(r'[+-]?k\d+', entry)

                val = 0.0
                for tok in tokens:
                    sign = -1 if tok.startswith('-') else 1
                    idx = int(tok.lstrip('+-')[1:])
                    val += sign * k[idx]

                Kt[i, j] = val

        return Kt