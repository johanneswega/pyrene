from pyrene.standard.packages import *
from dataclasses import dataclass

@dataclass
class Plotter():
    """class to plot any type of experimental data"""

    ### define initialization arguments ###
    labels : list = None
    colors : list = None
    ylabel : str = None
    xlabel : str = None
    yticks : bool = True
    yscale : str = 'linear'
    xscale : str = 'linear'
    ylim : list = None
    xlim : list = None
    figsize : tuple = (8, 5)
    marker : list = None
    alphas : list = None
    tight_layout : bool = True
    zeroline : bool = True
    fill : list = None
    waterfall : float = 0.0
    outside : bool = False
    fontsize : int = 10
    savefig : str = None
    slave : bool = False
    steady_state : list = None

    def plot_data(self, master_ax=None):
        """creates figure and axis object"""

        # make own figure if master 
        if not master_ax:
            self.fig, self.ax = plt.subplots(1, 1, figsize=self.figsize)
            ax = self.ax
        # otherwise plot data in respective axis object of the experimental class
        else:
            ax = master_ax

        # create labels
        if not self.labels:
            self.labels = [r'%s'%(s[:s.find('.')]) for s in self.files]

        # create fill
        if not self.fill:
            self.fill = [False for _ in self.files]

        # create colors
        if not self.colors:
            from pyrene.standard.misc import rainbow
            self.colors = rainbow(self.files)

        # create markers 
        if not self.marker: 
            self.marker = ['-' for _ in self.files]

        # create markers 
        if not self.alphas: 
            self.alphas = [1 for _ in self.files]

        if self.zeroline:
            ax.axhline(y=0, color='k')

        # loop through files and plot
        for i in range(len(self.files)):
            if self.ma[i]:
                if self.wn:
                        self.x_ma[i] = 1e4/self.x_ma[i]
                if self.marker[i]:
                    ax.plot(self.x[i], self.y[i] - self.waterfall * i, self.marker[i], color=self.colors[i], alpha=0.5)
                    ax.plot(self.x_ma[i], self.y_ma[i] - self.waterfall * i, self.marker[i], label=self.labels[i], color=self.colors[i])
                if self.fill[i]:
                    ax.fill_between(self.x_ma[i], -self.waterfall * i, self.y_ma[i] - self.waterfall * i, alpha=0.1, color=self.colors[i])
            else:
                if self.marker[i]:
                    ax.plot(self.x[i], self.y[i] - self.waterfall * i, self.marker[i], label=self.labels[i], color=self.colors[i], alpha=self.alphas[i])
                if self.fill[i]:
                    ax.fill_between(self.x[i], -self.waterfall * i, self.y[i] - self.waterfall * i, alpha=0.1, color=self.colors[i])

        # plot steady state spectra
        if self.steady_state:
            self.plot_steady_state(ax)


    def plot_steady_state(self, ax):
        """method to plot exported steady state spectra on axis object"""
        for i in range(len(self.steady_state)):
            # read data
            data = np.loadtxt(self.steady_state[i][0], delimiter=',', skiprows=1)
            s_wl = data[:,0]
            s_wn = data[:,1]
            s = data[:,2]
            # multiply intensity by lambda^4 if emission
            if 'em' in self.steady_state[i][0]:
                s = s*s_wl**4
                print(self.steady_state[i][0], ' has been multiplied by lambda^4 to convert to SE.')
            if self.IR==False:
                s = s[(s_wl>self.steady_state[i][1][0]) & (s_wl<self.steady_state[i][1][1])]
                s_wn = s_wn[(s_wl>self.steady_state[i][1][0]) & (s_wl<self.steady_state[i][1][1])]
                s_wl = s_wl[(s_wl>self.steady_state[i][1][0]) & (s_wl<self.steady_state[i][1][1])]
            else:
                s = s[(s_wn>self.steady_state[i][1][0]) & (s_wn<self.steady_state[i][1][1])]
                s_wl = s_wl[(s_wn>self.steady_state[i][1][0]) & (s_wn<self.steady_state[i][1][1])]
                s_wn = s_wn[(s_wn>self.steady_state[i][1][0]) & (s_wn<self.steady_state[i][1][1])]                    
            s = self.steady_state[i][2]*s/np.nanmax(s)
            if not self.wn:
                ax.fill_between(s_wl, 0, s, color=self.steady_state[i][3], label=self.steady_state[i][4], alpha=0.1)
            else:
                ax.fill_between(s_wn, 0, s, color=self.steady_state[i][3], label=self.steady_state[i][4], alpha=0.1) 

    def show_plot(self, ax, title=None):
        """creates axis labels etc. and saves figure"""

        # create secondary axis if wn 
        if self.wn and not self.IR:
            ax2 = ax.secondary_xaxis(location='top', functions=(lambda x: 1e4/x, lambda x: 1e4/x))
            ax.set_xlabel(r'$\tilde{\nu} / 10^3\,\text{cm}^{-1}$')
            ax2.set_xlabel(r'$\lambda / \text{nm}$')
            ax.invert_xaxis()
        elif self.IR:
            ax.set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
        else:
            ax.set_xlabel(r'$\lambda / \text{nm}$')

        if self.zeroline:
            if self.waterfall:
                for i in range(len(self.files)):
                    ax.axhline(y=-self.waterfall*i, color='k', alpha=0.5)
            else:
                ax.axhline(y=0, color='k')

        if self.ylabel:
            ax.set_ylabel(self.ylabel)
        if self.xlabel:
            ax.set_xlabel(self.xlabel)  
        if self.xlim:
            ax.set_xlim(self.xlim)
        if self.ylim:
            ax.set_ylim(self.ylim)
        if self.yscale:
            ax.set_yscale(self.yscale)
        if title:
            ax.set_title(title)
        if self.xscale:
            ax.set_xscale(self.xscale)
        if not self.yticks: 
            ax.set_yticks([])  
        if self.outside:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=self.fontsize)
        else:
            ax.legend()

    def get_delay_labels(self, delays, title=False):
        """method to compute delay labels"""
        labels = []
        for i in range(len(delays)):
            if title:
                lab += r'$\Delta t = $'
            if self.experiment == 'femto':
                if np.abs(delays[i])<1:
                    lab = r'%.3g fs'%(delays[i]*1e3)
                elif 1<=np.abs(delays[i])<1e3:
                    lab = r'%.3g ps'%(delays[i])
                elif 1e3<=np.abs(delays[i])<1e6:
                    lab = r'%.3g ns'%(delays[i]/1e3)
                else:
                    lab = r'%.3g µs'%(delays[i]/1e6)
            if self.experiment == 'nano':
                if np.abs(delays[i])<1:
                    lab = r'%.3g ps'%(delays[i]*1000)
                elif 1<=np.abs(delays[i])<1000:
                    lab = r'%.3g ns'%(delays[i])
                else:
                    lab = r'%.3g µs'%(delays[i]/1000) 
            labels.append(lab)
        return labels
        
    def save_fig(self):
        if self.tight_layout:
            self.fig.tight_layout()
        if self.savefig: 
            self.fig.savefig(self.savefig, transparent=True)