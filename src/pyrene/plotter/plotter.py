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

    def plot_data(self):
        """creates figure and axis object"""
        self.fig, self.ax = plt.subplots(1, 1, figsize=self.figsize)

        # create markers 
        if not self.marker: 
            self.marker = ['-' for _ in self.files]

        # create markers 
        if not self.alphas: 
            self.alphas = [1 for _ in self.files]

        if self.zeroline:
            self.ax.axhline(y=0, color='k')

        # loop through files and plot
        for i in range(len(self.files)):
            if self.ma[i]:
                if self.wn:
                        self.x_ma[i] = 1e4/self.x_ma[i]
                self.ax.plot(self.x[i], self.y[i] - self.waterfall * i, self.marker[i], color=self.colors[i], alpha=0.5)
                self.ax.plot(self.x_ma[i], self.y_ma[i] - self.waterfall * i, self.marker[i], label=self.labels[i], color=self.colors[i])
            else:
                self.ax.plot(self.x[i], self.y[i] - self.waterfall * i, self.marker[i], label=self.labels[i], color=self.colors[i], alpha=self.alphas[i])


    def show_plot(self):
        """creates axis labels etc. and saves figure"""

        # create secondary axis if wn 
        if self.wn and not self.IR:
            self.ax2 = self.ax.secondary_xaxis(location='top', functions=(lambda x: 1e4/x, lambda x: 1e4/x))
            self.ax.set_xlabel(r'$\tilde{\nu} / 10^3\,\text{cm}^{-1}$')
            self.ax2.set_xlabel(r'$\lambda / \text{nm}$')
            self.ax.invert_xaxis()
        elif self.IR:
            self.ax.set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
        else:
            self.ax.set_xlabel(r'$\lambda / \text{nm}$')

        if self.zeroline:
            if self.waterfall:
                for i in range(len(self.files)):
                    self.ax.axhline(y=-self.waterfall*i, color='k', alpha=0.5)
            else:
                self.ax.axhline(y=0, color='k')

        if self.ylabel:
            self.ax.set_ylabel(self.ylabel)
        if self.xlabel:
            self.ax.set_xlabel(self.xlabel)  
        if not self.yticks: 
            self.ax.set_yticks([])  
        if self.outside:
            self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=self.fontsize)
        else:
            self.ax.legend()
        if self.tight_layout:
            self.fig.tight_layout()
        if self.savefig: 
            self.fig.savefig(self.savefig, transparent=True)