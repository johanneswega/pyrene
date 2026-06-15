from pyrene.standard.packages import *
from dataclasses import dataclass

@dataclass
class Plotter():
    """class to plot any type of experimental data"""

    ### define initialization arguments ###
    labels : list = None
    colors : list = None
    ylabel : str = ''
    xlabel : str = ''
    figsize : tuple = (8, 5)
    marker : list = None

    def plot_data(self):
        """creates figure and axis object"""
        self.fig, self.ax = plt.subplots(1, 1, figsize=self.figsize)

        # create markers 
        if not self.marker: 
            self.marker = ['-' for _ in self.files]

        # loop through files and plot
        for i in range(len(self.files)):
            self.ax.plot(self.x[i], self.y[i], self.marker[i], label=self.labels[i], color=self.colors[i])
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_xlabel(self.xlabel)    