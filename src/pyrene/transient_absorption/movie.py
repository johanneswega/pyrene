from pyrene.standard.packages import *
from pyrene.data_reader.read import DataReader
from pyrene.data_reader.export import DataExporter
from pyrene.plotter.plotter import Plotter
from dataclasses import dataclass

@dataclass
class Movie(DataReader, Plotter):
    """class to render movies of TA spectra"""
    # time in seconds of the movie 
    time : float = 10.0
    # wether to normalize each frame on its maximum 
    normall : bool = False
    # name of mp4 file
    movname : str = 'movie.mp4'
    # whether to plot previous temporal slices
    before : bool = False

    # automatically call read_data method from parent datareader after init
    def __post_init__(self):
        ### read data ###
        self.two_dim = True
        # to cut initial data frame correctly
        self.contour = True

        # get axis label
        if self.norm or self.devide or self.normall:
            self.ylabel = r'norm. $\Delta A$'
        else:
            self.ylabel = r'$\Delta A / 10^{-3}$'

        # initialize movenorm list
        if self.norm or self.devide:
            self.movnorm = [None for _ in self.files] 

        # read data for initial file
        self.read_data()
        # get delays to be plotted in movie
        self.delays = self.y[0]
        self.y_cuts = None
        if self.norm:
            self.devide = self.movnorm
        if not self.normall:
            self.norm = None
        else:
            self.norm = [True for _ in self.files]
        self.norm_at = None
        self.contour = False
        self.slicing = None

        # generate figure
        self.fig, self.ax = plt.subplots(1, 1, figsize=self.figsize)

    def animate(self, i):
        # clear axis 
        self.ax.clear()
        # get title 
        time_title = self.get_delay_labels([self.delays[i]], title=True)
        print("dt = %s  |   %i/%i"%(time_title[0][14:], i+1, len(self.delays)))
        self.alphas = [1.0 for _  in self.files]
        ma = self.ma
        lab = self.labels
        self.delay = [self.delays[i] for _ in self.files]
        self.read_data()
        self.plot_data(master_ax=self.ax)
        
        # plot previous time slices if wanted
        if self.before:
            self.no_labels = True
            self.ma = [False for _ in self.files]
            for j in range(len(self.delays[self.delays<self.delays[i]])):
                self.alphas = [0.05 for _  in self.files]
                self.delay = [self.delays[j] for _ in self.files]
                self.read_data()
                self.plot_data(master_ax=self.ax)
            self.no_labels = False
            self.labels = lab
            self.ma = ma

        self.show_plot(self.ax, title=time_title[0])
        self.fig.tight_layout()

    def render(self):
        anim = FuncAnimation(self.fig, func=self.animate, frames=len(self.delays), interval=1)
        Writer = writers['ffmpeg']
        writer = Writer(fps=round(len(self.delays)/self.time), metadata={'artist': 'Me'}, bitrate=2500)
        anim.save(self.movname, writer)
        plt.show()