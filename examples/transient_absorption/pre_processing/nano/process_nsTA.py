from pyrene.transient_absorption.pre_processing import pixel_to_lambda, load_nsTA_file
from pyrene.transient_absorption import Contour, Spectra, Kinetics
from pyrene.standard.misc import save, load
from pyrene.standard.packages import *
import shutil

### set parameters ###
# pixel limits for pixel to lambda conversion 
limits = [0, 500]
# fit parameters (scale, shift)
p0 = [0.75, 357]
# wether to do automatic fit for pixel to lambda or pick points by hand
auto = True
# scale for contour plots
scale = [-20, 20]
# time delay from which < to calculate background 
tbg = -1.0

### analyze data ###
# make results folder 
if os.path.exists('results'):
    shutil.rmtree('results')
os.mkdir('results')

# pixel to lambda
if not 'wl.txt' in os.listdir():
    pixel_to_lambda(lim=limits, p0=p0, auto=auto, nano=True)

# load data
t, dA, nsamples = load_nsTA_file(file = [i for i in os.listdir() if 'SAMPLE' in i][0])
wl = np.loadtxt('wl.txt') 
save('results/dA.npy', t, wl, dA)

# determine time zero
if not 't0.txt' in os.listdir():
    from pyrene.transient_absorption.pre_processing.pixel_to_lambda import onclick
    t0_near = []
    c = Contour(files=['results/dA.npy'], titles=['select point near time zero!'], extend=['both'], scale=[scale], wn=False, nlevels=[31],
                yscale='symlog', lines=[True], zeroline=False, experiment='nano')
    cid = c.fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, [], t0_near, c.ax, c.fig))
    c.show()
    # zoom on time zero
    c = Contour(files=['results/dA.npy'], titles=['select precise time zero!'], extend=['both'], scale=[scale], wn=False, nlevels=[31],
                yscale='linear', lines=[True], ylim=[t0_near[0]-2, t0_near[0]+2], zeroline=False, experiment='nano')
    t0 = []
    cid = c.fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, [], t0, c.ax, c.fig))
    c.show()
    np.savetxt('t0.txt', t0)

# load time zero and correct
t0 = np.loadtxt('t0.txt')
t -= t0
save('results/dA.npy', t, wl, dA)
c = Contour(files=['results/dA.npy'], titles=['raw data'], extend=['both'], scale=[scale], wn=False, nlevels=[31],
            yscale='symlog', ylim=[-10, 500e3], lines=[True], zeroline=True, experiment='nano', savefig='results/raw_data.png')
c.show()

# background correct
from pyrene.standard.misc import find_index
bg = np.mean(dA[:find_index(t, tbg), :], axis=0)
from pyrene.standard.figures import figure_spectum
fig, ax = figure_spectum(wn=False, figsize=(5, 4))
ax.plot(wl, bg, '-b')
ax.set_title(rf'averaged background for $\Delta t < {tbg}$ ps')
ax.set_ylabel(r'$\Delta A / 10^{-3}$')
fig.tight_layout()
fig.savefig('results/background.png', transparent=True)
plt.show()

# subtract background and plot 
for i in range(len(t)):
    dA[i, :] -= bg
save('results/dA.npy', t, wl, dA)
c = Contour(files=['results/dA.npy'], titles=['background subtracted'], extend=['both'], scale=[scale], wn=False, nlevels=[31],
            yscale='symlog', ylim=[-10, 500e3], lines=[True], zeroline=True, experiment='nano', savefig='results/background_subtracted.png')
c.show()

# plot counting statistics
fig,ax = plt.subplots(ncols=1,nrows=1, figsize=(8,5))
ax.plot(t+t0, nsamples, '.r')
ax.set_xlabel(r'$\Delta t$ / ns')
ax.set_ylabel(r'samples')
ax.set_yscale('log')
ax.set_xscale('symlog')
ax.set_title('Counting Statistics')
fig.savefig('results/counting_statitics.png')
fig.tight_layout()
plt.show()

# plot overview plot
s = Spectra(files=['results/dA.npy'], overview=True, outside=True, savefig='results/overview.png', experiment='nano')
s.show()

# get and save error matrix
tx, RMS, nsamples = load_nsTA_file(file = [i for i in os.listdir() if 'SAMPLE' in i][0], error=True)
Err = RMS / np.sqrt(nsamples)[:, None]  
np.save('results/Err.npy', np.array([t, wl, Err], dtype=object))
c = Contour(files=['results/Err.npy'], titles=['background subtracted'], extend=['both'], scale=[(-4, 4)], wn=False, nlevels=[31],
            yscale='symlog', ylim=[-10, 500e3], lines=[True], zeroline=True, experiment='nano', savefig='results/background_subtracted.png')
c.show()

# convert to txt
from pyrene.standard.misc import convert_to_txt
convert_to_txt('results/dA.npy', experiment='nano', t_cuts=None, wl_cuts=None)
shutil.move("wavelength.txt", "results/wavelength.txt")
shutil.move("time.txt", "results/time.txt")
shutil.move("TA.txt", "results/TA.txt")