from pyrene.transient_absorption.pre_processing import pixel_to_lambda, load_fsTA_file
from pyrene.transient_absorption import Contour, Spectra
from pyrene.standard.misc import save, load
from pyrene.standard.packages import *
import shutil

### set parameters ###
# pixel limits for pixel to lambda conversion 
limits = [0, 500]
# fit parameters (scale, shift)
p0 = [0.947, 274]
# wether to do automatic fit for pixel to lambda or pick points by hand
auto = True
# scale for contour plots
scale = [-20, 20]
# wether to exclude any scan 
exclude = []
# time delay from which < to calculate background 
tbg = -1.0
# set true if already averaged and background subtracted
already_averaged = True
# tlim for chrip contour plot
tlim_chirp = [-1.0, 1.5]

### analyze data ##
if not already_averaged:
    # make results folder 
    if os.path.exists('results'):
        shutil.rmtree('results')
    os.mkdir('results')

    # pixel to lambda
    if not 'wl.txt' in os.listdir():
        pixel_to_lambda(lim=limits, p0=p0, auto=auto)

    # load data
    t, df = load_fsTA_file(file = [i for i in os.listdir() if 'SAMPLE' in i][0])
    wl = np.loadtxt('wl.txt') 

    # save the individual scans as .npy files and plot
    os.mkdir('results/scans')
    scan_files = []
    for i in range(len(df)):
        save('results/scans/scan' + f"{i}" + '.npy', t, wl, df[i])
        scan_files.append('results/scans/scan' + f"{i}" + '.npy')
    # split scan files in chuncks of three and plot contours
    chunks = [scan_files[i:i+3] for i in range(0, len(df), 3)]
    titles = ['scan %i'%i for i in range(len(df))]
    titles = [titles[i:i+3] for i in range(0, len(titles), 3)]
    for i in range(len(chunks)):
        c = Contour(files=chunks[i], titles=titles[i], extend=['both'], scale=[scale], wn=False, figsize=(len(chunks[i])*6, 4.5), nlevels=[31],
                    yscale='symlog', lines=[True], savefig=f'results/raw_scans_{i}.png', zeroline=False)
        c.show()
    # exclude any scans if wanted
    df = [df[i] for i in range(len(df)) if not i in exclude]

    # average scans and plot
    dA = np.mean(df, axis=0)
    save('results/scans/average.npy', t, wl, dA)
    c = Contour(files=['results/scans/average.npy'], titles=['average'], extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31],
                yscale='symlog', lines=[True], savefig=f'results/average.png', zeroline=False)
    c.show()

    # calculate background and plot
    from pyrene.standard.misc import find_index
    from pyrene.standard.figures import figure_spectum
    bg = np.mean(dA[:find_index(t, tbg), :], axis=0)
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
    save('results/scans/average_minus_bg.npy', t, wl, dA)
    c = Contour(files=['results/scans/average_minus_bg.npy'], titles=['average background subtracted'], 
                extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31],
                yscale='symlog', lines=[True], savefig=f'results/average_minus.png', zeroline=False)
    c.show()
else:
    t, wl, dA = load('results/scans/average_minus_bg.npy')
    from pyrene.fitting.fit_functions import koboyashi
    # do manual chirp correction if no OKE file present
    manuel = False
    if not 'ok.txt' in os.listdir():
        c = Contour(files=['results/scans/average_minus_bg.npy'], titles=['select points for chirp correction'], 
                    extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=tlim_chirp,
                    yscale='linear', lines=[True], zeroline=False)
        from pyrene.transient_absorption.pre_processing.pixel_to_lambda import onclick
        x = []
        y = []
        cid = c.fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, x, y, c.ax, c.fig))
        c.show()
        x = np.array(x)
        y = np.array(y)
        p, pcov = curve_fit(koboyashi, x, y)
        manuel = True
        # write fit_params to file
        with open('ok.txt', 'w') as file:
            # Write each value in the list as a separate line in the file
            for item in p:
                file.write(str(item) + '\n') 

    p = np.loadtxt('ok.txt')  
    # plot chirp fit 
    c = Contour(files=['results/scans/average_minus_bg.npy'], titles=['chirp fit'], 
                extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=tlim_chirp,
                yscale='linear', lines=[True], zeroline=False)
    c.ax.plot(wl, koboyashi(wl, *p), '--k')
    if manuel:
        c.ax.plot(x, y, 'wo', markeredgecolor='k')
    c.show()

    # correct for chirp 
    for i in range(len(wl)):
        t0 = koboyashi(wl[i], *p)
        dA[:,i] = np.interp(t, t-t0, dA[:,i])   
    save('results/dA.npy', t, wl, dA)
    # plot chirp corrected 
    c = Contour(files=['results/dA.npy'], titles=['chirp corrected data'], 
                extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=[-0.5, 1.5],
                yscale='linear', lines=[True], zeroline=True)
    c.show()
    # plot on whole scale 
    c = Contour(files=['results/dA.npy'], titles=['final TA spectrum'], 
                extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=[-1.0, 1800],
                yscale='symlog', lines=[True], zeroline=True)
    c.show()
    # plot overview plot
    s = Spectra(files=['results/dA.npy'], overview=True, outside=True)
    s.show()