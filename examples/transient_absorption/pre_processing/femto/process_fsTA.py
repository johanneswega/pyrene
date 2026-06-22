from pyrene.transient_absorption.pre_processing import pixel_to_lambda, load_fsTA_file
from pyrene.transient_absorption import Contour, Spectra, Kinetics
from pyrene.standard.misc import save, load
from pyrene.standard.packages import *
import shutil

### set parameters ###
# pixel limits for pixel to lambda conversion 
limits = [0, 500]
# fit parameters (scale, shift)
p0 = [0.947, 274]
# wether to do automatic fit for pixel to lambda or pick points by hand
auto = False
# scale for contour plots
scale = [-20, 20]
# wether to exclude any scan 
exclude = []
# time delay from which < to calculate background 
tbg = -1.0
# set true if already averaged and background subtracted
already_averaged = False
# tlim for chrip contour plot
tlim_chirp = [-0.5, 2.0]
# scale for OKE plot
scale_OKE = [0, 100]
# scatter
scatter = [560, 630]
# wavelengths to analyze in OKE measurement
wl_look = np.arange(350, 700, 10)
wl_look = np.concatenate([wl_look[wl_look<scatter[0]], wl_look[wl_look>scatter[1]]])
# whether to save data for individual scans
save_scans = False

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
    save('results/average_minus_bg.npy', t, wl, dA)
    c = Contour(files=['results/average_minus_bg.npy'], titles=['average background subtracted'], 
                extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31],
                yscale='symlog', lines=[True], savefig=f'results/average_minus.png', zeroline=False)
    c.show()

    # extract error matrix
    t, df = load_fsTA_file(file = [i for i in os.listdir() if 'SAMPLE' in i][0], error=True)
    Err = np.mean(df, axis=0)
    save('results/Error.npy', t, wl, Err)

t, wl, dA = load('results/average_minus_bg.npy')

### chirp correction ###
from pyrene.fitting.fit_functions import koboyashi
if not 'ok.txt' in os.listdir():
    # do manual chirp correction if no OKE file present
    if not any('OKE' in f for f in os.listdir()):
        c = Contour(files=['results/average_minus_bg.npy'], titles=['select points for chirp correction'], 
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
        # write fit_params to file
        with open('ok.txt', 'w') as file:
            # Write each value in the list as a separate line in the file
            for item in p:
                file.write(str(item) + '\n') 
        # plot chirp fit 
        c = Contour(files=['results/average_minus_bg.npy'], titles=['chirp fit'], 
                    extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=tlim_chirp,
                    yscale='linear', lines=[True], zeroline=False, savefig='results/chirp_fit.png')
        c.ax.plot(wl, koboyashi(wl, *p), '--k')
        c.ax.plot(x, y, 'wo', markeredgecolor='k')
        c.show()
    else:
        # do chirp correction via OKE
        # load OKE data
        tOKE, df = load_fsTA_file(file = [i for i in os.listdir() if 'OKE' in i][0])
        # average / background subtract etc.
        OKE = np.mean(df, axis=0)
        from pyrene.standard.misc import find_index
        bg = np.mean(OKE[:find_index(tOKE, tbg), :], axis=0)
        OKE -= bg
        OKE = np.abs(OKE)
        # save OKE file
        if os.path.exists('results/OKE'):
            shutil.rmtree('results/OKE')
        os.mkdir('results/OKE')
        save('results/OKE/OKE.npy', tOKE, wl, OKE)
        # extract OKE error 
        tOKE, df = load_fsTA_file(file = [i for i in os.listdir() if 'OKE' in i][0], error=True)
        OKE_error = np.mean(df, axis=0)
        save('results/OKE/Error.npy', tOKE, wl, OKE_error)
        # plot OKE trace
        c = Contour(files=['results/OKE/OKE.npy'], titles=['OKE'], cmap=['Blues'], white=[2], abs_white=[True], ylabel=r'$|\Delta A| / 10^{-3}$',
                    extend=['both'], scale=[scale_OKE], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=tlim_chirp, savefig='results/OKE/OKE.png',
                    yscale='linear', lines=[True], zeroline=False)        
        c.show()
        # plot kinetic traces
        from pyrene.standard.misc import rainbow
        from pyrene.fitting.fit_functions import gauss_conv_bi_exp
        colors = rainbow(wl_look)
        x = []
        y = []
        fwhm = []
        for i in range(len(wl_look)):
            k = Kinetics(files=['results/OKE/OKE.npy'], colors=[colors[i]], labels=[r'$\lambda_\text{probe} = %.3g\,\text{nm}$'%(wl_look[i])]
                         ,wavelength=[wl_look[i]], markersize=[3], marker=['o'], alphas=[0.4], xlim=tlim_chirp, show_only_for_short_time=True)
            trace = OKE[:, find_index(wl, wl_look[i])]
            p0 = [np.max(trace), 0.1, 0.2*np.max(trace), 3.0, tOKE[trace==np.max(trace)][0], 0.2]
            p, perr, chi2 = k.fit(file_index=0, model=gauss_conv_bi_exp, show_res=False, error='results/OKE/Error.npy', p0=p0)
            if len(p)>0 and 0.9<chi2<8.0:
                x.append(wl_look[i])
                y.append(p[-2])
                fwhm.append(p[-1])
            k.show()
        
        # plot time resolution
        fig, ax = plt.subplots(1,1)
        fwhm = np.array(fwhm)*1e3
        ax.plot(x, fwhm, 'ob', alpha=0.5)
        ax.axhline(y=np.mean(fwhm), linestyle='--', color='g', label=r'FWHM$_\text{mean} = %.3g$ fs'%(np.mean(fwhm)))
        ax.axhline(y=np.min(fwhm), linestyle='--', color='r',  label=r'FWHM$\text{min} = %.3g$ fs'%(np.min(fwhm)))
        ax.set_xlabel(r'$\lambda / \text{nm}$')
        ax.set_ylabel(r'FWHM / fs')
        ax.legend()
        ax.set_title('time resolution')
        fig.tight_layout()
        np.savetxt('results/OKE/timeresolution.txt', np.column_stack([x, fwhm]), header='wavelength / nm, time resolution / fs', delimiter=',')
        fig.savefig('results/OKE/timeresolution.png', transparent=True)
        plt.show(block=True)
        # write fit_params to file
        p, pcov = curve_fit(koboyashi, x, y)
        with open('ok.txt', 'w') as file:
            # Write each value in the list as a separate line in the file
            for item in p:
                file.write(str(item) + '\n')  
        # plot extracted points 
        c = Contour(files=['results/OKE/OKE.npy'], titles=['OKE'], cmap=['Blues'], white=[2], abs_white=[True], ylabel=r'$|\Delta A| / 10^{-3}$',
                    extend=['both'], scale=[scale_OKE], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=tlim_chirp, savefig='results/OKE/OKE_fit.png',
                    yscale='linear', lines=[True], zeroline=False)     
        c.ax.plot(wl, koboyashi(wl, *p), '--k')
        c.ax.plot(x, y, 'wo', markeredgecolor='k')  
        c.show()

# load and fit chirp
p = np.loadtxt('ok.txt')  
# plot chirp fit 
c = Contour(files=['results/average_minus_bg.npy'], titles=['chirp fit'], 
            extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=tlim_chirp,
            yscale='linear', lines=[True], zeroline=False, savefig='results/chirp_fit.png')
c.ax.plot(wl, koboyashi(wl, *p), '--k')
c.show()

# correct for chirp 
for i in range(len(wl)):
    t0 = koboyashi(wl[i], *p)
    dA[:,i] = np.interp(t, t-t0, dA[:,i])   
save('results/dA.npy', t, wl, dA)
# plot chirp corrected 
c = Contour(files=['results/dA.npy'], titles=['chirp corrected data'], 
            extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=[-0.5, 1.5],
            yscale='linear', lines=[True], zeroline=True, savefig='results/chirp_corrected.png')
c.show()
# plot on whole scale 
c = Contour(files=['results/dA.npy'], titles=['final TA spectrum'], 
            extend=['both'], scale=[scale], wn=False, figsize=(6, 4.5), nlevels=[31], ylim=[-1.0, 1800],
            yscale='symlog', lines=[True], zeroline=True, savefig='results/final_spec.png')
c.show()
# plot overview plot
s = Spectra(files=['results/dA.npy'], overview=True, outside=True, savefig='results/overview.png')
s.show()

# convert to txt
from pyrene.standard.misc import convert_to_txt
convert_to_txt('results/dA.npy', experiment='femto', t_cuts=None, wl_cuts=None)
shutil.move("wavelength.txt", "results/wavelength.txt")
shutil.move("time.txt", "results/time.txt")
shutil.move("TA.txt", "results/TA.txt")

if not save_scans:
    shutil.rmtree('results/scans')