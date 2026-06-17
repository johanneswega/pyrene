import numpy as np
import scipy.constants as sc
from matplotlib import pyplot as plt

# function to convert nm to kk and vice versa
def nm_to_kk(val):
    return (1/val)*10**4

# function to convert cm-1 to fs^-1
def cm_to_fs(val):
    return val*2.998e+10*10**(-15)

# function to convert cm-1 to ps^-1
def cm_to_ps(val):
    return (sc.c * (val*100))*10**(-12)

# FT from ps to cm-1
def FT_ps_to_cm(t, E):
    # Do the FFT
    fhat = np.fft.fftshift(np.fft.fft(E))
    freq = np.fft.fftshift(np.fft.fftfreq(len(t), d=t[1]-t[0]))
    f = freq[freq > 0]
    wn = (1 / (sc.c / f * 10**-12)) * 10**-2
    FT = np.abs(fhat[freq > 0])
    return wn, FT

# print fit parameters with erros
def print_params(p, pcov):
    print('')
    print('### Optimized Fit Parameters ###')
    print('')
    for i in range(len(p)):
        err = pcov[i,i]**(0.5)
        print('%i. %.3g +- %.3g'%(i+1, p[i], err))
    print('')

def sort_with_respect(sort_after, *args):
    # Sort all lists based on sort_after
    sorted_lists = sorted(zip(sort_after, *args))  
    
    # Unpack sorted tuples
    sorted_arrays = list(zip(*sorted_lists))  

    # Convert to numpy arrays and return
    return tuple(np.array(arr).tolist() for arr in sorted_arrays)
    
# Moving Average Filter
def moving_average(x, window_size):
    window = np.ones(window_size) / window_size
    return np.convolve(x, window, mode='valid')

# Function to calculate rainbow list
def rainbow(liste, r=None):
    if r!=None:
        return plt.cm.rainbow_r(np.linspace(0,1,len(liste))).tolist()
    else:
        return plt.cm.rainbow(np.linspace(0,1,len(liste))).tolist()
    
def gist_rainbow(liste, r=None):
    if r!=None:
        return plt.cm.gist_rainbow_r(np.linspace(0,1,len(liste))).tolist()
    else:
        return plt.cm.gist_rainbow(np.linspace(0,1,len(liste))).tolist()
    
# function to find index
def find_index(x, find):
    return np.argmin(np.abs(x - find))

# function for reading theoretical IR spectrum
def get_peaks(fname):
    x = []
    y = []
    with open(fname, "r") as file: 
        for line in file:
            # only read lines starting with #
            if line.startswith('#'):
                # split line
                line = line.split()
                # if second entry can be converted to a number add it
                try:
                    x.append(float(line[1]))
                    y.append(float(line[2]))
                except ValueError:
                    continue
    return np.array(x), np.array(y)

# load theoretical IR
def load_theo_IR(fname, cuts, fwhm, scale):
    # load peak information
    wnp, Ip = get_peaks(fname)

    wnp = scale*wnp

    # take only those peaks that lie within the range of cuts of exp. 
    Ip = Ip[(wnp>cuts[0]) & (wnp<cuts[1])]
    wnp = wnp[(wnp>cuts[0]) & (wnp<cuts[1])]
    
    # make frequency axis
    wn = np.linspace(cuts[0], cuts[1], 5000)
    I = np.zeros(len(wn))

    # construct spectrum
    from VautheyLab.fit_functions import lorentzian

    for i in range(len(Ip)):
        I += lorentzian(wn, Ip[i], wnp[i], fwhm)
    return (1/wn)*10**4, wn, I

# function to export theoIR to .txt as comp steady state in scripts
def export_theoIR(fname, cuts, fwhm, scale):
    head = r'wavenlength / µm,    wavenumber / cm-1,  absorbance'
    wl, wn, A = load_theo_IR(fname, cuts, fwhm, scale)  
    np.savetxt('TRIR_theoretical.txt', 
                np.column_stack([wl, wn, A]),
                header=head, delimiter=',')  

# calculate theoretical FTIR spectrum 
def calc_theo_TRIR(fGS, fES, cuts, fwhm, scale):
    head = r'wavenlength / µm,    wavenumber / cm-1,  absorbance'
    wl, wn, AGS = load_theo_IR(fGS, cuts, fwhm, scale)
    wl, wn, AES = load_theo_IR(fES, cuts, fwhm, scale)
    dA = AES - AGS
    np.savetxt('TRIR_theoretical.txt', 
                np.column_stack([wl, wn, dA]),
                header=head, delimiter=',')
    
# get sorted file after integer i.e. 2 before 10
def get_sorted_files(path):
    # Get all file names
    files = os.listdir(path)

    # get rid of .DS_store if present
    files = [i for i in files if not '.DS' in i]

    # Sort numerically (natural order)
    files_sorted = sorted(files, key=lambda x: int(os.path.splitext(x)[0]))
    return np.array([path + '/' + i for i in files_sorted])

# function to convert shida csv to txt to be used in absorption script
def shida_to_txt(shida, ma=False, units='wn'):
    data = np.loadtxt(shida, delimiter=',')
    if units=='wn':
        wn = data[:, 0]
        wl = nm_to_kk(wn)
    else:
        wl = data[:, 0]
        wn = nm_to_kk(wl)        
    A = data[:, 1]
    if ma!=False:
        wl = moving_average(wl, ma)
        wn = moving_average(wn, ma)
        A = moving_average(A, ma)
    # sort with respect to wavelengths
    #A = np.flip(A)
    #wn = np.flip(wn)
    wl, wn, A = sort_with_respect(wl, wn, A)
    np.savetxt('%s.txt'%(shida[:shida.find('.')]), 
                np.column_stack([wl, wn, A]),
                header='wavenlength / nm,    wavenumber / 10^3 cm-1,  absorbance', delimiter=',')
    
def get_data(file, delimiter=',', skiprows=1, usecols=[0,1], xcuts=None):
    data = np.loadtxt(file, delimiter=delimiter, skiprows=skiprows, usecols=usecols)
    x = data[:,0]
    y = data[:,1]
    x, y = sort_with_respect(x, y)
    if xcuts!=None:
        y = y[(x>xcuts[0]) & (x<xcuts[1])]
        x = x[(x>xcuts[0]) & (x<xcuts[1])]
    return x, y

def yline(ax, y, color='k', linestyle='--', alpha=0.3):
    ax.axhline(y=y, color=color, linestyle=linestyle, alpha=alpha)

def xline(ax, x, color='k', linestyle='--', alpha=0.3):
    ax.axvline(x=x, color=color, linestyle=linestyle, alpha=alpha)


def save(fname, t, wl, dA):
    np.save(fname, np.array([t, wl, dA], dtype=object))

def load(fname, t_cuts=None, wl_cuts=None):
    t, wl, dA = np.load(fname, allow_pickle=True)

    if wl_cuts!=None:
        dA = dA[:,(wl>wl_cuts[0])&(wl<wl_cuts[1])]
        wl = wl[(wl>wl_cuts[0])&(wl<wl_cuts[1])]

    if t_cuts!=None:
        dA = dA[(t>t_cuts[0])&(t<t_cuts[1]), :]
        t = t[(t>t_cuts[0])&(t<t_cuts[1])]
    
    return t, wl, dA

def interpolate_TA(file1, file2):
    # file for time axis
    t1, wl1, dA1 = load(file1)
    # file to be interpolated
    t2, wl2, dA2 = load(file2)
    # make empty dA frame
    dAint = np.zeros((len(t1), len(wl1)))   
    # go through wavelengths and interpolate
    for i in range(len(wl1)):
        # find index 
        index = find_index(wl2, wl1[i])
        dAint[:, i] = np.interp(t1, t2, dA2[:,index])
    # save interpolated file
    save('%s_int.npy'%(file2[:file2.find('.')]), t1, wl1, dAint)

# function to interpolate UVvis TA and TRIR on the same time grid
def interpolate_TA_TRIR(file1, file2):
    # load TRIR
    tIR, wnIR, dAIR = load_pdat(file1)
    # load vis
    tvis, wlvis, dAvis = load(file2)
    # for each wavelength we want to interplote the vis dA on the tIR axis
    dAint = np.zeros((len(tIR), len(wlvis)))
    for i in range(len(wlvis)):
        # get trace
        trace = dAvis[:, i]
        # interpolate
        int_trace = np.interp(tIR, tvis, trace)
        # add to array
        dAint[:, i] = int_trace
    save('%s_int.npy'%(file2[:file2.find('.')]), tIR, wlvis, dAint)

# function to interpolate two TRIR spectra on the same time grid
def interpolate_TRIR(file1, file2):
    # load TRIR
    tIR, wnIR, dAIR = load_pdat(file1)
    # load vis
    t, wn, dA = load_pdat(file2)
    # for each wavelength we want to interplote the vis dA on the tIR axis
    dAint = np.zeros((len(tIR), len(wn)))
    for i in range(len(wn)):
        # get trace
        trace = dA[:, i]
        # interpolate
        int_trace = np.interp(tIR, t, trace)
        # add to array
        dAint[:, i] = int_trace
    save_pdat('%s_int.pdat'%(file2[:file2.find('.')]), tIR, wn, dAint)

def save_pdat(name, t, wl, dA):
    pdat = np.zeros((len(t)+1, len(wl)+1))
    pdat[1:,0] = t
    pdat[0, 1:] = wl
    pdat[1:, 1:] = dA
    np.savetxt(name, pdat, header='ps*nm', delimiter=',')    

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA    

def average_pdats(files):
    dAs = []
    for i in range(len(files)):
        t, wl, dA = load_pdat(files[i])
        dAs.append(dA)
    dAs = np.array(dAs)
    dAmean = np.mean(dAs, axis=0)
    save_pdat('mean.pdat', t, wl, dAmean)

def npy_to_pdat(file):
    t, wl, dA = load(file)
    pdat = np.zeros((len(t)+1, len(wl)+1))
    pdat[1:,0] = t
    pdat[0, 1:] = wl
    pdat[1:, 1:] = dA
    np.savetxt('%s.pdat'%(file[:file.find('.')]), pdat, header='ps*nm', delimiter=',')

def pdat_to_npy(file):
    t, wl, dA = load_pdat(file)  
    np.save('%s.npy'%(file[:file.find('.')]), np.array([t, wl, dA], dtype=object))

def txt_to_npy(t_file, wl_file, dA_file):
    t = np.loadtxt(t_file, skiprows=1)
    wl = np.loadtxt(wl_file, skiprows=1)
    dA = np.loadtxt(dA_file, skiprows=1)
    save('dA_from_txt.npy', t, wl, dA)

def convert_to_txt(fname, experiment='femto', t_cuts=None, wl_cuts=None):
    if '.npy' in fname:
        t, l, dA = load(fname, t_cuts=t_cuts, wl_cuts=wl_cuts)
        # Save array a to a text file
        np.savetxt('wavelength.txt', l, header='wavelength / nm')

    elif '.pdat' in fname:
        t, l, dA = load_pdat(fname)
        # Save array a to a text file
        np.savetxt('wavenumber.txt', l, header='wavenumber / cm-1')

    if experiment=='nano':
        headert = 'time / ns'
    else:
        headert = 'time / ps'

    # Save array b to a text file
    np.savetxt('time.txt', t, header = headert)

    # Save array c to a text file
    np.savetxt('TA.txt', dA, header='Transient Absorption / mOD')        

def export(fname, spectra=None, kinetics=None, experiment='femto', t_cuts=None, wl_cuts=None):
    # load data 
    t, wl, dA = load(fname, t_cuts=t_cuts, wl_cuts=wl_cuts)
    wn = (1/wl)*10**4

    # go through all spectra and export 
    if not spectra==None:
        for i in range(len(spectra)):
            if experiment=='femto':
                if spectra[i] > 1000:
                    label = '%s_spectrum_at_%3g_ns.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])]/1000)
                elif np.abs(spectra[i]) < 1:
                    label = '%s_spectrum_at_%3g_fs.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])]*1000)
                else:
                    label = '%s_spectrum_at_%3g_ps.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])])
            else:
                if spectra[i] > 1000:
                    label = '%s_spectrum_at_%3g_us.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])]/1000)
                elif np.abs(spectra[i]) < 1:
                    label = '%s_spectrum_at_%3g_ps.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])]*1000)  
                else:
                    label = '%s_spectrum_at_%3g_ns.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])])           
                
            np.savetxt(label, np.column_stack([wl, wn, dA[find_index(t, spectra[i]), :]]),
                        header='wavelength / nm,     wavenumber / kK,        TA / mOD', delimiter=',')
        
    # go through all kinetics and export
    if not kinetics==None:
        for i in range(len(kinetics)):
            label = '%s_kinetics_at_%3g_nm.txt'%(fname[:fname.find('.')], wl[find_index(wl, kinetics[i])])
            if experiment=='femto':
                header = 'time / ps,      TA / mOD'
            else:
                header = 'time / ns,      TA / mOD'
            np.savetxt(label, np.column_stack([t, dA[:, find_index(wl, kinetics[i])]]),
            header=header, delimiter=',')
