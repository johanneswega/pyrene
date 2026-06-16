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
