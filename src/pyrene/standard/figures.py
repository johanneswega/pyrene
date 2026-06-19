import numpy as np
from matplotlib import pyplot as plt
np.seterr(divide='ignore')
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))

def figure_spectum(figsize=(8,5), IR=False, wn=True):
    """Creates standard spectrum figure (wavenumber vs. y) with double axis (wavelength on top) if IR==False
        If IR == True: no double axis"""
    fig, ax = plt.subplots(1,1,figsize=figsize)
    if not IR:
        if wn:
            ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
            ax.invert_xaxis()    
            ax.axhline(y=0, color='k')
            axsec = ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
            axsec.set_xlabel(r'$\lambda / $ nm') 
        else:
            ax.set_xlabel(r'$\lambda / $ nm')
    else:
        ax.set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
    ax.axhline(y=0, color='k')
    return fig, ax

def figure_kinetics(figsize=(8,5), experiment='femto'):
    """Creates standard kinetics figure (time vs. y)"""
    fig, ax = plt.subplots(1,1,figsize=figsize)
    if experiment=='femto':
        ax.set_xlabel(r'$\Delta t$ / ps')
    else:
        ax.set_xlabel(r'$\Delta t$ / ns')
    return fig, ax

def figure_contour(t, wl, dA, scale, figsize=(8,5), experiment='femto',
                   ylim=None, lines=True, nlines=25):
    """Creates standard contour figure (time vs. wavenumber vs. y)"""
    wn = (1/wl)*10**4
    # make figure
    fig, ax = plt.subplots(1,1,figsize=figsize)
    if lines==True:
        ax.contour(wn, t, dA, levels=nlines, linewidths=0.12, colors='k')
    D = ax.pcolormesh(wn, t, dA, vmin=scale[0], vmax=scale[1], cmap="RdBu_r")
    cbar = plt.colorbar(D, ax=ax)
    cbar.set_label(r'$\Delta{A} / 10^{-3}$')  
    ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
    ax.invert_xaxis()   
    axsec = ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
    axsec.set_xlabel(r'$\lambda / $ nm') 
    if experiment=='femto':
        ax.set_ylabel(r'$\Delta t / \text{ps}$')
    if ylim!=None:
        ax.set_ylim(ylim)
    return fig, ax

def solvchrom_figures(wns, solvs, cols, label='Stokes', stokes=False, lim=None, save=False):
    """method to make typical solvatochromism correlation figures"""
    from pyrene.standard.solvents import solvent_dic
    fig1, ax1 = plt.subplots(1,1,figsize=(5, 3.5))
    fig2, ax2 = plt.subplots(1,1,figsize=(5, 3.5))
    fig3, ax3 = plt.subplots(1,1,figsize=(5, 3.5))
    for i in range(len(solvs)):
        n = solvent_dic[solvs[i]][0]
        er = solvent_dic[solvs[i]][1]
        fe = (2*(er - 1)/(2*er+1))
        fn = (2*(n**2 - 1)/(2*n**2+1))
        df = fe - fn
        ax1.plot(df, wns[i], 'o', color=cols[i])
        ax2.plot(fe, wns[i], 'o', color=cols[i])
        ax3.plot(fn, wns[i], 'o', color=cols[i])
    if not stokes:
        text = r'$\Tilde{\nu}_{\text{max, %s}} / 10^3$ cm$^{-1}$'%(label)
    else:
        text = r'$\Delta \Tilde{\nu}_{\text{Stokes}} /$ cm$^{-1}$'
    ax1.set_ylabel(text)
    ax2.set_ylabel(text)
    ax3.set_ylabel(text)
    ax1.set_xlabel(r'$\Delta f = f(\varepsilon_r) - f(n^2)$')
    ax2.set_xlabel(r'$f(\varepsilon_r)$')
    ax3.set_xlabel(r'$f(n^2)$')
    if lim!=None:
        ax1.set_ylim(lim)
        ax2.set_ylim(lim)
        ax3.set_ylim(lim)
    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    if save!=False:
        fig1.savefig('%s_df.svg'%label, transparent=True)
        fig2.savefig('%s_f_epsr.svg'%label, transparent=True)
        fig3.savefig('%s_f_n2.svg'%label, transparent=True)
    return n, er, fe, fn, df