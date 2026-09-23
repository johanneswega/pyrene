from pyrene.standard.packages import *
from pyrene.standard.misc import moving_average as ma
from pyrene.standard.misc import find_index

# function to get data 
def get_data(file, f_type):
    if f_type=='r_com':
        data = np.loadtxt(file, skiprows=24)
        t = data[:,0]
        r_com = data[:,1]
        return t, r_com
    if f_type=='shells':
        data = np.loadtxt(file, skiprows=25)
        t = data[:,0]
        n = data[:,1]
        return t, n       

# function to find frame corresponding to certain property
def find_frame(t, y, target, dt, nstxout_compressed):

    # Find the frame closest to the target angle
    i = np.argmin(np.abs(y - target))

    # dt = Simulation time step in ps

    # nstxout_compressed Coordinates saved every xxxx MD steps

    # Time between trajectory frames
    frame_dt = dt * nstxout_compressed  

    # Corresponding trajectory frame
    frame = int(round(t[i] / frame_dt))
    print('')
    print(f'Time = {t[i]:.2f} ps')
    print(f'Frame = {frame}')
    print(f'property = {y[i]:.2f}')
    print('')
    return int(round(t[i]))

# function to plot / compare trajectories
def plot_trajectory(files, colors, labels, f_type, figsize=(6, 3.5), find=None, ylim=None,
                    xlim=False, MA=False, MA_npoints=10, outside=True, linewidth=1, save=None,
                    dt=0.002, nstxout_compressed=1000):
    # make figure
    fig, ax = plt.subplots(1,1,figsize=figsize)
    # get data and plot
    for i in range(len(files)):
        x, y = get_data(files[i], f_type)
        if MA==False:
            ax.plot(x, y, '-', color=colors[i], label=labels[i], linewidth=linewidth)
        else:
            ax.plot(x, y, '-', color=colors[i], alpha=0.2, linewidth=linewidth)
            ax.plot(ma(x,MA_npoints), ma(y,MA_npoints), '-', color=colors[i], label=labels[i], linewidth=linewidth)
        if find:
            times_found = []
            for j in range(len(find)):
                ti = find_frame(x, y, find[j], dt, nstxout_compressed)
                times_found.append(ti)
            np.savetxt('times_found.txt', np.column_stack([times_found, find]), header='time / ps, rcom/ nm')
    if f_type=='r_com':
        ax.set_ylabel(r'$r_{\text{com}} / $ nm')
        ax.set_xlabel(r'$t / $ ps')
    if xlim!=None:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if outside==True:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    else:
        ax.legend()
    fig.tight_layout()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()

# function to plot / compare trajectories using a stack
def plot_trajectory_stack(files, colors, labels, f_type, figsize=None, save=None, ylim=False,
                          xlim=False, MA=False, MA_npoints=10, outside=True, linewidth=1):
    # make figure
    if figsize==None:
        figsize = (6, len(files)*3)
    fig, ax = plt.subplots(len(files), 1, figsize=figsize, sharex=True, sharey=True)
    # get data and plot
    for i in range(len(files)):
        x, y = get_data(files[i], f_type)
        if MA==False:
            ax[i].plot(x, y, '-', color=colors[i], label=labels[i], linewidth=linewidth)
        else:
            ax[i].plot(x, y, '-', color=colors[i], alpha=0.2, linewidth=linewidth)
            ax[i].plot(ma(x,MA_npoints), ma(y,MA_npoints), '-', color=colors[i], label=labels[i], linewidth=linewidth)
        if f_type=='r_com':
            ax[i].set_ylabel(r'$r_{\text{com}} / $ nm')
        if outside==True:
            ax[i].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        else:
            ax[i].legend()
    ax[-1].set_xlabel(r'$t / $ ps')
    if xlim!=None:
        ax[-1].set_xlim(xlim)
    if ylim!=False:
        ax[-1].set_ylim(ylim)
    fig.tight_layout()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()

# histogram plotter
def histogram(files, colors, labels, f_type, nbins, width, figsize=(6, 3.5), outside=True, ylim=False, save=None, find=None, 
              dt=0.002, nstxout_compressed=1000, xticks=None):
    # make figure
    fig, ax = plt.subplots(1,1,figsize=figsize)
    # get data and plot
    for i in range(len(files)):
        x, y = get_data(files[i], f_type)
        if f_type == 'shells':
            # Discrete integer values: 0, 1, 2, ...
            values, counts = np.unique(y.astype(int), return_counts=True)
            ax.bar(values,counts, width=width, color=colors[i], edgecolor='k', align='center',label=labels[i])
        else:
            ax.hist(y, nbins, color=colors[i],width=width, edgecolor='k', label=labels[i])
        if find:
            for j in range(len(find)):
                find_frame(x, y, find[j], dt, nstxout_compressed)
    if f_type=='r_com':
        ax.set_xlabel(r'$r_{\text{com}} / $ nm')
    if f_type=='shells':
        ax.set_xlabel(r'$n$')
    ax.set_ylabel(r'$N$')
    if outside==True:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    else:
        ax.legend()
    if ylim!=False:
        ax.set_ylim(ylim)
    if xticks:
        ax.set_xticks(xticks)
    fig.tight_layout()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()

# histogram plotter
def histogram_stack(files, colors, labels, f_type, nbins, width, figsize=None, outside=True, ylim=False, save=None, prob=False):
    # make figure
    if figsize==None:
        figsize = (6, len(files)*3)
    fig, ax = plt.subplots(len(files), 1, figsize=figsize, sharex=True, sharey=True)
    # get data and plot
    for i in range(len(files)):
        x, y = get_data(files[i], f_type)
        if prob:
            ax[i].hist(y, nbins, color=colors[i], width=width, edgecolor='k', label=labels[i], density=True)
            ax[i].set_ylabel(r'$p$')
        else:
            ax[i].hist(y, nbins, color=colors[i], width=width, edgecolor='k', label=labels[i])
            ax[i].set_ylabel(r'$N$')
        if outside==True:
            ax[i].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        else:
            ax[i].legend()
    if f_type=='r_com':
        ax[-1].set_xlabel(r'$r_{\text{com}} / $ nm')
    if f_type=='shells':
        ax[-1].set_xlabel(r'$n$')
    if ylim!=False:
        ax[-1].set_ylim(ylim)
    fig.tight_layout()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()

# function to plot RDF
def plot_rdf(files,  colors, labels, figsize=(6, 3.5), outside=True, ylim=False, xlim=False, save=None):
    # make figure
    fig, ax = plt.subplots(1,1, figsize=figsize)
    ax.axhline(y=1, color='k', linestyle='--', alpha=0.2)
    for i in range(len(files)):
        data = np.loadtxt(files[i], skiprows=25)
        r = data[:,0]
        g = data[:,1]
        ax.plot(r, g, '-', color=colors[i], label = labels[i])
    if ylim:
        ax.set_ylim(ylim)
    if xlim: 
        ax.set_xlim(xlim)
    ax.set_ylabel(r'$g(r)$')
    ax.set_xlabel(r'$r_\text{com} / \text{nm}$')
    fig.tight_layout()
    ax.legend()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()

# function to plot pmf
def plot_pmf(files,  colors, labels, sub, figsize=(6, 3.5), outside=True, ylim=False, xlim=False, save=None):
    # make figure
    fig, ax = plt.subplots(1,1, figsize=figsize)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.2)
    for i in range(len(files)):
        data = np.loadtxt(files[i], skiprows=17)
        r = data[:,0]
        # convert from kJ/mol to eV
        w = data[:,1] * 0.01036
        w -= sub[i]
        ax.plot(r, w, '-', color=colors[i], label = labels[i])
    if ylim:
        ax.set_ylim(ylim)
    if xlim: 
        ax.set_xlim(xlim)
    ax.set_ylabel(r'$w(r)$ / eV')
    ax.set_xlabel(r'$r_\text{com} / \text{nm}$')
    fig.tight_layout()
    ax.legend()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()