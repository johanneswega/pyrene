from pyrene.standard.packages import *

def conversion(p, pixel, A_TA, wl_ref, A_ref):
    """function to be fitted"""
    wl = p[0]*pixel+p[1]
    A_int = np.interp(wl_ref, wl, A_TA)
    return A_ref - A_int

def onclick(event, x, y, ax, fig):
    """function for interactive plot"""
    if event.button == 1:  # Left mouse button click
        x.append(event.xdata)
        y.append(event.ydata)
        # Add marker at the clicked coordinates
        ax.scatter(event.xdata, event.ydata, color='k', marker='x')
        # Refresh the plot
        fig.canvas.draw()

def pixel_to_lambda(WL_file='WL.dat', Ho_file='HOLMIUM.dat', lim=[0, 500], p0=[0.947, 274], auto=True):
    # read data 
    WL = np.loadtxt(WL_file)    
    Ho_TA = np.loadtxt(Ho_file)
    pixel = Ho_TA[:,0] 
    I_Ho = Ho_TA[:,2]
    I_WL = WL[:,2]
    # make negative values positive if present
    I_Ho[I_Ho<0] = -1*I_Ho[I_Ho<0]
    I_WL[I_WL<0] = -1*I_WL[I_WL<0]
    # calculate absorbance
    A_TA = -np.log10(I_Ho/I_WL)
    # reference absorption spectrum
    Ho_ref = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'HOLM_ref.csv'), skiprows=1, delimiter=',', usecols=[0,1])
    wl_ref = Ho_ref[:,0]
    A_ref = Ho_ref[:,1]
    A_ref = A_ref[wl_ref>300]
    wl_ref = wl_ref[wl_ref>300]

    # make figure
    if auto:
        # least squares fit
        res = least_squares(conversion, x0=p0,
                            args=(pixel[(pixel<lim[1])&(pixel>lim[0])], 
                                A_TA[(pixel<lim[1])&(pixel>lim[0])], wl_ref, A_ref))
        p = res.x
        scale = p[0]
        shift = p[1]
    else:
        # make first interactive figure
        fig, ax = plt.subplots(2,1,figsize=(12, 7))
        ax[0].plot(pixel, A_TA, '-r', label='Absorption TA setup')
        ax[0].set_xlabel('pixel')
        ax[0].set_ylim([0.02, 2.0])
        ax[0].set_ylabel('absorbance')
        ax[1].plot(wl_ref, A_ref, '-b', label='Lit. Absorption spectrum')
        ax[1].set_xlabel('wavelength / nm')
        ax[1].set_ylabel('absorbance')
        ax[0].legend()
        ax[1].legend()
        ax[0].set_title('Select peaks in TA setup spectrum!')
        x_TA = []
        y_TA = []
        fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, x_TA, y_TA, ax[0], fig))
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        fig.tight_layout()
        plt.show()
        x_TA = np.array(x_TA)
        y_TA = np.array(y_TA)

        # make second interactive figure
        fig, ax = plt.subplots(2,1,figsize=(12, 7))
        ax[0].plot(pixel, A_TA, '-r', label='Absorption TA setup')
        ax[0].plot(x_TA, y_TA, 'xk')
        ax[0].set_ylim([0.02, 2.0])
        ax[0].set_xlabel('pixel')
        ax[0].set_ylabel('absorbance')
        ax[1].plot(wl_ref, A_ref, '-b', label='Lit. Absorption spectrum')
        ax[1].set_xlabel('wavelength / nm')
        ax[1].set_ylabel('absorbance')
        ax[0].legend()
        ax[1].legend()
        ax[0].set_title('Select the same peaks in lit. spectrum!')
        x_lit = []
        y_lit = []
        fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, x_lit, y_lit, ax[1], fig))
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        fig.tight_layout()
        plt.show()
        x_lit = np.array(x_lit)
        y_lit = np.array(y_lit)

        # do linear regression
        fig, ax = plt.subplots(1,1)
        ax.plot(x_TA, x_lit, 'ob')
        scale, shift = np.polyfit(x_TA, x_lit, 1)
        wlfine = np.linspace(0, 520, 100)
        ax.plot(wlfine, scale*wlfine + shift, '--k', label = r'px = scale $\cdot \lambda$ + shift')
        ax.set_ylabel('wavelength / nm')
        ax.set_xlabel('pixel')
        ax.set_title(r'scale = %.3g, shift = %.3g'%(scale, shift))
        fig.tight_layout()
        plt.show()

    # plot results
    fig, ax = plt.subplots(1,1)
    title = 'scale = %.3g, shift = %.3g'%(scale,shift)

    # write results to file 
    wl = scale*pixel+shift
    # Open the file in write mode
    with open('wl.txt', 'w') as file:
        # Write each value in the list as a separate line in the file
        for item in wl:
            file.write(str(item) + '\n')

    ax.plot(scale*pixel+shift, A_TA, '-b', linewidth=1.5, label='TA setup')
    ax.plot(wl_ref, A_ref, '-r', linewidth=1.5, label='Ref.')
    ax.set_xlabel(r'$\lambda /$ nm')
    ax.set_ylabel(r'norm. Abs.')
    ax.set_xticks(np.linspace(300,800,6))
    ax.set_yticks([])
    ax.set_title(title)
    ax.legend()
    ax2 = ax.secondary_xaxis("top", functions=(lambda x: (x-shift)/scale,lambda x: x*scale + shift))
    ax2.set_xlabel('pixel')
    ax2.set_xticks(np.linspace(0,500,6))
    fig.tight_layout()
    fig.savefig('results/pixel_to_lambda.png')
    plt.show()