import IPython

def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.constants as sc
    from numpy import sin, cos, tan, arcsin, arccos, arctan, pi as sin, cos, tan, arcsin, arccos, arctan, pi
    from numpy import log as ln 
    from numpy import log10 as lg
    from numpy import exp as exp 

    def wn(wl):
        return (1/wl)*10**(7)

    def wl(wn):
        return (1/wn)*10**(7)

    def sprint(a):
        print("%.5g"%a)

    def plot(x, y, mark=None, label=None, ylab=None, xlab=None, xlog=None, ylog=None):
        plt.plot(x, y, mark, label=label)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        if xlog==True:
            plt.xscale('log')
        if ylog==True:
            plt.yscale('log')
        if label!=None:
            plt.legend()
        plt.tight_layout()
        plt.show()

    IPython.start_ipython(
        argv=[],
        user_ns=locals()
    )
