from pyrene.steady_state import Absorption

a = Absorption(
    files=['Data/FTIR_2.0'],
    x_cuts=[(1125, 1275)],
    baseline_at=[(1275)],
    norm=[True],
    IR=True,
    colors=['r'],
    labels=['your molecule'])

a.plot_calculated_spectrum(logfile='Data/Gaussian_output_IR.log', cuts=[1000, 1400], fwhm=10, 
                           scale=0.95, color='tomato', label='B3YLP-def2SVP')
a.show()