from pyrene.steady_state import Absorption

# concentration in M
conc = 24e-6
# pathlength in cm
pathlength = 1

# initialize absorption class
a = Absorption(files=['Data/abs_file1.csv'],
                x_cuts=[(350, 800)],
                colors=['r'],
                c=[conc], 
                l=[pathlength],
                baseline_at=[700],
                labels=['your molecule'])

a.plot_calculated_spectrum(logfile='Data/Gaussian_output_TDDFT.log', cuts=[15, 40], fwhm=2, norm=False,
                           shift=-4, color='gray', label='TD-DFT B3YLP-def2SVP')
a.show()