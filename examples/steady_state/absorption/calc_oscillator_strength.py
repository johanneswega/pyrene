from pyrene.steady_state import Absorption

# concentration in M
conc = 24e-6
# pathlength in cm
pathlength = 1

# initialize absorption class
a = Absorption(files=['abs_file1.csv'],
                x_cuts=[(350, 800)],
                colors=['r'],
                c=[conc], 
                l=[pathlength],
                baseline_at=[700],
                labels=['your molecule'], savefig='Strickler.png')

# use calc_oscillator_strength(limits, n, nu0, file_index)
# integration limits as list in nm
# refractive index 
# center frequency in kk
a.calc_oscillator_strength([16, 22.5], 1.421, 18.8, 0)
a.show()