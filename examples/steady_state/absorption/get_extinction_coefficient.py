from pyrene.steady_state import Absorption

a = Absorption(
    files=['Data/abs_file1.csv'],
    x_cuts=[(350, 800)], c=[24e-6], l=[1],
    wn=False,
    baseline_at=[700],
    colors=['r'],
    labels=['your molecule'])

# find wavelength / extinction coefficient at maximum absorbance 
a.find(0)
# find absorbance / extinction coefficient at specific wavelength
a.find(0, where=409)

a.show()