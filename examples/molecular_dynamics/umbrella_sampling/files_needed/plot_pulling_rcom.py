from pyrene.standard.packages import *
from pyrene.molecular_dynamics.analysis import *

plot_trajectory(files=['pull/com_dist.xvg'], colors=['b'], ylim=[0.1, 2.5], save='pull_trajectory.svg',
                find=np.arange(0.35, 2.05, 0.05).tolist(),
                labels=[r'$\text{BPe}^{\bullet -} \cdots \text{DPA}^{\bullet +}$ in DCM'], f_type='r_com')