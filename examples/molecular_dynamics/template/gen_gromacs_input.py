import os 
import shutil
import re
from scipy.constants import Avogadro as Na

# specify molname and solvent
molname = 'BPe'
solv = 'ACN'

### solvent paramters taken from ###
### From GROMACS Molecular Dynamics Simulations to Electronic Absorption Spectra: A Tutorial for Small Molecules in Organic Solvents ###
### Picarelli et al. J. Chem. Educ. (2026) ###
# solvent density / g/mL
rho = 0.776
# isothermal compressibility / bar^-1
alpha = '1.070e-4'
# molar mass in g/mol
M = 41.05

# Length of cubic simulation box in nm
L = 5
# Temperature 
T = '298.15'

# Energy Minimizatin Parameters #
integrator_min = 'steep' 
nsteps_min = '50000'
emstep_min = '0.1'
emtol_min = '500'
cutoff_scheme_min = 'Verlet'
ns_type_min = 'grid'
nstlist_min  = '10'          
rcoulomb_min = '1.2'         
rvdw_min = '1.2'        
rlist_min = '1.2'
coulombtype_min = 'PME'        
pme_order_min = '4'         
fourierspacing_min = '0.12'   
pbc_min = 'xyz'       

# NVT Equilibration Parameters #
integrator_nvt = 'md' 
dt_nvt = '0.001'
nsteps_nvt = '500000'
continuation_nvt = 'no'
nstenergy_nvt = '100'
nstlog_nvt = '100'
constraint_algorithm_nvt = 'lincs'
constraints_nvt = 'h-bonds'
lincs_iter_nvt = '1'
lincs_order_nvt = '4'
rlist_nvt = '1.2'
rcoulomb_nvt = '1.2'
rvdw_nvt = '1.2'
coulombtype_nvt = 'PME'
pme_order_nvt = '4'
fourierspacing_nvt = '0.16'
tcoupl_nvt = 'v-rescale'
tc_grps_nvt = 'system'
tau_t_nvt = '0.1'
pcoupl_nvt = 'no'
pbc_nvt = 'xyz'  
DispCorr_nvt = 'EnerPres'
gen_vel_nvt = 'yes'
gen_seed_nvt = '-1'

# NPT Equilibration Parameters #
integrator_npt = 'md'
dt_npt = '0.001'
nsteps_npt = '500000'
continuation_npt = 'yes'
nstenergy_npt = '100'
nstlog_npt = '100'
nstxout_compressed_npt = '100'
constraint_algorithm_npt = 'lincs' 
constraints_npt = 'h-bonds'
lincs_iter_npt = '1'
lincs_order_npt = '4'
cutoff_scheme_npt = 'Verlet' 
rcoulomb_npt = '1.2'
rvdw_npt = '1.2'
rlist_npt = '1.2'
coulombtype_npt = 'PME'
pme_order_npt = '4'
fourierspacing_npt = '0.16'
tcoupl_npt = 'v-rescale'
tc_grps_npt = 'system'
tau_t_npt = '0.5' 
pcoupl_npt = 'C-rescale'
pcoupltype_npt = 'isotropic'
tau_p_npt = '5.0'
ref_p_npt = '1.0'
pbc_npt = 'xyz'
DispCorr_npt = 'EnerPres'
gen_vel_npt = 'no'

# Production Parameters #
integrator_prod = 'md'
dt_prod = '0.002'       
nsteps_prod = '25000000'    
continuation_prod = 'yes'
nstenergy_prod = '1000'            
nstlog_prod = '1000'      
nstxout_compressed_prod = '1000'      
constraint_algorithm_prod  = 'lincs' 
constraints_prod = 'h-bonds'     
lincs_iter_prod = '1'          
lincs_order_prod = '4'          
ns_type_prod = 'grid'       
nstlist_prod = '10'          
rcoulomb_prod = '1.2'         
rvdw_prod = '1.2'         
rlist_prod = '1.2'         
coulombtype_prod = 'PME'         
pme_order_prod = '4'           
fourierspacing_prod = '0.16'      
tcoupl_prod = 'Nose-Hoover'   
tc_grps_prod = 'system'   
tau_t_prod = '0.5'              
pcoupl_prod = 'C-rescale'     
pcoupltype_prod = 'isotropic'             
tau_p_prod = '5.0'                   
ref_p_prod = '1.0'         
pbc_prod = 'xyz'         
DispCorr_prod = 'EnerPres'    
gen_vel_prod = 'no' 

# calculate using solvent density how many solvent molecules to add
# legth of the box in / m 
Lm = L*1e-9
# volume in m^3
V = Lm**3
# convert to mL
V = V*10**6
# calculate the required number of ACN molecules
N = round((rho*Na*V)/M)

### STEP 0: make box script ###
content = f'''#!/bin/bash

# make box 
gmx editconf -f {molname}.gro -bt cubic -box {L} {L} {L} -o box_{molname}.gro -center {L/2} {L/2} {L/2}

# add solvent
gmx insert-molecules -f box_{molname}.gro -ci {solv}.gro -nmol {N} -o box_{molname}_{solv}.gro -try 10000'''

with open('step0_make_box.sh', 'w') as file:
    file.write(content)

### STEP 1: Energy minimization ###

# assemble min.mdp file 
if not os.path.exists('min'):
    os.mkdir('min')

content = f'''; Run parameters
integrator            = {integrator_min}       ; steepest descent 
nsteps                = {nsteps_min}      ; max minimization steps
emstep                = {emstep_min}         ; step size
emtol                 = {emtol_min}       ; minimization stop energy

; Neighbor searching
cutoff-scheme          = {cutoff_scheme_min}
ns_type                = {ns_type_min}       ; search neighboring grid cells
nstlist                = {nstlist_min}          ; update freq for neighbor list and long range forces
rcoulomb               = {rcoulomb_min}        ; short-range electrostatic cutoff (in nm)
rvdw                   = {rvdw_min}       ; short-range van der Waals cutoff (in nm)
rlist                  = {rlist_min}        ; Short-range neighbor list

; Electrostatics
coulombtype            = {coulombtype_min}        ; Particle Mesh Ewald for long-range electrostatics
pme_order              = {pme_order_min}          ; cubic interpolation
fourierspacing         = {fourierspacing_min}       ; grid spacing for FFT

; Periodic boundary conditions
pbc                    = {pbc_min}        ; 3-D PBC
'''

with open('min.mdp', 'w') as file:
    file.write(content)

content=f'''#!/bin/bash

# energy minimization preprocessor
gmx grompp -f min.mdp -c box_{molname}_{solv}.gro -p topol.top -o min/min

# peform energy minimization
gmx mdrun -v -deffnm min/min

# export potential energy
echo 10 | gmx energy -f min/min.edr -o min/ener.xvg

# plot potential energy
python min/plot_epot.py'''

with open('step1_energy_minimization.sh', 'w') as file:
    file.write(content)

# make plot file for energy
content_epot = r'''from pyrene.standard.packages import *
from pyrene.standard.misc import moving_average as ma

data = np.loadtxt('min/ener.xvg', skiprows=24)
step = data[:,0]
Epot = data[:,1]

fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
ax.plot(step, Epot, '-b')
ax.set_ylabel(r'$E_{pot}$ / kJ $\cdot$ mol$^{-1}$')
ax.set_xlabel(r'step')
fig.tight_layout()
fig.savefig('min/Epot.png', dpi=200)'''

with open('min/plot_epot.py', 'w') as file:
    file.write(content_epot)

### Step 2 NVT Equilibration ###
# generate nvt.mdp file and folder
if not os.path.exists('nvt'):
    os.mkdir('nvt')

content = f'''; Run parameters
integrator            = {integrator_nvt}          ; leap-frog integrator
dt                    = {dt_nvt}      ; time step length 
nsteps                = {nsteps_nvt}      ; number of steps
continuation          = {continuation_nvt}          ; continue from coordinates but generate new velocities
 
; Output control
nstenergy             = {nstenergy_nvt}         ; steps between saving energies          
nstlog                = {nstlog_nvt}         ; steps between saving log

; Bond parameters
constraint_algorithm  = {constraint_algorithm_nvt} 
constraints           = {constraints_nvt}     ; H bonds constrained
lincs_iter            = {lincs_iter_nvt}           ; accuracy of LINCS
lincs_order           = {lincs_order_nvt}           ; also related to accuracy

; Neighborsearching
rlist       = {rlist_nvt}      ; short-range neighborlist cutoff (in nm)
rcoulomb    = {rcoulomb_nvt}       ; short-range electrostatic cutoff (in nm)
rvdw        = {rvdw_nvt}       ; short-range van der Waals cutoff (in nm)

; Electrostatics
coulombtype           = {coulombtype_nvt}        ; Particle Mesh Ewald for long-range electrostatics
pme_order             = {pme_order_nvt}           ; cubic interpolation
fourierspacing        = {fourierspacing_nvt}        ; grid spacing for FFT

; Temperature coupling
tcoupl                = {tcoupl_nvt}   ; modified Berendsen thermostat
tc-grps               = {tc_grps_nvt}      ; couple all the system together
tau_t                 = {tau_t_nvt}         ; time constant, in ps
ref_t                 = {T}      ; reference temperature in K

; Pressure coupling 
pcoupl                = {pcoupl_nvt}

; Periodic boundary conditions
pbc                   = {pbc_nvt}         ; 3-D PBC

; Dispersion correction
DispCorr              = {DispCorr_nvt}    ; account for cut-off vdW scheme

; Velocity generation
gen_vel               = {gen_vel_nvt}
gen_temp              = {T}      ; temperature for Maxwell distribution
gen_seed              = {gen_seed_nvt}          ; generate a random seed'''

with open('nvt.mdp', 'w') as file:
    file.write(content)

content=f'''#!/bin/bash

# NVT equilibration preprocessor
gmx grompp -f nvt.mdp -c min/min.gro -p topol.top -o nvt/nvt

# NVT equilibration
gmx mdrun -v -deffnm nvt/nvt

# export potential energy
echo 11 | gmx energy -f nvt/nvt.edr -o nvt/ener.xvg 

# export temperature
echo 15 | gmx energy -f nvt/nvt.edr -o nvt/temp.xvg 

# plot energy
python nvt/plot_epot.py

# plot temperature
python nvt/plot_T.py'''

with open('step2_NVT_equilibration.sh', 'w') as file:
    file.write(content)

# make plot file for energy
content_T = r'''from pyrene.standard.packages import *
from pyrene.standard.misc import moving_average as ma

data = np.loadtxt('nvt/temp.xvg', skiprows=24)
t = data[:,0]
T = data[:,1]

fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
ax.plot(t, T, '-r')
ax.plot(ma(t,50), ma(T,50), '-', color='orange')
ax.axhline(y=298.15, linestyle='--', color='k')
ax.set_ylabel(r'$T$ / K')
ax.set_xlabel(r'$t$ / ps')
fig.tight_layout()
fig.savefig('nvt/NPT_T.png', dpi=200)'''

# make plot file for energy
content_epot = r'''from pyrene.standard.packages import *
from pyrene.standard.misc import moving_average as ma

data = np.loadtxt('nvt/ener.xvg', skiprows=24)
step = data[:,0]
Epot = data[:,1]

fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
ax.plot(step, Epot, '-b')
ax.set_ylabel(r'$E_{pot}$ / kJ $\cdot$ mol$^{-1}$')
ax.set_xlabel(r'step')
fig.tight_layout()
fig.savefig('nvt/NPT_Epot.png', dpi=200)'''

with open('nvt/plot_epot.py', 'w') as file:
    file.write(content_epot)

with open('nvt/plot_T.py', 'w') as file:
    file.write(content_T)

### Step 3 NPT Equilibration ###
# generate npt.mdp file and folder
if not os.path.exists('npt'):
    os.mkdir('npt')

content = f'''; Run parameters
integrator            = {integrator_npt}          ; leap-frog integrator
dt                    = {dt_npt}       ; time step length 
nsteps                = {nsteps_npt}      ; number of steps
continuation          = {continuation_npt}
 
; Output control
nstenergy             = {nstenergy_npt}         ; steps between saving energies          
nstlog                = {nstlog_npt}         ; steps between saving log
nstxout-compressed    = {nstxout_compressed_npt}         ; steps between saving compressed coords

; Bond parameters
constraint_algorithm  = {constraint_algorithm_npt} 
constraints           = {constraints_npt}     ; H bonds constrained
lincs_iter            = {lincs_iter_npt}           ; accuracy of LINCS
lincs_order           = {lincs_order_npt}           ; also related to accuracy

; Neighborsearching
cutoff-scheme         = {cutoff_scheme_npt}
rcoulomb              = {rcoulomb_npt}         ; short-range electrostatic cutoff (in nm)
rvdw                  = {rvdw_npt}         ; short-range van der Waals cutoff (in nm)
rlist                 = {rlist_npt}         ; Short-range neighbor list

; Electrostatics
coulombtype           = {coulombtype_npt}         ; Particle Mesh Ewald for long-range electrostatics
pme_order             = {pme_order_npt}           ; cubic interpolation
fourierspacing        = {fourierspacing_npt}       ; grid spacing for FFT

; Temperature coupling
tcoupl                = {tcoupl_npt}   ; modified Berendsen thermostat
tc-grps               = {tc_grps_npt}   
tau_t                 = {tau_t_npt}         ; time constant, in ps
ref_t                 = {T}     ; reference temperature in K

; Pressure coupling is on
pcoupl                  = {pcoupl_npt}     ; Pressure coupling on in NPT
pcoupltype              = {pcoupltype_npt}             ; uniform scaling of box vectors
tau_p                   = {tau_p_npt}                   ; time constant, in ps
ref_p                   = {ref_p_npt}                   ; reference pressure, in bar
compressibility = {alpha}      ; isothermal compressibility, bar^-1 for {solv} 
;refcoord_scaling = com

; Periodic boundary conditions
pbc                   = {pbc_npt}         ; 3-D PBC

; Dispersion correction
DispCorr              = {DispCorr_npt}    ; account for cut-off vdW scheme

; Velocity generation
gen_vel               = {gen_vel_npt} 
'''

with open('npt.mdp', 'w') as file:
    file.write(content)

content=f'''#!/bin/bash

# NPT equilibration preprocessor
gmx grompp -f npt.mdp -c nvt/nvt.gro -p topol.top -o npt/npt

# NPT equilibration
gmx mdrun -v -deffnm npt/npt

# export potential energy
echo 11 | gmx energy -f npt/npt.edr -o npt/ener.xvg
# plot potential energy
python npt/plot_epot.py

# export temperature
echo 15 | gmx energy -f npt/npt.edr -o npt/temp.xvg
# plot temperature
python npt/plot_T.py

# export density
echo 23 | gmx energy -f npt/npt.edr -o npt/density.xvg
# plot density
python npt/plot_rho.py

# export pressure
echo 17 | gmx energy -f npt/npt.edr -o npt/pressure.xvg
# plot pressure
python npt/plot_p.py
'''

with open('step3_NPT_equilibration.sh', 'w') as file:
    file.write(content)

# make plot file for energy
content_T = r'''from pyrene.standard.packages import *
from pyrene.standard.misc import moving_average as ma

data = np.loadtxt('npt/temp.xvg', skiprows=24)
t = data[:,0]
T = data[:,1]

fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
ax.plot(t, T, '-r')
ax.plot(ma(t,50), ma(T,50), '-', color='orange')
ax.axhline(y=298.15, linestyle='--', color='k')
ax.set_ylabel(r'$T$ / K')
ax.set_xlabel(r'$t$ / ps')
fig.tight_layout()
fig.savefig('npt/NPT_T.png', dpi=200)'''

# make plot file for energy
content_epot = r'''from pyrene.standard.packages import *
from pyrene.standard.misc import moving_average as ma

data = np.loadtxt('npt/ener.xvg', skiprows=24)
step = data[:,0]
Epot = data[:,1]

fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
ax.plot(step, Epot, '-b')
ax.set_ylabel(r'$E_{pot}$ / kJ $\cdot$ mol$^{-1}$')
ax.set_xlabel(r'step')
fig.tight_layout()
fig.savefig('npt/NPT_Epot.png', dpi=200)'''

# make plot file for density
content_rho = r'''from pyrene.standard.packages import *
from pyrene.standard.misc import moving_average as ma

data = np.loadtxt('npt/density.xvg', skiprows=24)
t = data[:,0]
rho = data[:,1]

fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
ax.plot(t, rho, '-g')
ax.plot(ma(t,50), ma(rho,50), '-', color='limegreen')
ax.axhline(y=np.mean(rho[-200:-1]), linestyle='--', color='orange')
ax.axhline(y=%.3g, linestyle='--', color='k')
ax.set_ylim([%.3g, %.3g])
ax.set_ylabel(r'$\rho$ / kg$\cdot$m$^{-3}$')
ax.set_xlabel(r'$t$ / ps')
fig.tight_layout()
fig.savefig('npt/NPT_rho.png', dpi=200)'''%(rho*1000, rho*1000 - 300, rho*1000 + 300)

# make plot file for pressure
content_P = r'''from pyrene.standard.packages import *
from pyrene.standard.misc import moving_average as ma

data = np.loadtxt('npt/pressure.xvg', skiprows=24)
t = data[:,0]
p = data[:,1]

fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
ax.plot(t, p, '-', color='purple')
ax.axhline(y=1, linestyle='--', color='k')
ax.set_ylabel(r'$p$ / bar')
ax.set_xlabel(r'$t$ / ps')
fig.tight_layout()
fig.savefig('npt/NPT_P.png', dpi=200)'''

with open('npt/plot_rho.py', 'w') as file:
    file.write(content_rho)

with open('npt/plot_epot.py', 'w') as file:
    file.write(content_epot)

with open('npt/plot_T.py', 'w') as file:
    file.write(content_T)

with open('npt/plot_P.py', 'w') as file:
    file.write(content_P)

### Step 4 Production ###
if not os.path.exists('prod'):
    os.mkdir('prod')

content = f''' ; Run parameters
integrator            = {integrator_prod}        ; leap-frog integrator
dt                    = {dt_prod}       ; time step length 
nsteps                = {nsteps_prod}    ; number of steps
continuation          = {continuation_prod}
; Output control
nstenergy             = {nstenergy_prod}        ; steps between saving energies          
nstlog                = {nstlog_prod}      ; steps between saving log
nstxout-compressed    = {nstxout_compressed_prod}       ; steps between saving compressed coords

; Bond parameters
constraint_algorithm  = {constraint_algorithm_prod} 
constraints           = {constraints_prod}     ; H bonds constrained
lincs_iter            = {lincs_iter_prod}           ; accuracy of LINCS
lincs_order           = {lincs_order_prod}           ; also related to accuracy

; Neighborsearching
ns_type               = {ns_type_prod}        ; search neighboring grid cells
nstlist               = {nstlist_prod}          ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb              = {rcoulomb_prod}         ; short-range electrostatic cutoff (in nm)
rvdw                  = {rvdw_prod}        ; short-range van der Waals cutoff (in nm)
rlist                 = {rlist_prod}         ; Short-range neighbor list

; Electrostatics
coulombtype           = {coulombtype_prod}         ; Particle Mesh Ewald for long-range electrostatics
pme_order             = {pme_order_prod}           ; cubic interpolation
fourierspacing        = {fourierspacing_prod}        ; grid spacing for FFT

; Temperature coupling
tcoupl                = {tcoupl_prod}   ; thermostat
tc-grps               = {tc_grps_prod}   
tau_t                 = {tau_t_prod}        ; time constant, in ps
ref_t                 = {T}      ; reference temperature in K

; Pressure coupling is on
pcoupl                  = {pcoupl_prod}     ; Pressure coupling on in NPT
pcoupltype              = {pcoupltype_prod}             ; uniform scaling of box vectors
tau_p                   = {tau_p_prod}                   ; time constant, in ps
ref_p                   = {ref_p_prod}                   ; reference pressure, in bar
compressibility = {alpha}      ; isothermal compressibility, bar^-1 for {solv}
;refcoord_scaling        = com

; Periodic boundary conditions
pbc                   = {pbc_prod}         ; 3-D PBC

; Dispersion correction
DispCorr              = {DispCorr_prod}    ; account for cut-off vdW scheme

; Velocity generation
gen_vel               = {gen_vel_prod} 
'''

with open('prod.mdp', 'w') as file:
    file.write(content)

content = 'copy MD files on cluster and run JOB_Prod.sh'
with open('step4_production', 'w') as file:
    file.write(content)

# make JOB script
content = f'''#!/bin/bash
#SBATCH --job-name={molname}_in_{solv}
#SBATCH --output=slurm-%J.out
#SBATCH --error=slurm-%J.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00-5:00:00
#SBATCH --partition=shared-cpu

module load GCC/11.3.0 CUDA/11.4.1 OpenMPI/4.1.4 GROMACS/2023.1

gmx grompp -f prod.mdp -c npt/npt.gro -p topol.top -o prod/prod
gmx mdrun -nt 32 -v -deffnm prod/prod'''

with open('JOB_Prod.sh', 'w') as file:
    file.write(content)

### Step 5 Analysis ### 
# make script to center trajectory on the molcule and generate index file
if not os.path.exists('analysis'):
    os.mkdir('analysis')

content = '''#!/bin/bash

# center trajectory on molecule
gmx trjconv -s prod/prod.tpr -f prod/prod.xtc -o prod/prod-mol.xtc -pbc mol -center 

# generate index file
gmx make_ndx -f prod/prod.gro '''

with open('step5_analysis.sh', 'w') as file:
    file.write(content) 

### Make the file executable ##
os.chmod('step0_make_box.sh', 0o755)
os.chmod('step1_energy_minimization.sh', 0o755)
os.chmod('step2_NVT_equilibration.sh', 0o755)
os.chmod('step3_NPT_equilibration.sh', 0o755)
os.chmod('step5_analysis.sh', 0o755)

### Make style file to plot trajectory ###
content = ''' # General display settings
display depthcue off
color Display Background white
axes location Off
display projection Orthographic

# Load the molecule
mol load gro prod/prod.gro xtc prod/prod-mol.xtc

# Delete all existing representations
set molid [molinfo top]
set nreps [molinfo $molid get numreps]

for {set i [expr {$nreps - 1}]} {$i >= 0} {incr i -1} {
    mol delrep $i $molid
}

# Representation for RES
mol representation Licorice
mol color type
mol material Opaque
mol selection {resname RES}
mol addrep top

# Solvent RES
mol representation points 2
mol color type
mol material Opaque
mol selection {resname RES}
mol addrep top'''

#with open('vmd_style.sh', 'w') as file:
    #file.write(content)

# generate a times file
content = f''' ### NVT Equilibration ###
total simulation time = {float(nsteps_nvt) * float(dt_nvt)} ps
time step = {1000 * float(dt_nvt)} fs
results saved every = {float(nstenergy_nvt) * float(dt_nvt)} ps
total number of frames = {float(nsteps_nvt) / float(nstenergy_nvt)}

### NPT Equilibration ###
total simulation time = {float(nsteps_nvt) * float(dt_npt)} ps
time step = {1000 * float(dt_npt)} fs
results saved every = {float(nstenergy_npt) * float(dt_npt)} ps
total number of frames = {float(nsteps_npt) / float(nstenergy_npt)}

### Production ###
total simulation time = {float(nsteps_prod) * float(dt_prod) / 1000} ns
time step = {1000 * float(dt_prod)} fs
results saved every = {float(nstenergy_prod) * float(dt_prod)} ps
total number of frames = {float(nsteps_prod) / float(nstenergy_prod)}
'''
with open('sim_time.txt', 'w') as file:
    file.write(content) 
