### README - Perfoming Basic Molecular Dynamics Simulations with Gromacs ###

1.) Using veloxchem generate .gro and .itp files for your molecule. See veloxchem folder 
for further info. To use velox chem conda environment: 
conda activate vlxchem 

You want to optimize the geometry of your molecule with gaussian, then calculate RESP charges 
with veloxchem. Then use the veloxchem pipelines to generate the force field files. 

If your molecule contains rotable bonds you also have to perform a relaxed scan calculation with Gaussian
then with velox chem you can reparametrize the force field (for this see reparam FF folder in example/molecular_dynamics/reparam_FF)

2.) Then you can look for your solvent in the GAFF force field to get its itp and gro file for detailed instructions on how 
to set up the topol.top file look at examples/molecular_dynamics/notes 

[ defaults ]
; nbfunc        comb-rule       gen-pairs        fudgeLJ   fudgeQQ
1               2               yes             0.500000  0.833333

[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon      
 c3       c3          0.00000  0.00000   A     3.39771e-01   4.51035e-01
...

#include "Py_Et_Py.itp"
#include "ACN.itp"

[ system ]
Py_Et_Py in ACN

[ molecules ]
; Compound        nmols
PYR               1
ACN               1434

should look something like this your topol file. In particular make sure no [defaults] or [atomtypes] in either itp file and no 
duplictes in [atomtypes]

3.) Run the and fill in the information in the gen_input_gromacs.py file. This will generate: 

make_box.sh 
run_energy_min.sh
run_nvt.sh
run_npt.sh 

which can be run via the termal for the making the box steps, energy min, and NVT and NPT equilibration

4.) The script also generate the necessary files for the production and the files needed for the server. 
