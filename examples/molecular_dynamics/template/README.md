# README — Performing Basic Molecular Dynamics Simulations with GROMACS

This README describes the basic workflow for setting up and running molecular dynamics (MD) simulations with **GROMACS**, using **VeloxChem** and **Gaussian** to generate and parameterize the molecular force field.

---

## 1. Generate the Molecular Force Field

Use **VeloxChem** to generate the `.gro` and `.itp` files for your molecule. See the `examples/molecular_dynamics/generate_FF_files` folder for further information and examples.

### Recommended workflow

1. **Optimize the molecular geometry with Gaussian.**
2. **Calculate RESP charges with VeloxChem.**
3. Use the **VeloxChem pipelines** to generate the required force-field files (`.gro` and `.itp`).

### Molecules with rotatable bonds

If your molecule contains **rotatable bonds**, you should also perform a **relaxed potential-energy scan** with Gaussian.

The resulting scan can then be used with VeloxChem to **reparameterize the force field**.

For an example, see:

```text
examples/molecular_dynamics/reparam_FF
```

---

## 2. Set Up the Solvent and Topology

Next, identify your solvent in the **GAFF force field** and obtain the corresponding `.itp` and `.gro` files or generate them yourself using the veloxchem workflow above.

For detailed instructions on constructing the `topol.top` file, see:

```text
examples/molecular_dynamics/notes/
```

A typical `topol.top` file should contain sections similar to the following:

```text
[ defaults ]
; nbfunc        comb-rule       gen-pairs        fudgeLJ   fudgeQQ
1               2               yes              0.500000  0.833333

[ atomtypes ]
; name   bond_type     mass     charge   ptype   sigma         epsilon
 c3      c3            0.00000  0.00000  A       3.39771e-01   4.51035e-01
...

#include "Molecule.itp"
#include "ACN.itp"

[ system ]
Molecule in ACN

[ molecules ]
; Compound        nmols
MOL               1
ACN               1434
```

### Important checks

In particular, make sure that:

* Neither of the included `.itp` files contains a `[ defaults ]` section.
* Neither of the included `.itp` files contains an `[ atomtypes ]` section if these definitions are already present in `topol.top`.
* There are **no duplicate atom-type definitions** in `[ atomtypes ]`.
* The molecule names in `[ molecules ]` match the names defined in the corresponding `.itp` files.
* The number of solvent molecules corresponds to the actual contents of the simulation box.

---

## 3. Generate the GROMACS Simulation Scripts

Fill in the required simulation parameters in:

```text
gen_input_gromacs.py
```

Running this script generates the following files:

```text
step0_make_box.sh
step1_energy_minimization.sh
step2_NVT_equilibration.sh
step3_NPT_equilibration.sh
step4_production
step5_analysis.sh
```

### Running the simulations

The first stages can be run directly from the terminal:

1. **Box construction**
2. **Energy minimization**
3. **NVT equilibration**
4. **NPT equilibration**

The generated scripts automate these steps using GROMACS.

You will likely need to **manually adapt `step0_make_box.sh`** if you are simulating a system more complicated than a single chromophore dissolved in a solvent.

### Production run

For the production simulation, the script also generates the necessary **JOB file** for submitting the calculation to the cluster.

Before submitting the job, check that the generated parameters and resource requirements are appropriate for your cluster setup.

Finally, there is also a `sim_time.txt`that summarizes all the simulation time used.

---

## 4. Visualize the Trajectories

For interactive visualization of trajectories and simulation boxes, I recommend **NGLView**. See:

```text
examples/molecular_dynamics/ngl/
```

for further information.

---

## 5. Analyze the Simulations

Examples of different analysis routines can be found in:

```text
examples/molecular_dynamics/analysis/
```
