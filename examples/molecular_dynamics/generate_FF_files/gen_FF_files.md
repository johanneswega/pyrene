# Generate Force Field Files for Your Molecule with VeloxChem

This guide describes how to generate **GROMACS topology and coordinate files** for a molecule using [VeloxChem](https://veloxchem.org/).

The main goal is to generate:

* `.gro` — the molecular coordinates/topology information used by GROMACS
* `.itp` — the molecular force-field parameters and topology
* `.top` — VeloxChem also generates a `.top` file, but in we will use the `.itp` file when assembling the complete system topology of actual project

The workflow is:

1. Optimize the molecular geometry with a quantum-chemistry program.
2. Export the optimized geometry as an `.xyz` file.
3. Inspect the optimized geometry.
4. Calculate **RESP atomic charges** with VeloxChem.
5. Generate the GROMACS force-field files.
6. If necessary, **reparameterize rotatable bonds**.

---

## 1. Start VeloxChem

First, activate the Conda environment containing VeloxChem:

```bash
conda activate vlxchem
```

Then start Jupyter Lab:

```bash
jupyter lab
```

Open the VeloxChem notebook that you will use for the force-field generation.

---

## 2. Optimize the Molecular Geometry

Before generating the force field, you need a reasonable **optimized molecular geometry**.

We first perform a geometry optimization using a quantum-chemistry program. In our workflow, we use **Gaussian** for this step.

---

## 3. Extract the Optimized Geometry in XYZ Format

VeloxChem requires the optimized geometry in **XYZ format**.

If your optimized geometry is currently stored in a Gaussian .log file, you can extract it using the following code:



```python
import cclib

# name of folder that contains gaussian .log file
mol = 'BPe_anion'

elements = [
    None, "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    # etc.
]

data = cclib.io.ccread(f"{mol}/{mol}.log")

coords = data.atomcoords[-1]

with open(f"{mol}/optimized.xyz", "w") as f:
    f.write(f"{data.natom}\n")
    f.write("Optimized geometry\n")

    for atomic_number, (x, y, z) in zip(data.atomnos, coords):
        element = elements[atomic_number]
        f.write(f"{element} {x:.8f} {y:.8f} {z:.8f}\n")
```

## 4. Inspect the Optimized Geometry

Before continuing with the force-field generation, it is useful to visually inspect the optimized structure.

We can use Veloxchem which uses **py3Dmol** to look at the structure. Also be sure to set the charge and multiplicity of the molecule correctly for the subsequent RESP charge calculation.


```python
import veloxchem as vlx
molecule = vlx.Molecule.read_xyz_file(f"{mol}/optimized.xyz")
molecule.set_charge(-1)         
molecule.set_multiplicity(2)    
molecule.show(atom_indices=True)
```


<div id="3dmolviewer_1789572738668471"  style="position: relative; width: 400px; height: 300px;">
        <p id="3dmolwarning_1789572738668471" style="background-color:#ffcccc;color:black">3Dmol.js failed to load for some reason.  Please check your browser console for error messages.<br></p>
        </div>
<script>

var loadScriptAsync = function(uri){
  return new Promise((resolve, reject) => {
    //this is to ignore the existence of requirejs amd
    var savedexports, savedmodule;
    if (typeof exports !== 'undefined') savedexports = exports;
    else exports = {}
    if (typeof module !== 'undefined') savedmodule = module;
    else module = {}

    var tag = document.createElement('script');
    tag.src = uri;
    tag.async = true;
    tag.onload = () => {
        exports = savedexports;
        module = savedmodule;
        resolve();
    };
  var firstScriptTag = document.getElementsByTagName('script')[0];
  firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
});
};

if(typeof $3Dmolpromise === 'undefined') {
$3Dmolpromise = null;
  $3Dmolpromise = loadScriptAsync('https://cdn.jsdelivr.net/npm/3dmol@2.5.4/build/3Dmol-min.js');
}

var viewer_1789572738668471 = null;
var warn = document.getElementById("3dmolwarning_1789572738668471");
if(warn) {
    warn.parentNode.removeChild(warn);
}
$3Dmolpromise.then(function() {
viewer_1789572738668471 = $3Dmol.createViewer(document.getElementById("3dmolviewer_1789572738668471"),{backgroundColor:"white"});
viewer_1789572738668471.zoomTo();
	viewer_1789572738668471.addModel("34\n\nC              0.323503000000         2.862637000000         0.000000000000\nC             -0.931462000000         3.535525000000         0.000000000000\nC             -2.120066000000         2.837175000000         0.000000000000\nC             -2.139444000000         1.427687000000         0.000000000000\nC             -0.900992000000         0.711825000000         0.000000000000\nC             -0.900992000000        -0.711826000000         0.000000000000\nC             -2.139444000000        -1.427687000000         0.000000000000\nC             -3.361798000000        -0.681332000000         0.000000000000\nC             -3.361798000000         0.681332000000         0.000000000000\nC             -2.120066000000        -2.837176000000         0.000000000000\nC             -0.931462000000        -3.535525000000         0.000000000000\nC              0.323504000000        -2.862637000000         0.000000000000\nC              0.341958000000        -1.427134000000         0.000000000000\nC              1.584928000000        -0.722181000000         0.000000000000\nC              2.786302000000        -1.484612000000         0.000000000000\nC              2.754970000000        -2.871450000000         0.000000000000\nC              1.550050000000        -3.568861000000         0.000000000000\nC              0.341958000000         1.427134000000         0.000000000000\nC              1.584928000000         0.722181000000         0.000000000000\nC              2.786302000000         1.484612000000         0.000000000000\nC              2.754970000000         2.871450000000         0.000000000000\nC              1.550050000000         3.568861000000         0.000000000000\nH             -0.933730000000         4.629504000000         0.000000000000\nH             -3.073528000000         3.375240000000         0.000000000000\nH             -4.306969000000        -1.233635000000         0.000000000000\nH             -4.306969000000         1.233635000000         0.000000000000\nH             -3.073528000000        -3.375240000000         0.000000000000\nH             -0.933730000000        -4.629504000000         0.000000000000\nH              3.752418000000        -0.980529000000         0.000000000000\nH              3.698877000000        -3.426261000000         0.000000000000\nH              1.535240000000        -4.662217000000         0.000000000000\nH              3.752418000000         0.980529000000         0.000000000000\nH              3.698877000000         3.426261000000         0.000000000000\nH              1.535239000000         4.662217000000         0.000000000000\n");
	viewer_1789572738668471.addLabel("1",{"position": {"x": 0.323503, "y": 2.8626369999999994, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("2",{"position": {"x": -0.9314619999999999, "y": 3.535525, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("3",{"position": {"x": -2.120066, "y": 2.837175, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("4",{"position": {"x": -2.139444, "y": 1.4276869999999997, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("5",{"position": {"x": -0.9009919999999999, "y": 0.7118249999999999, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("6",{"position": {"x": -0.9009919999999999, "y": -0.7118259999999998, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("7",{"position": {"x": -2.139444, "y": -1.4276869999999997, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("8",{"position": {"x": -3.3617979999999994, "y": -0.6813319999999999, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("9",{"position": {"x": -3.3617979999999994, "y": 0.6813319999999999, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("10",{"position": {"x": -2.120066, "y": -2.837176, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("11",{"position": {"x": -0.9314619999999999, "y": -3.535525, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("12",{"position": {"x": 0.323504, "y": -2.8626369999999994, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("13",{"position": {"x": 0.341958, "y": -1.427134, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("14",{"position": {"x": 1.584928, "y": -0.722181, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("15",{"position": {"x": 2.786302, "y": -1.484612, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("16",{"position": {"x": 2.75497, "y": -2.87145, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("17",{"position": {"x": 1.5500499999999997, "y": -3.5688609999999996, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("18",{"position": {"x": 0.341958, "y": 1.427134, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("19",{"position": {"x": 1.584928, "y": 0.722181, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("20",{"position": {"x": 2.786302, "y": 1.484612, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("21",{"position": {"x": 2.75497, "y": 2.87145, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("22",{"position": {"x": 1.5500499999999997, "y": 3.5688609999999996, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("23",{"position": {"x": -0.9337299999999998, "y": 4.629504, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("24",{"position": {"x": -3.073528, "y": 3.37524, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("25",{"position": {"x": -4.306969, "y": -1.2336349999999998, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("26",{"position": {"x": -4.306969, "y": 1.2336349999999998, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("27",{"position": {"x": -3.073528, "y": -3.37524, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("28",{"position": {"x": -0.9337299999999998, "y": -4.629504, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("29",{"position": {"x": 3.7524179999999996, "y": -0.9805289999999999, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("30",{"position": {"x": 3.698877, "y": -3.426261, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("31",{"position": {"x": 1.53524, "y": -4.662217, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("32",{"position": {"x": 3.7524179999999996, "y": 0.9805289999999999, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("33",{"position": {"x": 3.698877, "y": 3.426261, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.addLabel("34",{"position": {"x": 1.535239, "y": 4.662217, "z": 0.0}, "alignment": "center", "fontColor": 0, "backgroundColor": 16777215, "backgroundOpacity": 0.0, "fontSize": 16});
	viewer_1789572738668471.setViewStyle({"style": "outline", "width": 0.05});
	viewer_1789572738668471.setStyle({"stick": {}, "sphere": {"scale": 0.25}});
	viewer_1789572738668471.zoomTo();
viewer_1789572738668471.render();
});
</script>


## 5. Calculate RESP Atomic Charges

The next step is to calculate **RESP (Restrained Electrostatic Potential) charges** for the molecule.

### Why do we need partial atomic charges?

Molecular mechanics force fields represent electrostatic interactions using **partial atomic charges** assigned to individual atoms.

For example, in a typical molecular-dynamics simulation, the electrostatic interaction between two atoms is calculated from their partial charges. Therefore, assigning physically reasonable charges is essential for obtaining realistic intermolecular interactions.

### Why RESP charges?

RESP charges are obtained by fitting atomic partial charges to reproduce the molecule's **quantum-mechanical electrostatic potential (ESP)**.

The RESP procedure adds restraints to the charge-fitting process. These restraints help avoid unrealistic charge values and generally produce a more chemically reasonable and transferable charge distribution.

### RESP vs. CHELP/CHelpG

Another common approach is **CHELP** or **CHelpG** (Charges from Electrostatic Potentials using a Grid-based method), which also determines atomic charges by fitting them to the quantum-mechanical electrostatic potential.

The main difference is that **RESP introduces restraints during the fitting procedure**, whereas the traditional CHELP/CHelpG approach is primarily an unconstrained ESP fit with a particular grid definition.

This matters because a straightforward ESP fit can sometimes produce very large or otherwise undesirable charges, particularly for atoms in chemically equivalent environments or for atoms whose charges are poorly determined by the ESP.

RESP was developed to address some of these issues and to produce charge sets that are more suitable for molecular-mechanics force fields.

### RESP and GAFF

For **GAFF (General AMBER Force Field)** workflows, RESP-derived charges are commonly used and are part of the established Amber parameterization methodology.

Therefore, if your goal is to generate a GAFF-compatible force field for a small organic molecule, **RESP charges are generally the appropriate choice**.

> **In short:**
> **ESP → fit charges to reproduce the electrostatic potential**
> **RESP → ESP fitting + restraints to obtain more chemically reasonable charges**
> **GAFF → commonly used together with RESP/Amber-style charge derivation**

---

## 6. Calculate RESP Charges with VeloxChem

The RESP charges can be calculated directly using VeloxChem on your laptop. Depending on the size of your molecule and the level of theory used, this calculation can take some time.


```python
basis = vlx.MolecularBasis.read(molecule, "6-31G*", ostream=None)
if molecule.get_multiplicity()>1:
    scf_drv = vlx.ScfUnrestrictedDriver()
else:
    scf_drv = vlx.ScfRestrictedDriver()
scf_results = scf_drv.compute(molecule, basis)
resp_drv = vlx.RespChargesDriver()
resp_charges = resp_drv.compute(molecule, basis, scf_results)
```

                                                                                                                              
                                                Self Consistent Field Driver Setup                                            
                                               ====================================                                           
                                                                                                                              
                       Wave Function Model             : Spin-Unrestricted Hartree-Fock                                       
                       Initial Guess Model             : Superposition of Atomic Densities                                    
                       Convergence Accelerator         : Two Level Direct Inversion of Iterative Subspace                     
                       Max. Number of Iterations       : 50                                                                   
                       Max. Number of Error Vectors    : 10                                                                   
                       Convergence Threshold           : 1.0e-06                                                              
                       ERI Screening Threshold         : 1.0e-12                                                              
                       Linear Dependence Threshold     : 1.0e-06                                                              
                                                                                                                              
    * Info * Starting Reduced Basis SCF calculation...                                                                        
    * Info * ...done. SCF energy in reduced basis set: -839.859720829338 a.u. Time: 84.21 sec.                                
                                                                                                                              
                                                                                                                              
                   Iter. | Hartree-Fock Energy | Energy Change | Gradient Norm | Max. Gradient | Density Change               
                   --------------------------------------------------------------------------------------------               
                      1      -840.145227329775    0.0000000000      0.25025755      0.01491915      0.00000000                
                      2      -840.152779827787   -0.0075524980      0.05745526      0.00305159      0.08625317                
                      3      -840.153815498841   -0.0010356711      0.03775993      0.00245335      0.03724136                
                      4      -840.154637919131   -0.0008224203      0.03207002      0.00234949      0.02767490                
                      5      -840.155465467912   -0.0008275488      0.02945823      0.00203725      0.02857658                
                      6      -840.158609255986   -0.0031437881      0.01633039      0.00109502      0.15154322                
                      7      -840.159350637517   -0.0007413815      0.00845324      0.00054324      0.09324613                
                      8      -840.159380664484   -0.0000300270      0.00538831      0.00036358      0.01040481                
                      9      -840.159430709435   -0.0000500450      0.00344665      0.00022266      0.01910810                
                     10      -840.159444659401   -0.0000139500      0.00227353      0.00023146      0.00825157                
                     11      -840.159450966954   -0.0000063076      0.00151827      0.00011528      0.00516347                
                     12      -840.159455859172   -0.0000048922      0.00097986      0.00007191      0.00688423                
                     13      -840.159457452144   -0.0000015930      0.00045238      0.00003111      0.00365730                
                     14      -840.159457721596   -0.0000002695      0.00029821      0.00001584      0.00171315                
                     15      -840.159457837770   -0.0000001162      0.00013263      0.00000890      0.00107118                
                     16      -840.159457855171   -0.0000000174      0.00008586      0.00000494      0.00040605                
                     17      -840.159457860875   -0.0000000057      0.00004263      0.00000317      0.00017653                
                     18      -840.159457862547   -0.0000000017      0.00002403      0.00000253      0.00010301                
                     19      -840.159457863158   -0.0000000006      0.00001515      0.00000190      0.00005227                
                     20      -840.159457863513   -0.0000000004      0.00000961      0.00000133      0.00004604                
                     21      -840.159457863669   -0.0000000002      0.00000527      0.00000091      0.00003640                
                     22      -840.159457863735   -0.0000000001      0.00000312      0.00000044      0.00002951                
                     23      -840.159457863753   -0.0000000000      0.00000201      0.00000025      0.00001592                
                     24      -840.159457863758   -0.0000000000      0.00000143      0.00000014      0.00000491                
                     25      -840.159457863762   -0.0000000000      0.00000078      0.00000009      0.00000780                
                                                                                                                              
                   *** SCF converged in 25 iterations. Time: 960.70 sec.                                                      
                                                                                                                              
                   Spin-Unrestricted Hartree-Fock:                                                                            
                   -------------------------------                                                                            
                   Total Energy                       :     -840.1594578638 a.u.                                              
                   Electronic Energy                  :    -2431.2934503009 a.u.                                              
                   Nuclear Repulsion Energy           :     1591.1339924371 a.u.                                              
                   ------------------------------------                                                                       
                   Gradient Norm                      :        0.0000007770 a.u.                                              
                                                                                                                              
                                                                                                                              
                   Ground State Information                                                                                   
                   ------------------------                                                                                   
                   Charge of Molecule            : -1.0                                                                       
                   Multiplicity (2S+1)           :  2                                                                         
                   Magnetic Quantum Number (M_S) :  0.5                                                                       
                   Expectation value of S**2     :  1.5610                                                                    
                                                                                                                              
                                                                                                                              
                                                    RESP Charges Driver Setup                                                 
                                                   ===========================                                                
                                                                                                                              
                                             Number of Conformers         :  1                                                
                                             Number of Layers             :  4                                                
                                             Points per Square Angstrom   :  1.0                                              
                                             Total Number of Grid Points  :  1120                                             
                                                                                                                              
                                                                                                                              
                                                         First Stage Fit                                                      
                                                        -----------------                                                     
                                                                                                                              
                                             Restraint Strength           :  0.0005                                           
                                             Restrained Hydrogens         :  No                                               
                                             Max. Number of Iterations    :  50                                               
                                             Convergence Threshold (a.u.) :  1e-06                                            
                                                                                                                              
                                          *** Charge fitting converged in 13 iterations.                                      
                                                                                                                              
                                            No. | Atom |  Constraints | Charges (a.u.)                                        
                                           --------------------------------------------                                       
                                              1     C                      0.034601                                           
                                              2     C                     -0.234599                                           
                                              3     C                     -0.213211                                           
                                              4     C                     -0.007702                                           
                                              5     C                      0.011521                                           
                                              6     C                      0.011520                                           
                                              7     C                     -0.007702                                           
                                              8     C                     -0.216509                                           
                                              9     C                     -0.216509                                           
                                             10     C                     -0.213210                                           
                                             11     C                     -0.234598                                           
                                             12     C                      0.034601                                           
                                             13     C                      0.066231                                           
                                             14     C                     -0.048213                                           
                                             15     C                     -0.204709                                           
                                             16     C                     -0.203311                                           
                                             17     C                     -0.192738                                           
                                             18     C                      0.066232                                           
                                             19     C                     -0.048213                                           
                                             20     C                     -0.204711                                           
                                             21     C                     -0.203310                                           
                                             22     C                     -0.192740                                           
                                             23     H                      0.117203                                           
                                             24     H                      0.121140                                           
                                             25     H                      0.122212                                           
                                             26     H                      0.122212                                           
                                             27     H                      0.121140                                           
                                             28     H                      0.117203                                           
                                             29     H                      0.122857                                           
                                             30     H                      0.120026                                           
                                             31     H                      0.105201                                           
                                             32     H                      0.122857                                           
                                             33     H                      0.120026                                           
                                             34     H                      0.105201                                           
                                           --------------------------------------------                                       
                                                   Total Charge  : -1.000000                                                  
                                                                                                                              
                                                           Fit Quality                                                        
                                                          -------------                                                       
                                           Relative Root-Mean-Square Error  :  0.010121                                       
                                                                                                                              
                                             *** No refitting in second stage needed.                                         
                                                                                                                              
                                           Reference:                                                                         
                                           J. Phys. Chem. 1993, 97, 10269-10280.                                              
                                                                                                                              



```python
import numpy as np
print(resp_charges)
np.savetxt(f'{mol}/RESP_charges.txt', resp_charges)
print("Atom     RESP charge")
print(20 * "-")

for label, resp_charge in zip(molecule.get_labels(), resp_charges):
    print(f"{label :s} {resp_charge : 18.6f}")

print(20 * "-")

print(f"Total: {resp_charges.sum() : 13.6f}")
```

    [ 0.03460066 -0.23459908 -0.21321128 -0.00770227  0.01152052  0.01152034
     -0.00770199 -0.21650928 -0.21650888 -0.21320973 -0.23459814  0.03460064
      0.06623143 -0.04821295 -0.20470892 -0.20331122 -0.19273835  0.066232
     -0.04821302 -0.2047106  -0.20331043 -0.19273966  0.11720349  0.12113997
      0.12221195  0.12221185  0.12113995  0.11720313  0.12285744  0.12002628
      0.10520127  0.12285741  0.12002624  0.1052012 ]
    Atom     RESP charge
    --------------------
    C           0.034601
    C          -0.234599
    C          -0.213211
    C          -0.007702
    C           0.011521
    C           0.011520
    C          -0.007702
    C          -0.216509
    C          -0.216509
    C          -0.213210
    C          -0.234598
    C           0.034601
    C           0.066231
    C          -0.048213
    C          -0.204709
    C          -0.203311
    C          -0.192738
    C           0.066232
    C          -0.048213
    C          -0.204711
    C          -0.203310
    C          -0.192740
    H           0.117203
    H           0.121140
    H           0.122212
    H           0.122212
    H           0.121140
    H           0.117203
    H           0.122857
    H           0.120026
    H           0.105201
    H           0.122857
    H           0.120026
    H           0.105201
    --------------------
    Total:     -1.000000


## 7. Generate the Force-Field Files

Once the optimized geometry and RESP charges have been obtained, we can use the **VeloxChem force-field pipeline** to generate the GROMACS-compatible files:


```python
RES_NAME = 'BPE'
ff_gen = vlx.MMForceFieldGenerator()
ff_gen.partial_charges = np.loadtxt(f'{mol}/RESP_charges.txt')
ff_gen.create_topology(molecule)
print(ff_gen.rotatable_bonds)

ff_gen.write_gromacs_files(f"{mol}/{RES_NAME}", f"{mol}/{RES_NAME}")
```

    * Info * Sum of partial charges is not a whole number.                                                                    
    * Info * Compensating by removing -1.000e-06 from the largest charge.                                                     
                                                                                                                              
    * Info * Using GAFF (v2.11) parameters.                                                                                   
             Reference: J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman, D. A. Case, J. Comput. Chem. 2004,
             25, 1157-1174.
                                                                                                                              
    []


The pipeline generates the files required to describe your molecule in a molecular-dynamics simulation.

The important output files are:

### `.gro`

The `.gro` file contains the molecular coordinates in a format that can be read by GROMACS.

### `.itp`

The `.itp` file contains the molecular topology and force-field information, including things such as:

* atom definitions,
* atom types,
* partial charges,
* bonds,
* angles,
* dihedrals,
* and other interaction parameters.

---

## 8. Reparameterize Rotatable Bonds if Necessary

There is one important additional consideration: **rotatable bonds**.

If your molecule contains bonds that can freely rotate, the automatically generated force-field parameters may not always reproduce the correct rotational energetics.

### What does "reparameterization" mean?

Force-field parameters describe the energetic cost of different molecular configurations.

For a rotatable bond, the relevant parameter is typically the **dihedral potential**, which determines how the energy changes as the bond rotates.

A generic dihedral potential can be thought of as:

```text
Energy
  ^
  |       /\          /\
  |      /  \        /  \
  |_____/    \______/    \____> dihedral angle
```

The default force field provides a predefined dihedral potential. However, for some molecules, this potential may not accurately describe the rotational energy obtained from quantum mechanics.

**Reparameterization** means adjusting the force-field parameters—typically the torsional/dihedral parameters—so that the molecular-mechanics energy profile better reproduces a quantum-mechanical reference.

This is particularly important when a rotatable bond controls:

* the molecule's preferred conformation,
* the relative stability of different conformers,
* or the overall flexibility of the molecule.

Reparameterization is a somewhat more advanced, but VeloxChem provides tools to make it relatively straightforward.

For an example of how to perform this procedure, see:

```text
examples/molecular_dynamics/reparam_FF
```
