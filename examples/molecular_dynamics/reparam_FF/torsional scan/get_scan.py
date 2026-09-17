import re

# ============================================================
# USER SETTINGS
# ============================================================

INPUT_XYZ = "scan.xyz"
INPUT_LOG = "Py_Et_Py.log"

OUTPUT_XYZ = "1-2-3-4.xyz"
OUTPUT_ENERGY = "scan_energy.dat"

# Dihedral atoms
DIHEDRAL = (1, 2, 3, 4)

# Starting dihedral angle
START_DIHEDRAL = -48.836

# Dihedral increment
DIHEDRAL_STEP = 20.0

# Hartree -> kcal/mol
HARTREE_TO_KCAL = 627.509474


# ============================================================
# FUNCTIONS
# ============================================================

def parse_scan_xyz(filename):
    """
    Read the XYZ scan file.

    Expected format:

    56
     Scan point 1 scf energy: -1307.61027891
    C ...
    ...
    H ...

    Returns a list containing scan point, energy and coordinates.
    """

    with open(filename, "r") as f:
        lines = f.readlines()

    scans = []
    i = 0

    while i < len(lines):

        if not lines[i].strip():
            i += 1
            continue

        # Number of atoms
        try:
            natoms = int(lines[i].strip())
        except ValueError:
            i += 1
            continue

        if i + 1 >= len(lines):
            break

        # Scan point and energy
        match = re.search(
            r"Scan point\s+(\d+)\s+scf energy:\s+([-\d.Ee+]+)",
            lines[i + 1],
            re.IGNORECASE
        )

        if not match:
            i += 1
            continue

        point = int(match.group(1))
        energy = float(match.group(2))

        # Coordinates
        atoms = lines[i + 2:i + 2 + natoms]

        if len(atoms) != natoms:
            raise RuntimeError(
                f"Could not read {natoms} atoms "
                f"for scan point {point}"
            )

        scans.append({
            "point": point,
            "energy": energy,
            "atoms": atoms
        })

        i += natoms + 2

    return scans


def parse_iterations(filename):
    """
    Extract the final optimization iteration for each
    Gaussian scan point.

    Looks for lines such as:

    Step number   6 out of a maximum of 336 on scan point 15 out of 19
    """

    iterations = {}

    pattern = re.compile(
        r"Step number\s+(\d+).*?"
        r"scan point\s+(\d+)\s+out of\s+(\d+)",
        re.IGNORECASE
    )

    with open(filename, "r", errors="ignore") as f:

        for line in f:

            match = pattern.search(line)

            if match:

                iteration = int(match.group(1))
                scan_point = int(match.group(2))

                # Keep the last iteration encountered
                iterations[scan_point] = iteration

    return iterations


def calculate_dihedral(point):
    """
    Calculate the dihedral angle corresponding to
    a scan point.
    """

    angle = START_DIHEDRAL + (point - 1) * DIHEDRAL_STEP

    return angle 


# ============================================================
# MAIN
# ============================================================

def main():

    d1, d2, d3, d4 = DIHEDRAL

    print("Reading XYZ file...")
    scans = parse_scan_xyz(INPUT_XYZ)

    if not scans:
        raise RuntimeError(
            "No scan points were found in the XYZ file."
        )

    print(f"Found {len(scans)} scan points.")

    print("Reading Gaussian log file...")
    iterations = parse_iterations(INPUT_LOG)

    print(
        f"Dihedral: {d1}-{d2}-{d3}-{d4}"
    )

    print(
        f"Starting angle: {START_DIHEDRAL:.2f} degrees"
    )

    print(
        f"Step: {DIHEDRAL_STEP:.2f} degrees"
    )

    print()


    # ========================================================
    # WRITE FORMATTED XYZ
    # ========================================================

    with open(OUTPUT_XYZ, "w") as xyz_out:

        total_points = len(scans)

        for scan in scans:

            point = scan["point"]
            energy = scan["energy"]

            dihedral = calculate_dihedral(point)

            iteration = iterations.get(point, "?")

            xyz_out.write(
                f"{len(scan['atoms'])}\n"
            )

            xyz_out.write(
                f"Scan Cycle {point}/{total_points} ; "
                f"Dihedral {d1}-{d2}-{d3}-{d4} = "
                f"{dihedral:.2f} ; "
                f"Iteration {iteration} "
                f"Energy {energy:.8f}\n"
            )

            for atom in scan["atoms"]:
                xyz_out.write(atom)


    # ========================================================
    # CONVERT ENERGIES TO KCAL/MOL
    # ========================================================

    energies_kcal = []

    for scan in scans:

        energy_hartree = scan["energy"]

        energy_kcal = (
            energy_hartree * HARTREE_TO_KCAL
        )

        energies_kcal.append(energy_kcal)


    # ========================================================
    # CALCULATE RELATIVE ENERGIES
    # ========================================================

    minimum_energy = min(energies_kcal)

    relative_energies = [
        energy - minimum_energy
        for energy in energies_kcal
    ]


    # ========================================================
    # WRITE ENERGY FILE
    # ========================================================

    with open(OUTPUT_ENERGY, "w") as energy_out:

        energy_out.write(
            f"# Dihedral {d1}-{d2}-{d3}-{d4}\n"
        )

        energy_out.write(
            "# Columns: Dihedral(deg) "
            "Energy(Hartree) "
            "Energy(kcal/mol) "
            "Relative Energy(kcal/mol)\n"
        )

        for i, scan in enumerate(scans):

            point = scan["point"]

            dihedral = calculate_dihedral(point)

            energy_hartree = scan["energy"]

            energy_kcal = energies_kcal[i]

            relative_energy = relative_energies[i]

            energy_out.write(
                f"{dihedral:10.2f} "
                f"{energy_hartree:18.10f} "
                f"{energy_kcal:18.8f} "
                f"{relative_energy:18.8f}\n"
            )


    # ========================================================
    # SUMMARY
    # ========================================================

    print("Done!")
    print()
    print(f"Formatted XYZ : {OUTPUT_XYZ}")
    print(f"Energy file   : {OUTPUT_ENERGY}")
    print()
    print(
        f"Minimum energy: "
        f"{minimum_energy:.8f} kcal/mol"
    )


if __name__ == "__main__":
    main()