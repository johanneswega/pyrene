import os
import shutil
import subprocess

# ============================================================
# USER SETTINGS
# ============================================================

# Text file containing:
# time [ps]    rCOM [nm]
WINDOW_FILE = "times_found.txt"

# Pulling trajectory and topology
TRAJ = "pull/pull-mol.xtc"
TPR = "pull/pull.tpr"

# Index file
INDEX = "index.ndx"

# Output directory
OUTDIR = "umbrella"

# Files to copy into every umbrella directory
FILES_TO_COPY = [
    "index.ndx",
    "topol.top",
    "BPe.itp",
    "DCM.itp",
    "DPA.itp",
]

# Umbrella MDP template
MDP_TEMPLATE = "umbrella.mdp"


# ============================================================
# HELPER FUNCTION
# ============================================================

def run_command(command, cwd=None):
    print("\nRunning:")
    print(" ".join(command))

    result = subprocess.run(command, cwd=cwd)

    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with return code {result.returncode}"
        )


# ============================================================
# READ WINDOW DATA
# ============================================================

windows = []

with open(WINDOW_FILE, "r") as f:
    for line in f:
        line = line.strip()

        # Ignore blank lines and comments
        if not line or line.startswith("#"):
            continue

        parts = line.split()

        if len(parts) < 2:
            continue

        time_ps = float(parts[0])
        rcom_nm = float(parts[1])

        windows.append((time_ps, rcom_nm))


# Sort windows by rCOM
windows.sort(key=lambda x: x[1])


print("\nUmbrella windows found:")
for time_ps, rcom_nm in windows:
    print(f"  time = {time_ps:8.3f} ps    rCOM = {rcom_nm:.3f} nm")


# ============================================================
# CHECK INPUT FILES
# ============================================================

required_files = [
    WINDOW_FILE,
    TRAJ,
    TPR,
    INDEX,
    MDP_TEMPLATE,
]

required_files += FILES_TO_COPY

for filename in required_files:
    if not os.path.isfile(filename):
        raise FileNotFoundError(
            f"Required file not found: {filename}"
        )


# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# CREATE EACH UMBRELLA WINDOW
# ============================================================

for time_ps, rcom_nm in windows:

    # Example:
    # rCOM = 0.35 -> umbrella_035
    window_name = f"umbrella_{int(round(rcom_nm * 100)):03d}"

    window_dir = os.path.join(OUTDIR, window_name)

    print("\n" + "=" * 70)
    print(f"Setting up {window_name}")
    print(f"  time = {time_ps:.3f} ps")
    print(f"  target rCOM = {rcom_nm:.3f} nm")
    print("=" * 70)

    os.makedirs(window_dir, exist_ok=True)


    # --------------------------------------------------------
    # Copy common files
    # --------------------------------------------------------

    for filename in FILES_TO_COPY:

        destination = os.path.join(
            window_dir,
            os.path.basename(filename)
        )

        shutil.copy2(filename, destination)


    # --------------------------------------------------------
    # Create window-specific umbrella.mdp
    # --------------------------------------------------------

    with open(MDP_TEMPLATE, "r") as f:
        mdp = f.read()

    # Replace a placeholder:
    #
    #     TARGET_DISTANCE
    #
    # in umbrella.mdp with the actual rCOM.

    mdp = mdp.replace(
        "TARGET_DISTANCE",
        f"{rcom_nm:.6f}"
    )

    mdp_path = os.path.join(window_dir, "umbrella.mdp")

    with open(mdp_path, "w") as f:
        f.write(mdp)


    # --------------------------------------------------------
    # Extract starting structure from pulling trajectory
    # --------------------------------------------------------

    start_gro = os.path.join(window_dir, "start.gro")

    run_command([
        "bash",
        "-c",
        f'printf "2\\n0\\n" | gmx_mpi trjconv '
        f'-f "{TRAJ}" '
        f'-s "{TPR}" '
        f'-o "{start_gro}" '
        f'-dump "{time_ps:.6f}" '
        f'-center'
    ])


    # --------------------------------------------------------
    # Create SLURM job script
    # --------------------------------------------------------

    job_name = f"BPe_DPA_{int(round(rcom_nm * 100)):03d}"

    job_script = f"""#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --output=slurm-%J.out
#SBATCH --error=slurm-%J.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00-5:00:00
#SBATCH --partition=shared-cpu

module load GCC/11.3.0 CUDA/11.4.1 OpenMPI/4.1.4 GROMACS/2023.1

gmx grompp \\
    -f umbrella.mdp \\
    -c start.gro \\
    -p topol.top \\
    -n index.ndx \\
    -o umbrella.tpr

if [ $? -ne 0 ]; then
    echo "grompp failed."
    exit 1
fi

gmx mdrun \\
    -nt 32 \\
    -v \\
    -deffnm umbrella
"""

    job_path = os.path.join(window_dir, "umbrella_job.sh")

    with open(job_path, "w") as f:
        f.write(job_script)

    os.chmod(job_path, 0o755)

print("\n")
print("=" * 70)
print("All umbrella windows have been prepared.")
print("=" * 70)

print("\nDirectories:")
for time_ps, rcom_nm in windows:
    window_name = f"umbrella_{int(round(rcom_nm * 100)):03d}"
    print(f"  {OUTDIR}/{window_name}/")
