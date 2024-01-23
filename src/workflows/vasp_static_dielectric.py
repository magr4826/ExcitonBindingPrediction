# external imports
import os
import sys
import time
import warnings
import numpy as np
from pymatgen.io.vasp.outputs import Vasprun

# local imports
import src.utils.vasp_helper as vasp_helper
import src.utils.vasp_write as vasp_write

# start timing
start_time = time.time()

# base directory
base_dir = str(sys.argv[1])

# number of cores
ncores = int(sys.argv[2])

# get the id from the args that are called with this script
id = str(sys.argv[3])

# workflow directory
wf_dir = os.path.join(os.getcwd(), "vasp_static_dielectric")
if not os.path.exists(wf_dir):
    os.mkdir(wf_dir)
os.chdir(wf_dir)

# initialize structure for a vasp calculation
structure, name = vasp_helper.vasp_init_structure(id, base_dir)

# create the input
vasp_write.write_static_eps(structure, ncores, kppa=3000)

# run the calculation
os.system(f"mpirun -np {ncores} vasp_std > log.file")

with open("log.file") as f:
    if "generating k-point grid is not commensurate" in f.read():
        vasp_write.write_static_eps(structure, ncores, kppa=3000, isym0=True)
        os.system(f"mpirun -np {ncores} vasp_std > log.file")


# parse the output
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    vasprun = Vasprun("vasprun.xml")
    eps_static = vasprun.epsilon_static
    eps_ionic = vasprun.epsilon_ionic

# delete unwanted files ...
for f in os.listdir():
    if f not in [
        "INCAR",
        "KPOINTS",
        "POSCAR",
        "POTCAR",
        "OUTCAR",
        "vasprun.xml",
        "log.file",
    ]:
        os.remove(f)

# save them as numpy arrays
np.savetxt(f"eps_static_{id}.csv", eps_static)
np.savetxt(f"eps_ionic_{id}.csv", eps_ionic)

# save calculation time
with open("timing.txt", "a+") as f:
    f.write(
        f"{os.path.basename(__file__):<25}  {(time.time() - start_time):7.2f} s  {ncores} cores\n"
    )
