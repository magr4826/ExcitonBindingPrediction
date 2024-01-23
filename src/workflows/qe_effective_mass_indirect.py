"""
This workflow calculates the effective masses at the VBM and the CBM
in the cartesian directions. It also outputs estimates for the direct 
and the indirect gap, as it is the most accurate workflow.
"""

# external imports
import sys
import os
import time
import numpy as np
import xml.etree.ElementTree as ET
import shutil

# local imports
from src.utils.constants import *
import src.utils.qe_helper as qe_helper
import src.utils.qe_write as qe_write
import src.utils.qe_runner as qe_runner

# start timing
start_time = time.time()

# base directory
base_dir = str(sys.argv[1])

# number of cores
ncores = int(sys.argv[2])

# get the id from the args that are called with this script
id = str(sys.argv[3])

# workflow directory
wf_dir = os.path.join(os.getcwd(), "qe_effective_masses_indirect")
if not os.path.exists(wf_dir):
    os.mkdir(wf_dir)
os.chdir(wf_dir)

# initialize structure for a QuantumEspresso calculation
structure, name, ibrav = qe_helper.qe_init_structure(id, base_dir)

# check if convergence has been carried out
convergence_flag, calc_data_conv = qe_runner.qe_convergence_checker(id, wf_dir, ncores)

# the effective_mass calculation needs increased accuracy,
# so we increase the number of plane waves and k-points
# and run an nscf calculation
calc_data_conv.k_points_grid = [2 * kpt for kpt in calc_data_conv.k_points_grid]
calc_data_conv.pw_cutoff = 2 * calc_data_conv.pw_cutoff
calc_data_conv.calc_type = "nscf"
filename = qe_runner.qe_pw_run(calc_data_conv, qe_write.write_pw, ncores=ncores)

# find the k-points of indirect gap
xml_path = os.path.join("out", id + ".xml")
num_elec = qe_helper.qe_get_electrons(calc_data_conv)
(
    kpoints,
    eigenvalues,
    _,
    _,
) = qe_helper.qe_get_eigenvalues(xml_path, num_elec)

# this is the most accurate calculation, so we will also get the gaps from here
_, indirect_gap, direct_gap = qe_helper.qe_identify_manifolds(xml_path, tolerance=0.1)

# get the k-point of the VBM
valence_band = eigenvalues[:, int(num_elec / 2 - 1)]
vbm_index = np.argmax(valence_band)
vbm_k = kpoints[vbm_index]
vbm_k = [float(num) for num in vbm_k.split()]

# get the k-point of the CBM
conductance_band = eigenvalues[:, int(num_elec / 2)]
cbm_index = np.argmin(conductance_band)
cbm_k = kpoints[cbm_index]
cbm_k = [float(num) for num in cbm_k.split()]

# get the masses, occupations and energies
hole_mass_array = qe_runner.qe_emasses_run(calc_data_conv, vbm_k, ncores)
elec_mass_array = qe_runner.qe_emasses_run(calc_data_conv, cbm_k, ncores)

# save the masses
np.savetxt(
    os.path.join(f"hole_masses_{id}.csv"),
    hole_mass_array,
    header="occ,m_x,m_y,m_z,energy_gamma",
    comments="",
)

np.savetxt(
    os.path.join(f"elec_masses_{id}.csv"),
    elec_mass_array,
    header="occ,m_x,m_y,m_z,energy_gamma",
    comments="",
)

# Save the gaps
np.savetxt(f"indirect_gap_{id}.csv", [indirect_gap])
np.savetxt(f"direct_gap_{id}.csv", [direct_gap])

# delete the output folder at the end
if os.path.isdir("out"):
    shutil.rmtree("out")

# save calculation time
with open("timing.txt", "a+") as f:
    f.write(
        f"{os.path.basename(__file__):<25}  {(time.time() - start_time):7.2f} s  {ncores} cores\n"
    )
