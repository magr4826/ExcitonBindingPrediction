"""
This workflow calculates the local density of states (i.e., states, per energy, per volume)
and from those calculates the number of states localized at each atomic site.
"""

# external imports
import sys
import os
import time
from pymatgen.io.common import VolumetricData
import shutil
import itertools

# local imports
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
wf_dir = os.path.join(os.getcwd(), "qe_ldos_ranges")
if not os.path.exists(wf_dir):
    os.mkdir(wf_dir)
os.chdir(wf_dir)

# initialize structure for a QuantumEspresso calculation
structure, name, ibrav = qe_helper.qe_init_structure(id, base_dir)

# check if convergence has been carried out
convergence_flag, calc_data_conv = qe_runner.qe_convergence_checker(
    id, wf_dir, ncores, "qe_convergence_SG15"
)

# create nscf calculation based on scf calculation and run it
# double the plane-wave cutoff to improve the FFT grid resolution
# and thus, the real space resolution
calc_data_conv.pw_cutoff = calc_data_conv.pw_cutoff * 2
# calc_data_conv.k_points_grid = [4 * kpt for kpt in calc_data_conv.k_points_grid]
calc_data_conv.calc_type = "nscf"
filename_nscf = qe_runner.qe_pw_run(
    calc_data_conv,
    qe_write.write_pw,
    ncores=ncores,
    kwargs={"disk_io": "low"},
)

# read out VBM and CBM
num_elec = int(qe_helper.qe_get_electrons(calc_data_conv))
xml_path = os.path.join(os.getcwd(), "out", id + ".xml")
_, _, vbm, cbm = qe_helper.qe_get_eigenvalues(xml_path, num_elec)

energy_range = [0.8, 1, 1.2]
real_space_range = [0.8, 1, 1.2]

for id_en, id_r in itertools.product(
    range(len(energy_range)), range(len(real_space_range))
):
    energy = energy_range[id_en]
    radius = real_space_range[id_r]
    # write, run and parse the ldos of the below and above gap states
    vbm_ldos, structure = qe_runner.qe_ldos_run(
        calc_data_conv, vbm - energy, vbm, radius
    )
    cbm_ldos, structure = qe_runner.qe_ldos_run(
        calc_data_conv, cbm, cbm + energy, radius
    )

    # save the atomic site and the corresponding number of localized above-/below-gap states
    with open(f"{id}_ldos_E{id_en}_r{id_r}.csv", "w") as outfile:
        for idx in range(len(structure.labels)):
            outfile.write(
                " ".join(
                    [
                        structure.labels[idx],
                        str(vbm_ldos[idx]),
                        str(cbm_ldos[idx]),
                        "\n",
                    ]
                )
            )

# delete the output folder at the end
if os.path.isdir("out"):
    shutil.rmtree("out")
os.remove(os.path.join(wf_dir, "tmp.pp"))
os.remove(os.path.join(wf_dir, "ldos.cube"))

# save calculation time
with open("timing.txt", "a+") as f:
    f.write(
        f"{os.path.basename(__file__):<25}  {(time.time() - start_time):7.2f} s  {ncores} cores\n"
    )
