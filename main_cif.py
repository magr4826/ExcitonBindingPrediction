"""
This is the main script. Select here which materials you want to calculate and what workflows to carry out.
"""

# external imports
import os
import pickle
from pymatgen.io.cif import CifParser

# local imports
from api_key import *
import src.utils.basic_utils as basic_utils

"""
START OF USER INPUT SECTION
"""

# calculation setup ("local" or "batchjob")
# (local is mostly just used for testing and debugging purposes)
calc_setup = "batchjob"

# number of cores per job
ncores = 16

# requested job memory (does nothing when doing local calculation)
memory = 40000  # (mb)

# calculation type defined through script names (see "/src/workflows")
# arrays of workflows that will be done one after each other
# (convergence workflows should be carried out before any other calculations, if accurate results are needed)
script_name = [
    "qe_convergence_SG15",
    "qe_ldos",
    "qe_wannier",
    "qe_effective_mass_indirect",
    "vasp_static_dielectric",
]

# script_name = ["qe_ldos_ranges"]

# material IDs
# (for local calculation only use one material at a time!)

# full material list
material_id = ["KRrBr2"]

"""
END OF USER INPUT SECTION
"""

# name of the calculation folder
calc_folder = "calc"

# get the base directory so each script knows where the source folder is
base_dir = os.getcwd()

# loop over all materials
for mp in material_id:
    name = mp
    # start in the correct directory
    os.chdir(base_dir)
    print(
        mp,
        flush=True,
    )
    # setup everything
    if not os.path.exists(os.path.join(calc_folder, mp)):
        # create directory
        os.makedirs(os.path.join(calc_folder, mp))

    if not os.path.exists(os.path.join(calc_folder, mp, "structure.pckl")):
        structure = CifParser(os.path.join("cifs", mp + ".cif")).get_structures()[0]

        # save the structure to a file which the workflows can read
        f = open(os.path.join(calc_folder, mp, "structure.pckl"), "wb")
        pickle.dump([structure, name], f)
        f.close()

    # create the job file and start the job
    filename = basic_utils.start_calc(
        base_dir,
        calc_setup,
        mp,
        script_name,
        file_name="exc_pred_job",
        ncores=ncores,
        memory=memory,
        calc_folder=calc_folder,
    )
