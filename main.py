"""
This is the main script. Select here which materials you want to calculate and what workflows to carry out.
"""

# external imports
import os
import pickle

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
ncores = 8

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

# material IDs
# (for local calculation only use one material at a time!)

# full material list
material_id = [
    "mp-66",
    "mp-111",
    "mp-149",
    "mp-463",
    "mp-661",
    "mp-682",
    "mp-804",
    "mp-856",
    "mp-909",
    "mp-1087",
    "mp-1138",
    "mp-1143",
    "mp-1253",
    "mp-1265",
    "mp-1342",
    "mp-1415",
    "mp-1500",
    "mp-1550",
    "mp-1672",
    "mp-1784",
    "mp-1985",
    "mp-2064",
    "mp-2160",
    "mp-2472",
    "mp-2490",
    "mp-2534",
    "mp-2605",
    "mp-2758",
    "mp-8062",
    "mp-20351",
    "mp-22862",
    "mp-22870",
    "mp-22875",
    "mp-22898",
    "mp-22903",
    "mp-22916",
    "mp-23037",
    "mp-23155",
    "mp-23167",
    "mp-23193",
    "mp-23202",
    "mp-23251",
    "mp-23268",
    "mp-23291",
    "mp-23703",
    "mp-28077",
    "mp-611517",
    "mp-612118",
    "mp-732337",
    "mp-995214",
    "mp-1228891",
]

"""
END OF USER INPUT SECTION
"""

# name of the calculation folder
calc_folder = "calc"

# get the base directory so each script knows where the source folder is
base_dir = os.getcwd()

# loop over all materials
for mp in material_id:
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
        structure, name, metal_check = basic_utils.get_structure(api_key, mp)

        # if the compound is a metal in DFT/the Materials Project, skip it
        if metal_check:
            continue

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
