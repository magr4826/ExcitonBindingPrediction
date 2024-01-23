"""
Basic utilities which are generally useful in for HT calculations,
but not specific to any single code
"""

# external imports
import os
import textwrap
import warnings
from mp_api.client import MPRester
from pymatgen.io.vasp.inputs import Kpoints


# functions which convert units
def ha2ev(val):
    return val * 13.6057039763 * 2


def ev2ha(val):
    return val / 13.6057039763 / 2


def au2ang(val):
    return val * 0.529177249


def ang2au(val):
    return val / 0.529177249


def get_structure(api_key, material_id):
    """
    Gets structure from Materials Project with given material_id.
    INPUT:
        api_key:        Your Materials Project API key
        material_id:    The Materials Project id of the material
    OUTPUT:
        structure:      The structure of the material
        name:           The sum formula of the material in the Materials Project
        metal_check:    Whether the material is a metal or not
    """
    metal_check = False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with MPRester(api_key) as m:
            docs = m.summary.search(
                material_ids=[material_id],
                fields=["formula_pretty", "structure", "is_metal"],
            )
            name = docs[0].formula_pretty
            structure = docs[0].structure
            if docs[0].is_metal:
                print("Material is a metal!", flush=True)
                metal_check = True

    return structure, name, metal_check


def get_kpt_grid(structure, kppa):
    """
    Generates the k-point-grid for a structure with given k-point density
    INPUT:
        structure:      The structure for which the k-point-grid should be generated
        kppa:           The k-point density in points per inverse Angstrom^3. If equal to 0, return a gamma-only grid
    OUTPUT:
        kgrid:          A 3x1 list which specifies the number of subdivision in reciprocal space in each dimension
                        (i.e., the normal way to format k-grids). The k-grid is rounded up to the nearest even number
                        to ensure good parallelization
    """
    kpts = Kpoints.automatic_density(structure=structure, kppa=kppa, force_gamma=True)
    if kppa == 0:  # for gamma-only calculations
        kgrid = [1, 1, 1]
    else:
        kgrid = [i if i % 2 == 0 else i + 1 for i in kpts.kpts[0]]
    return kgrid


# YOU NEED TO ADJUST THE BATCHJOB PART OF THIS FUNCTION FOR THE SUPERCOMPUTER
# THAT YOU USE, ELSE THIS MAY NOT WORK... THE LOCAL SETUP SHOULD WORK.
def start_calc(
    base_dir,
    calc_setup,
    material_id,
    script_names,
    file_name,
    ncores,
    memory,
    calc_folder="calc",
):
    """
    Writes the .lsf/.bash file for a given workflow.
    INPUT:
        base_dir:       path to the main folder
        calc_setup:     local or batchjob calculation
        material_id:    Materials Project ID
        script_names:   arrays with names of workflows that will be processed in the given order (must be in workflow dictionary)
        file_name:      name of the .lsf/.bash file
        ncores:         number of cores requested
        memory:         maximum required memory
        calc_folder:    name of the calculation folder
    OUTPUT:
        file_name:      name of the .lsf/.bash file
    """
    # go to the calculation folder
    os.chdir(os.path.join(base_dir, calc_folder, material_id))

    # if it is a batchjob, write the BSUB-lsf file
    # (our cluster uses the IBM Spectrum LSF system)
    # CHANGE THIS ACCORDING TO YOUR OWN SETUP
    if calc_setup == "batchjob":
        f = open(file_name + ".lsf", "w")
        f.write(
            textwrap.dedent(
                f"""\
            #!/bin/csh
            #BSUB -q BatchXL
            #BSUB -R "mem > {memory} span[hosts=1]"
            #BSUB -oo lsf_qe_%J.log
            #BSUB -eo lsf_qe_%J.log
            #BSUB -J Exc_B
            #BSUB -n {ncores}
            date
            pwd
            """
            )
        )
        # we don't do 'conda activate base' here as it is already loaded in your terminal on the cluster
        f.write(f"setenv PYTHONPATH {base_dir}\n")
        f.write(f"echo {material_id}\n")
        for wf in script_names:
            f.write(f"echo {wf}\n")
            f.write(
                f"python3 {base_dir}/src/workflows/{wf}.py {base_dir} {ncores} {material_id}\n"
            )
        f.write("date")
        f.close()
        os.system(f"bsub < {file_name}.lsf")
    # if it is a local calculation, write a bash file
    elif calc_setup == "local":
        f = open(file_name + ".bash", "w")
        f.write(
            textwrap.dedent(
                f"""\
            #!/bin/bash
            date
            pwd
            """
            )
        )
        f.write("conda activate base\n")
        f.write(f"export PYTHONPATH={base_dir}\n")
        f.write(f"echo {material_id}\n")
        for wf in script_names:
            f.write(f"echo {wf}\n")
            f.write(
                f"python3 {base_dir}/src/workflows/{wf}.py {base_dir} {ncores} {material_id}\n"
            )
        f.write("date")
        f.close()
        os.system(f"bash {file_name}.bash")
    else:
        raise Exception(
            f"""QUITTING: No job submitted! Wrong calc_setup given, 
            only "local" and "batchjob" are currently supported!"""
        )
    return file_name
