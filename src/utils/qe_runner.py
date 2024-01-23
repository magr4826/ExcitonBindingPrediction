"""
Functions that can start pw.x calculations
"""

# external imports
import os
import pickle
import shutil
import textwrap
import numpy as np
import xml.etree.ElementTree as ET
from pymatgen.io.common import VolumetricData
import itertools

# local imports
from src.utils.calc_data_class import calc_data
import src.utils.qe_write as qe_write
import src.utils.qe_helper as qe_helper
from src.utils.constants import *
import src.utils.basic_utils as basic_utils


def qe_pw_run(
    calc_data,
    input_name,
    ncores,
    kwargs={},
):
    """
    This function writes and starts a pw.x calculation and does some rudimentary error handling.
    INPUT:
        calc_data:      The calculation that should be run
        input_name:     Which function should be used to write the input (usually a function in qe_write)
        ncores:         How many cores should be used by the calculation
        kwargs:         Keyword arguments to be passed to the function input_name
    OUTPUT:
        filename:       name of the input file which was executed
    """

    # write the input file for pw.x and start the calculation
    filename = input_name(calc_data, **kwargs)
    os.system(f"mpirun -np {ncores} pw.x -inp {filename}.in > {filename}.out")

    # if there is an error, try to use the more stable paro diagonalization
    # (errors can happen with nscf calculations, we never observed one in scf calculations)
    with open(filename + ".out") as f:
        if "Error" in f.read():
            kwargs["diagonalization"] = "paro"
            filename = input_name(calc_data, **kwargs)
            os.system(f"mpirun -np {ncores} pw.x -inp {filename}.in > {filename}.out")
            open("paro.txt", "a").close()  # create a file so we know that paro is used

    # check if the "eigenvalues not converged" appear
    # but we need to separate two cases
    # 1st case: last iteration of a scf calculation
    # 2nd case: nscf calculation
    with open(filename + ".out") as f:
        # parse the output file
        if "nscf" in filename:
            out_str = f.read()  # keep the whole output file
        else:
            out_str = f.read()  # only we the part after the last iteration
            out_str = out_str.split("iteration #")[
                -1
            ]  # we only keep the part after the last iteration

        # check if the eigenvalues did not converge
        if "eigenvalues not converged" in out_str:
            kwargs["diagonalization"] = "paro"
            filename = input_name(calc_data, **kwargs)
            os.system(f"mpirun -np {ncores} pw.x -inp {filename}.in > {filename}.out")
            open("paro.txt", "a").close()  # create a file so we know that paro is used

    return filename


def qe_convergence_checker(id, wf_dir, ncores, conv_dir="qe_convergence_SG15"):
    """
    Check if convergence for a material has been performed
    If yes:     Copy output of converged scf calculation to workflow dir
                Return converged convergence parameters
    If no:      Perform one scf with default parameters in workflow dir
                Return default (unconverged) convergence parameters
    This function should be called at the start of most QE workflows
    INPUT:
        id:             Materials Project id of the material
        wf_dir:         Directory of the workflow which called this function
        ncores:         Number of cores used (in case convergence hasnt been run before)
        conv_dir:       Name of the workflow directory in which the convergence was run (in case there are various versions)
    OUTPUT:
        conv_flag:      Boolean, whether convergence ran before or not
        calc_data_conv: calc_data of the converged calculation, if it ran before, or the default parameters otherwise
    """

    # check if convergence has been succesfully carried out
    conv_flag = False
    conv_dir = os.path.join(wf_dir, os.pardir, conv_dir)
    if os.path.exists(conv_dir):
        if os.path.exists(os.path.join(conv_dir, "conv_calc_data.pckl")):
            conv_flag = True
            print(
                f"\nConvergence has been carried out, using converged parameters and results from\n{conv_dir}",
                flush=True,
            )
        else:
            print(
                f"\nConvergence has started but not successfully finished, starting calculation with initial values for\n{conv_dir}",
                flush=True,
            )
    else:
        print(
            f"\nConvergence has not been started. Starting calculation with initial values for\n{conv_dir}",
            flush=True,
        )

    if conv_flag:
        # load the converged calc-data object
        f = open(os.path.join(conv_dir, "conv_calc_data.pckl"), "rb")
        calc_data_conv = pickle.load(f)
        f.close()

        # create the necessary folders
        if not os.path.exists(os.path.join(wf_dir, "out")):
            os.mkdir(os.path.join(wf_dir, "out"))
        if not os.path.exists(os.path.join(wf_dir, "out", id + ".save")):
            os.mkdir(os.path.join(wf_dir, "out", id + ".save"))

        # copy the results of the scf-calculation
        shutil.copy(
            os.path.join(conv_dir, "out", id + ".xml"),
            os.path.join(wf_dir, "out", id + ".xml"),
        )  # not sure if this is even necessary
        shutil.copy(
            os.path.join(conv_dir, "out", id + ".save", "data-file-schema.xml"),
            os.path.join(wf_dir, "out", id + ".save", "data-file-schema.xml"),
        )
        shutil.copy(
            os.path.join(conv_dir, "out", id + ".save", "charge-density.dat"),
            os.path.join(wf_dir, "out", id + ".save", "charge-density.dat"),
        )
        return conv_flag, calc_data_conv

    else:
        base_dir = os.path.join(wf_dir, os.pardir, os.pardir, os.pardir)
        structure, name, ibrav = qe_helper.qe_init_structure(id, base_dir)
        calc_data_conv = calc_data(
            structure,
            name,
            id=id,
            ibrav=ibrav,
            calc_type="scf",
            kppa=1500,
            pseudo=os.path.join(base_dir, "pseudo", "SG15"),
        )

        qe_pw_run(calc_data_conv, qe_write.write_pw, ncores, kwargs={})
        return conv_flag, calc_data_conv


def qe_emasses_run(calc_data, kpt0, ncores):
    """
    Run a calculation to get the effective masses around a certain k-point
    INPUT:
        calc_data:      calc_data which stores the information about the material
        kpt0:           k-point around which the effective mass should be calculated
        ncores:         Number of cores which should be used
    OUTPUT:
        output:         Array of occupations, effective masses (in x,y,z-direction), and energies at kpt0
    """

    # define the finite difference star and run the calculation
    k_points = textwrap.dedent(
        f"""\
        7
        {kpt0[0]} {kpt0[1]} {kpt0[2]} 1
        {kpt0[0]+0.01} {kpt0[1]} {kpt0[2]} 1
        {kpt0[0]-0.01} {kpt0[1]} {kpt0[2]} 1 
        {kpt0[0]} {kpt0[1]+0.01} {kpt0[2]} 1
        {kpt0[0]} {kpt0[1]-0.01} {kpt0[2]} 1 
        {kpt0[0]} {kpt0[1]} {kpt0[2]+0.01} 1
        {kpt0[0]} {kpt0[1]} {kpt0[2]-0.01} 1 
        """
    )
    calc_data.k_points_grid = k_points
    _ = qe_pw_run(calc_data, qe_write.write_emass, ncores=ncores)

    # postprocessing to get the effective masses using finite difference

    # load and parse the energies, occupations, k-points and alat (kpts is normalized in units of 2 pi/alat)
    path = os.path.join("out", calc_data.identifier + ".save", "data-file-schema.xml")
    tree = ET.parse(path)
    energies = tree.findall(".//{*}eigenvalues")
    kpts = tree.findall(".//{*}k_point")
    kpts = kpts[0 : len(energies)]
    occ = tree.findall(".//{*}occupations")
    occ = occ[1:]
    alat = tree.find(".//{*}atomic_structure")
    alat = float(alat.attrib["alat"]) * au
    for i in range(len(kpts)):
        kpts[i] = [float(kpt) for kpt in kpts[i].text.split()]
        energies[i] = [float(en) for en in energies[i].text.split()]
        occ[i] = [float(oc) for oc in occ[i].text.split()]
    occ = occ[0]
    emasses = []

    # perform the actual finite difference calculation
    for i in range(len(occ)):
        emass_x = (
            Ha
            / hbar**2
            * (energies[1][i] - 2 * energies[0][i] + energies[2][i])
            / (((kpts[1][0] - kpts[0][0]) * 2 * np.pi / alat) ** 2)
        ) ** (-1) / me
        emass_y = (
            Ha
            / hbar**2
            * (energies[3][i] - 2 * energies[0][i] + energies[4][i])
            / (((kpts[1][0] - kpts[0][0]) * 2 * np.pi / alat) ** 2)
        ) ** (-1) / me
        emass_z = (
            Ha
            / hbar**2
            * (energies[5][i] - 2 * energies[0][i] + energies[6][i])
            / (((kpts[1][0] - kpts[0][0]) * 2 * np.pi / alat) ** 2)
        ) ** (-1) / me
        emasses.append([emass_x, emass_y, emass_z])

    # output results
    emasses = np.array(emasses)
    energies_kpt0 = np.array_split(energies[0], len(energies[0]))
    energies_kpt0 = [[basic_utils.ha2ev(energy[0])] for energy in energies_kpt0]
    emasses = np.append(emasses, energies_kpt0, axis=1)
    output = np.array([occ, emasses[:, 0], emasses[:, 1], emasses[:, 2], emasses[:, 3]])
    output = np.transpose(output)

    return output


def qe_ldos_run(calc_data, E1, E2):
    """
    Postprocesses a calculation to return the number of states around the ions in a structure
    INPUT:
        calc_data:      calc_data which stores the information about the structure
        E1:             Lower boundary for energy integration (in eV)
        E2:             Upper boundary for energy integration (in eV)
    OUTPUT:
        ldos:           Array of local number of states
        structure:      structure read from the output cube file (This should be the same as in the calc_data)
    """
    # write and run ldos calculation
    filename_ldos = qe_write.write_ldos(calc_data.identifier, E1, E2)
    os.system(f"pp.x -inp {filename_ldos}.in > {filename_ldos}.out")

    # parse ldos output
    ldos = []
    radius = 1  # radius of the sphere in which the number of states is summed up
    ildos = VolumetricData.from_cube("ldos.cube")
    structure = ildos.structure
    a = ildos.dim
    coords = []
    for x, y, z in itertools.product(*(list(range(i)) for i in a)):
        coords.append([x / a[0], y / a[1], z / a[2]])
    for site in structure:
        sites_dist = structure.lattice.get_points_in_sphere(coords, site.coords, radius)
        data = np.array(sites_dist, dtype=object)
        sum = 0
        for point in data:
            fcoords = point[0]
            supercell = point[3]
            fcoords = fcoords - supercell
            sum += ildos.value_at(fcoords[0], fcoords[1], fcoords[2])
        sum = sum * structure.volume / np.prod(a) * (basic_utils.ang2au(1) ** 3)
        ldos.append(sum)

    return ldos, structure
