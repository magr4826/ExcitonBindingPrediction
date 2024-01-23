"""
Functions that write input files for QE (but don't actually run them)
"""

# external imports
import textwrap
import numpy as np

# local imports
import src.utils.basic_utils as basic_utils
import src.utils.qe_helper as qe_helper


def write_pw(
    calc_data,
    occupations="smearing",
    disk_io="medium",
    diagonalization="david",
    symmorphic=False,
    kpoints_mode="automatic",
    nosym=False,
):
    """
    Writes the input file for a pw.x calculation
    INPUT:
        calc_data:      calc_data of the calculation which should be run
        For the other inputs, see the pw.x documentation
    OUTPUT:
        output_filename:       Name of the written input file
    """

    output_filename = (
        "pw_"
        + calc_data.identifier
        + "_"
        + calc_data.calc_type
        + "_k"
        + str(calc_data.kppa)
        + "_E"
        + str(calc_data.pw_cutoff)
    )
    pseudo = {}
    for elem in calc_data.structure.types_of_species:
        pseudo[elem.name] = elem.name + ".upf"
    k_points_grid = calc_data.k_points_grid
    control = {
        "calculation": calc_data.calc_type,
        "prefix": calc_data.identifier,
        "outdir": f"out/",
        "pseudo_dir": calc_data.pseudo,
        "disk_io": disk_io,
    }
    system = {
        "ecutwfc": calc_data.pw_cutoff,
        "occupations": occupations,
        "degauss": basic_utils.ev2ha(0.025) / 2,  # 25meV smearing to "emulate" 300K
        "smearing": "gaussian",
        "nbnd": int(
            np.ceil(
                np.max(
                    [
                        qe_helper.qe_get_electrons(calc_data),
                        8,
                    ]
                )
            )
        ),
        "force_symmorphic": symmorphic,
        "nosym": nosym,
    }
    if calc_data.ibrav != 0:
        system.update(calc_data.ibrav)
    electrons = {
        "conv_thr": 1e-10,
        "diagonalization": diagonalization,
        "diago_full_acc": True,
        "diago_thr_init": 5e-6,
    }
    input = qe_helper.qe_PWInput(
        calc_data.structure,
        pseudo=pseudo,
        kpoints_grid=k_points_grid,
        control=control,
        system=system,
        electrons=electrons,
        kpoints_mode=kpoints_mode,
    )
    input.write_file(output_filename + ".in")
    return output_filename


def write_emass(calc_data, diagonalization="paro"):
    """
    Writes the input file for an emass calculation
    INPUT:
        calc_data:      calc_data of the calculation which should be run
        For the other inputs, see the pw.x documentation
    OUTPUT:
        output_filename:       Name of the written input file
    """

    output_filename = (
        "pw_"
        + calc_data.identifier
        + "_emass"
        + "_k"
        + str(calc_data.kppa)
        + "_E"
        + str(calc_data.pw_cutoff)
    )
    pseudo = {}
    for elem in calc_data.structure.types_of_species:
        pseudo[elem.name] = elem.name + ".upf"
    k_points_grid = calc_data.k_points_grid
    control = {
        "calculation": calc_data.calc_type,
        "prefix": calc_data.identifier,
        "outdir": f"out/",
        "pseudo_dir": calc_data.pseudo,
        "disk_io": "nowf",
    }
    system = {
        "ecutwfc": calc_data.pw_cutoff,
        "occupations": "smearing",
        "degauss": basic_utils.ev2ha(0.025) / 2,  # 25meV smearing to "emulate" 300K
        "smearing": "gaussian",
        "nbnd": int(
            np.ceil(
                np.max(
                    [
                        qe_helper.qe_get_electrons(calc_data),
                        8,
                    ]
                )
            )
        ),
    }
    if calc_data.ibrav != 0:
        system.update(calc_data.ibrav)
    electrons = {
        "conv_thr": 1e-10,
        "diagonalization": diagonalization,
        "diago_full_acc": True,
        "diago_thr_init": 5e-6,
    }
    input = qe_helper.qe_PWInput(
        calc_data.structure,
        pseudo=pseudo,
        kpoints_grid=k_points_grid,
        control=control,
        system=system,
        electrons=electrons,
        kpoints_mode="tpiba",
    )
    input.write_file(output_filename + ".in")

    return output_filename


def write_ldos(identifier, E1, E2):
    """
    Writes the input file for a ldos postprocessing calculation
    Specifially, the calculation will return a cube file with the number of states per volume lying between E1 and E2
    INPUT:
        identifier:     identifier/prefix of the material (usually the Materials Project id)
        E1:             Lower boundary for energy integration
        E2:             Upper boundary for energy integration
    OUTPUT:
        filename_ldos:       Name of the written input file
    """

    filename_ldos = f"{identifier}_ldos"
    file = open(filename_ldos + ".in", "w")
    file.write(
        textwrap.dedent(
            f"""\
        &inputpp
            prefix = '{identifier}'
            outdir='out/',
            plot_num=10,
            emin = {E1},
            emax = {E2},
            delta_e = 0.1,
        /
        &plot
            iflag=3,
            output_format=6,
            fileout='ldos.cube',
        /
    """
        )
    )
    file.close()
    return filename_ldos
