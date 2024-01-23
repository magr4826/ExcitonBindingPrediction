# external imports
import os
import textwrap
from pymatgen.io.vasp.inputs import Kpoints, Poscar

# local import
import src.utils.vasp_helper as vasp_helper


def write_static_eps(structure, ncores, kppa=3000, isym0=False):
    # INCAR
    incar_str = textwrap.dedent(
        f"""\
        PREC     = Accurate 
        ENCUT    = 520 
        ALGO     = Normal 
        ISMEAR   = 0 
        SIGMA    = 0.01 
        EDIFF    = 1E-8 
        NELM     = 100 
        LREAL    = .FALSE. 
        NPAR     = {ncores:d} 
        KPAR     = 1 
        LEPSILON = .TRUE. 
        LPEAD    = .TRUE.
        IBRION   = 8 
        LWAVE    = .FALSE. 
        SYMPREC  = 1E-8 
        """
    )
    if isym0:
        incar_str = incar_str + "ISYM     = 0"
    with open("INCAR", "w") as f:
        f.write(incar_str)

    # KPOINTS
    kpt = Kpoints.automatic_density(structure=structure, kppa=kppa, force_gamma=True)
    kpt.kpts[0] = [i if i % 2 == 0 else i + 1 for i in kpt.kpts[0]]
    kpt.write_file("KPOINTS")

    # POSCAR (suppress the pymatgen warnings)
    poscar = Poscar(structure=structure)
    poscar.write_file("POSCAR")

    # POTCAR
    element_paths = [
        f"/usr/app-soft1/vasp/potpaw_PBE.54/{vasp_helper.potcar_dict[elem]:s}/POTCAR"
        for elem in poscar.site_symbols
    ]
    os.system(f"cat {' '.join(element_paths):s} > POTCAR")
