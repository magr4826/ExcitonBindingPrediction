"""
Functions that are
a) Only necessary for QE and 
b) Don't start calculations on their own (i.e. don't need qe_write)
"""

# external imports
import os
import re
import math
import spglib
import pickle
import numpy as np
import xml.etree.ElementTree as ET
from ase.io.espresso import ibrav_to_cell
from pymatgen.io.pwscf import PWInput
from pymatgen.core import Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import PeriodicSite, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher

# local imports
from src.utils.basic_utils import ang2au, ha2ev, ev2ha


class qe_PWInput(PWInput):
    """
    Initializes a PWSCF input file.
    Adapted version of the original pymatgen version by Miguel Marques (thanks again!)
    See original pymatgen code for documentation
    """

    def __str__(self):
        out = []
        site_descriptions = {}

        if self.pseudo is not None:
            site_descriptions = self.pseudo
        else:
            c = 1
            for site in self.structure:
                name = None
                for k, v in site_descriptions.items():
                    if site.properties == v:
                        name = k

                if name is None:
                    name = site.specie.symbol + str(c)
                    site_descriptions[name] = site.properties
                    c += 1

        def to_str(v):
            if isinstance(v, str):
                return f"'{v}'"
            if isinstance(v, float):
                return f"{str(v).replace('e', 'd')}"
            if isinstance(v, bool):
                if v:
                    return ".TRUE."
                return ".FALSE."
            return v

        for k1 in ["control", "system", "electrons", "ions", "cell"]:
            v1 = self.sections[k1]
            out.append(f"&{k1.upper()}")
            sub = []
            for k2 in sorted(v1.keys()):
                if isinstance(v1[k2], list):
                    n = 1
                    for l in v1[k2][: len(site_descriptions)]:
                        sub.append(f"  {k2}({n}) = {to_str(v1[k2][n - 1])}")
                        n += 1
                else:
                    sub.append(f"  {k2} = {to_str(v1[k2])}")
            if k1 == "system":
                if "ibrav" not in self.sections[k1]:
                    sub.append("  ibrav = 0")
                if "nat" not in self.sections[k1]:
                    sub.append(f"  nat = {len(self.structure)}")
                if "ntyp" not in self.sections[k1]:
                    sub.append(f"  ntyp = {len(site_descriptions)}")
            sub.append("/")
            out.append(",\n".join(sub))

        out.append("ATOMIC_SPECIES")
        for k, v in sorted(site_descriptions.items(), key=lambda i: i[0]):
            e = re.match(r"[A-Z][a-z]?", k).group(0)
            if self.pseudo is not None:
                p = v
            else:
                p = v["pseudo"]
            out.append(f"  {k}  {Element(e).atomic_mass:.4f} {p}")

        out.append("ATOMIC_POSITIONS crystal")
        if self.pseudo is not None:
            for site in self.structure:
                out.append(f"  {site.specie} {site.a:.8f} {site.b:.8f} {site.c:.8f}")
        else:
            for site in self.structure:
                name = None
                for k, v in sorted(site_descriptions.items(), key=lambda i: i[0]):
                    if v == site.properties:
                        name = k
                out.append(f"  {name} {site.a:.8f} {site.b:.8f} {site.c:.8f}")

        out.append(f"K_POINTS {self.kpoints_mode}")
        if self.kpoints_mode == "automatic":
            kpt_str = [f"{i}" for i in self.kpoints_grid]
            kpt_str.extend([f"{i}" for i in self.kpoints_shift])
            out.append(f"  {' '.join(kpt_str)}")
        elif (
            self.kpoints_mode == "crystal"
        ):  # just used for wannier calculation in our case
            out.append(self.kpoints_grid)
        elif (
            self.kpoints_mode == "tpiba"
        ):  # just used for effective mass calculation in our case
            out.append(self.kpoints_grid)
        elif self.kpoints_mode == "gamma":
            pass

        # Difference to original
        if "ibrav" not in self.sections["system"]:
            out.append("CELL_PARAMETERS angstrom")
            for vec in self.structure.lattice.matrix:
                out.append(f"  {vec[0]:.15f} {vec[1]:.15f} {vec[2]:.15f}")
        return "\n".join(out) + "\n\n"


def qe_standardize_cell(structure, symprec=0.1):
    """
    Obtain the standard primitive structure from the input structure.
    INPUT:
        structure:      Structure that is supposed to be processed
        symprec:        Precision for symmetry prediction. Lower the value to get more precision
                        at a higher risk of not identifying symmetries
    OUTPUT:
        standard_structure: Standardized primitive structure
    """
    # Atomic positions have to be specified by scaled positions for spglib.
    lattice = structure.lattice.matrix
    scaled_positions = structure.frac_coords
    numbers = [i.specie.Z for i in structure.sites]
    cell = (lattice, scaled_positions, numbers)
    lattice, scaled_positions, numbers = spglib.standardize_cell(
        cell, to_primitive=True, symprec=symprec
    )
    s = Structure(lattice, numbers, scaled_positions)
    standard_structure = s.get_sorted_structure()
    return standard_structure


def qe_get_ibrav(structure):
    """
    Transforms a structure into the QuantumEspresso standard format to make sure that the symmetry is detected correctly.
    INPUT:
        structure:      Structure that is analyzed
    OUTPUT:
        espresso_in:    dict containing the celldm and ibrav
        qe_structure:   Standardized structure
    """
    # first we determine the conventional structure
    sym = SpacegroupAnalyzer(structure, symprec=1e-5)
    std_struct = sym.get_conventional_standard_structure(international_monoclinic=False)

    # obtain the space group number
    spg = sym.get_space_group_number()

    # this is the structure in espresso input format
    espresso_in = {}
    espresso_in["celldm(1)"] = ang2au(std_struct.lattice.a)

    if spg in [
        195,
        198,
        200,
        201,
        205,
        207,
        208,
        212,
        213,
        215,
        218,
        221,
        222,
        223,
        224,
    ]:
        # simple cubic
        espresso_in["ibrav"] = 1

    elif spg in [196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228]:
        # face-centered cubic
        espresso_in["ibrav"] = 2

    elif spg in [197, 199, 204, 206, 211, 214, 217, 220, 229, 230]:
        # body-centered cubic
        espresso_in["ibrav"] = 3

    elif ((spg >= 168) and (spg <= 194)) or spg in [
        143,
        144,
        145,
        147,
        149,
        150,
        151,
        152,
        153,
        154,
        156,
        157,
        158,
        159,
        162,
        163,
        164,
        165,
    ]:
        # hexagonal
        espresso_in["ibrav"] = 4
        espresso_in["celldm(3)"] = std_struct.lattice.c / std_struct.lattice.a

    elif spg in [146, 148, 155, 160, 161, 166, 167]:
        # rhombohedral
        espresso_in["ibrav"] = 5
        aR = math.sqrt(std_struct.lattice.a**2 / 3 + std_struct.lattice.c**2 / 9)
        cosgamma = (2 * std_struct.lattice.c**2 - 3 * std_struct.lattice.a**2) / (
            2 * std_struct.lattice.c**2 + 6 * std_struct.lattice.a**2
        )
        espresso_in["celldm(1)"] = ang2au(aR)
        espresso_in["celldm(4)"] = cosgamma

    elif spg in [
        75,
        76,
        77,
        78,
        81,
        83,
        84,
        85,
        86,
        89,
        90,
        91,
        92,
        93,
        94,
        95,
        96,
        99,
        100,
        101,
        102,
        103,
        104,
        105,
        106,
        111,
        112,
        113,
        114,
        115,
        116,
        117,
        118,
        123,
        124,
        125,
        126,
        127,
        128,
        129,
        130,
        131,
        132,
        133,
        134,
        135,
        136,
        137,
        138,
    ]:
        # simple tetragonal
        espresso_in["ibrav"] = 6
        espresso_in["celldm(3)"] = std_struct.lattice.c / std_struct.lattice.a

    elif spg in [
        79,
        80,
        82,
        87,
        88,
        97,
        98,
        107,
        108,
        109,
        110,
        119,
        120,
        121,
        122,
        139,
        140,
        141,
        142,
    ]:
        # body-centered tetragonal
        espresso_in["ibrav"] = 7
        espresso_in["celldm(3)"] = std_struct.lattice.c / std_struct.lattice.a

    elif spg in [
        16,
        17,
        18,
        19,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        47,
        48,
        49,
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57,
        58,
        59,
        60,
        61,
        62,
    ]:
        # simple orthorhombic
        espresso_in["ibrav"] = 8
        espresso_in["celldm(2)"] = std_struct.lattice.b / std_struct.lattice.a
        espresso_in["celldm(3)"] = std_struct.lattice.c / std_struct.lattice.a

    elif spg in [20, 21, 35, 36, 37, 38, 39, 40, 41, 63, 64, 65, 66, 67, 68]:
        # base-centered orthorhombic
        espresso_in["ibrav"] = 9
        espresso_in["celldm(2)"] = std_struct.lattice.b / std_struct.lattice.a
        espresso_in["celldm(3)"] = std_struct.lattice.c / std_struct.lattice.a

    elif spg in [22, 42, 43, 69, 70]:
        # face-centered orthorhombic
        espresso_in["ibrav"] = 10
        espresso_in["celldm(2)"] = std_struct.lattice.b / std_struct.lattice.a
        espresso_in["celldm(3)"] = std_struct.lattice.c / std_struct.lattice.a

    elif spg in [23, 24, 44, 45, 46, 71, 72, 73, 74]:
        # body-centered orthorhombic
        espresso_in["ibrav"] = 11
        espresso_in["celldm(2)"] = std_struct.lattice.b / std_struct.lattice.a
        espresso_in["celldm(3)"] = std_struct.lattice.c / std_struct.lattice.a

    else:
        raise NotImplementedError(f"ibrav not defined for spg {spg}")

    # get espresso unit cell lattice
    cell = ibrav_to_cell(espresso_in)[-1]

    # get sites for the new lattice
    new_sites = []
    latt = Lattice(cell)
    for s in std_struct:
        new_s = PeriodicSite(
            s.specie,
            s.coords,
            latt,
            to_unit_cell=True,
            coords_are_cartesian=True,
            properties=s.properties,
        )
        # I set the tolerance to the position tolerance of the periodic site which is 1e-5.
        # The default tolerance for the is_periodic_image function is 1e-8 which sometimes
        # does not recognize that a set of sites is periodic.
        if not any(
            [
                new_s.is_periodic_image(ns, tolerance=new_s.position_atol)
                for ns in new_sites
            ]
        ):
            new_sites.append(new_s)

    # generate output structure
    qe_structure = Structure.from_sites(
        new_sites,
        to_unit_cell=True,
        validate_proximity=True,  # validate_proximity=True should be set to be save ...
    )

    # check if the new and old primitive cell still match
    matcher = StructureMatcher(primitive_cell=False)
    if matcher.fit(structure, qe_structure):
        return espresso_in, qe_structure
    else:
        raise Exception(
            f"""
                QUITTING: The input primitive cell and the new primitive cell for Quantum Espresso dont match! 
                The detection of periodic sites may have gone wrong! The calculation will now run with ibrav=0
                which may slow the calculation, as less symmetries may be detected!
                """
        )


def qe_get_cutoff(structure, pseudo):
    """
    Parse the suggested plane wave cutoffs from the pseudopotentials associated with a structure and find the largest one.
    This is useful for getting an initial guess for convergence.
    INPUT:
        structure:      structure for which the suggested cutoff is found
        pseudo:         Right now a dummy, in a later version one can specify the pseudopotential here
    OUTPUT:
        pw_cutoff:      Suggested plane wave cutoff in Ry
    """
    # get the suggested cutoff for all element in the structure
    pw_cutoff = []
    for elem in structure.types_of_species:
        pw_cutoff.append(
            60
        )  # the SG15 pseudopotentials are optimized so that every element has the same cutoff
        # in the future one can specify the cutoffs for other pseudopotential here, like the PseudoDojo ones
        # so this function is thus just a dummy for now

    # find the maximum cutoff
    pw_cutoff = np.max(pw_cutoff)

    return pw_cutoff


def qe_init_structure(id, base_dir):
    """
    Initialize a structure stored earlier in a .pckl for further calculations
    INPUT:
        id:         Materials Project id of the material
        base_dir:   base_dir of the entire project
    OUTPUT:
        structure:  Loaded and standardized structure
        name:       Sum formula of the structure
        ibrav:      dict containing celldm and ibrav
    """

    # set path, get structure and name of system automatically
    f = open(os.path.join(base_dir, "calc", id, "structure.pckl"), "rb")
    structure, name = pickle.load(f)
    f.close()
    structure = qe_standardize_cell(structure, symprec=0.1)
    try:
        ibrav, structure = qe_get_ibrav(structure)
    except Exception:
        ibrav = 0

    return structure, name, ibrav


def qe_get_electrons(calc_data):
    """
    Gets the total number of included pseudopotential electrons in a structure saved in a calc_data.
    INPUT:
        calc_data:      Calc_data containing the structure
    OUTPUT:
        electrons:      Number of electrons
    """
    electrons = 0
    for elem in calc_data.structure.species:
        with open(os.path.join(calc_data.pseudo, elem.name + ".upf")) as file:
            for line in file:
                if "z_valence" in line:
                    electrons += float(
                        re.findall(
                            "-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[+-]?\ *[0-9]+)?", line
                        )[0]
                    )
                    break
    return electrons


def qe_read_totenergy(id):
    """
    Read out total energy from the results of a QuantumEspresso calculation.
    INPUT:
        id:             Materials Project id
    OUTPUT:
        etot:           Total energy in Hartree
    """
    path = os.path.join("out", id)
    tree = ET.parse(f"{path}.xml")
    etot = tree.find(".//{*}etot").text
    etot = float(etot)
    return etot


def qe_get_eigenvalues(xml_path, num_elec):
    """
    Parses the eigenvalues as a array of size (num_kpt x num_bands).
    Also obtains the vbm and cbm.
    INPUT:
        xml_path:       Path of the xml-output of a pw.x calculation
        num_elec:       Number of electrons in the cell
    OUTPUT:
        kpoints:        list of kpoints
        eigenvalues:    array of eigenvalues
        vbm:            VBM energy in eV
        cbm:            CBM energy in eV
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()
    kpoints = []
    for kpt in root.iter("k_point"):
        kpoints.append(kpt.text)
    eigenvalues = []
    for eigs in root.iter("eigenvalues"):
        eigs_num = [float(num) for num in eigs.text.strip("\n").split()]
        eigenvalues.append(eigs_num)
    eigenvalues = np.array(eigenvalues)

    vbm = ha2ev(max(eigenvalues[:, int(num_elec / 2 - 1)]))
    cbm = ha2ev(min(eigenvalues[:, int(num_elec / 2)]))

    return kpoints, eigenvalues, vbm, cbm


def qe_identify_manifolds(xml_path, tolerance):
    """
    Identify separate manifolds in the valence band and returns an
    array of band indices after which the band structure could be cut
    INPUT:
        xml_path:       Path of the xml-output of a pw.x calculation
        tolerance:      Energy in eV of separation between bands before they are considered separate
    OUTPUT:
        manifold_ids:   Band indices of top band in each found manifold
        indirect_gap:   Indirect gap in eV
        direct_gap:     Direct gap in eV
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()
    kpoints = []
    weights = []
    for kpt in root.iter("k_point"):
        kpoints.append(kpt.text)
        weights.append(float(kpt.attrib["weight"]))
    eigenvalues = []
    for eigs in root.iter("eigenvalues"):
        eigs_num = [float(num) for num in eigs.text.strip("\n").split()]
        eigenvalues.append(eigs_num)
    num_elec = float(root.findall(".//{*}nelec")[0].text)
    eigenvalues = np.array(eigenvalues)

    manifold_ids = []
    for band_idx in range(int(num_elec / 2) - 1):
        truth = np.greater_equal.outer(
            eigenvalues[:, band_idx + 1],
            eigenvalues[:, band_idx] + ev2ha(tolerance),
        )
        truth = np.all(truth)
        if truth:
            manifold_ids.append(band_idx + 1)

    HOMO = ha2ev(max(eigenvalues[:, int(num_elec / 2 - 1)]))
    LUMO = ha2ev(min(eigenvalues[:, int(num_elec / 2)]))
    indirect_gap = LUMO - HOMO
    direct_gap = np.min(
        ha2ev(eigenvalues[:, int(num_elec / 2)] - eigenvalues[:, int(num_elec / 2 - 1)])
    )
    return manifold_ids, indirect_gap, direct_gap
