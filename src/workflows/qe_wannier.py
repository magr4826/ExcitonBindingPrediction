"""
This workflow calculates the MLWFs. From these, it calculates the average irreducible spread,
the mean distance of the MLWF centers to the closest ion, and the derivative of the valence band
"""

# external imports
import sys
import os
import re
import time
import textwrap
import numpy as np
import shutil
import pandas as pd

# local imports
from src.utils.calc_data_class import calc_data
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
wf_dir = os.path.join(os.getcwd(), "qe_wannier")

if not os.path.exists(wf_dir):
    os.mkdir(wf_dir)
os.chdir(wf_dir)

# initialize structure for a QuantumEspresso calculation
structure, name, ibrav = qe_helper.qe_init_structure(id, base_dir)

# check if convergence has been carried out
convergence_flag, calc_data_conv_gga = qe_runner.qe_convergence_checker(
    id, wf_dir, ncores
)

# do you want the plot data for the mlwfs? THIS TAKES A LONG TIME!!!
plot = False

# for how many iterations should the wannierization run?
num_iter = (
    500  # low number for faster code testing, increase the number to get better MLWFs
)

# do you want to delete the output folder at the end?
delete_out = True

# Dimensions of the k-grid as string, which must be put in the wannier90 input file
mp_size = " ".join(str(k) for k in calc_data_conv_gga.k_points_grid)

# get k-meshes for pw-nscf and wannier calculation
kmesh_wannier = os.popen(f"kmesh.pl {mp_size} wann").read()
kmesh_pw = os.popen(f"kmesh.pl {mp_size}").read()
kmesh_pw = kmesh_pw[17:]
calc_data_wann = calc_data_conv_gga
calc_data_wann.calc_type = "nscf"
calc_data_wann.k_points_grid = kmesh_pw
filename_nscf = qe_runner.qe_pw_run(
    calc_data_wann,
    qe_write.write_pw,
    ncores=ncores,
    kwargs={"disk_io": "low", "kpoints_mode": "crystal", "nosym": True},
)

# write the .pw2wan file
# Use the SCDM method for the starting point
f = open(id + "_pw2wan.in", "w")
f.write(
    textwrap.dedent(
        f"""\
    &inputpp
    outdir = 'out/'
    prefix = '{id}'
    seedname = '{id}'
    spin_component = 'none'
    write_mmn = .true.
    write_amn = .true.
    write_unk = {plot}
    scdm_proj = .true.
    scdm_entanglement = 'isolated'
    /"""
    )
)
f.close()

# identify valence manifolds
# this assumes there is always at least one band that has to be excluded
# If all valence bands are connected, this can lead to an issue
manifolds, _, _ = qe_helper.qe_identify_manifolds(f"out/{id}.xml", tolerance=0.1)
max_manifold = 0
semicore_bands = ""
if not manifolds == []:
    max_manifold = manifolds[-1]
    semicore_bands = f"1-{max_manifold}"

# generate the string of atomic sites that will be written in the Wannier90-Input
atoms_frac = ""
for atom in structure.sites:
    atoms_frac += atom.species_string + " "
    atoms_frac += " ".join(str(coord) for coord in atom.frac_coords)
    atoms_frac += "\n"
unit_cell_cart = str(structure.lattice)

# calculate the number of MLWFs (all valence bands - disconnected valence bands)
nbnd = int(qe_helper.qe_get_electrons(calc_data_wann) * 0.5)
num_wann = int(nbnd - max_manifold)

# write the _win file (the input file for the wannier90 calculation)
f = open(id + ".win", "w")
f.write(
    textwrap.dedent(
        f"""\
num_wann = {num_wann}
exclude_bands = {semicore_bands} {int(nbnd)+1}-{int(np.max([2*nbnd,8]))}
num_iter = {num_iter}

auto_projections =.true.

begin atoms_frac
{atoms_frac}end atoms_frac

geninterp = true
geninterp_alsofirstder = true
        
begin unit_cell_cart
{unit_cell_cart}
end unit_cell_cart
        
mp_grid : {mp_size}
        
begin kpoints
{kmesh_wannier}end kpoints

"""
    )
)
f.close()

# execute the wannierization
os.system(f"mpirun -np {ncores} wannier90.x -pp {id}")
os.system(f"mpirun -np {ncores} pw2wannier90.x -inp {id}_pw2wan.in > {id}_pw2wan.out")
os.system(f"mpirun -np {ncores} wannier90.x {id}")

# check if the wannierization worked properly
# some jobs may just fail at random...
# in these cases we throw an exception
# than the user just needs to restart the workflow
# for the materials in question
# we noticed two possible reasons for failure

# 1st possible problem: the schur decomposition fails
werr_flag = False
for f in os.listdir():
    if f.endswith(".werr"):
        werr_flag = True
if werr_flag:
    raise Exception("\nWANNIERIZATION FAILED (Schur error)...")

# 2st possible problem: the wannierization does not convergence
with open(id + ".wout", "r") as f:
    wout_str = f.read()
omega_d0 = re.search(r"Omega OD[ \t]+=[ \t]+\*+", wout_str)
if omega_d0:
    raise Exception("\nWANNIERIZATION FAILED (convergence problem)...")

# read and save mean gauge-invariant spread
omega_I_M = []
wout = open(id + ".wout", "r")
wout_lines = wout.readlines()
wout.close()
for i, line in enumerate(wout_lines):
    if "Omega I" in line:
        omega_I_curr = float(line[48:56])
        omega_I_M.append(omega_I_curr / num_wann)
np.savetxt(f"omega_irr_{id}.csv", omega_I_M)

# read the mlwf centers and calculate their distance from the nearest atom, then save these
wout = open(id + ".wout", "r")
wout_lines = wout.readlines()
wout.close()
mlwf_center_cart = []
mlwf_spread = []
regex = re.compile("\(.*\n")
flag = False
mlwf_idx = 0
for i, line in enumerate(wout_lines):
    if "Final State" in line:
        flag = True
        continue
    if flag == True:
        mlwf_idx += 1
        mlwf_center_str = regex.findall(line)[0]
        mlwf_center_str = mlwf_center_str.replace("(", "")
        mlwf_center_str = mlwf_center_str.replace(")", "")
        mlwf_center_str = mlwf_center_str.replace(",", "")
        mlwf_center_str.strip()
        mlwf_center_str = mlwf_center_str.split()
        mlwf_center_curr = [float(num) for num in mlwf_center_str[0:3]]
        mlwf_spread_curr = float(mlwf_center_str[3])
        mlwf_center_cart.append(mlwf_center_curr)
        mlwf_spread.append(mlwf_spread_curr)
        if mlwf_idx >= num_wann:
            flag = False

lat = structure.lattice
mlwf_center_frac = lat.get_fractional_coords(mlwf_center_cart)
distances = lat.get_all_distances(mlwf_center_frac, structure.frac_coords)
distance_min = np.min(distances, axis=1)
distance_mean = np.mean(distance_min)
np.savetxt(f"mlwf_distance_mean_{id}.csv", [distance_mean])

# create the plot data, if wanted
if plot:
    f = open(id + ".win", "a")
    f.write("wannier_plot = true\n")
    f.write("restart = plot\n")
    f.write("wannier_plot_supercell = 3\n")
    f.close()
    os.system(f"mpirun -np {ncores} wannier90.x {id}")

# create the k-point grid that will be interpolated
dummy_calc_data = calc_data(
    structure,
    name,
    id=id,
    ibrav=ibrav,
    calc_type="nscf",
    pw_cutoff=calc_data_wann.pw_cutoff,
    kppa=int(calc_data_wann.kppa),
    pseudo=os.path.join(base_dir, "pseudo", "SG15"),
)

# the simplest way to get a finer interpolation is to change this here
mp_size_interp = " ".join(str(1 * k) for k in dummy_calc_data.k_points_grid)

# generate the k-point string that will be put in the interpolation input
kmesh_wannier = os.popen(f"kmesh.pl {mp_size_interp} wann").read()
kmesh_wann_interp = kmesh_wannier.splitlines()
num_kpt_interp = len(kmesh_wann_interp)
for line_idx in range(0, len(kmesh_wann_interp)):
    kmesh_wann_interp[line_idx] = str(line_idx + 1) + " " + kmesh_wann_interp[line_idx]
kmesh_wann_interp = "\n".join(kmesh_wann_interp)

# write the interpolation input
f = open(id + "_geninterp.kpt", "w")
f.write(
    textwrap.dedent(
        f"""\
k-mesh for interpolation of {id}
crystal
{num_kpt_interp}
{kmesh_wann_interp}
"""
    )
)
f.close()

# run the interpolation
os.system(f"postw90.x {id}")

# parse the output, i.e., read the derivatives, put them in a useful format and save them
interp_df = pd.read_csv(
    f"{id}_geninterp.dat",
    comment="#",
    names=["kpt_idx", "kx", "ky", "kz", "E", "Ex", "Ey", "Ez"],
    delim_whitespace=True,
    skipinitialspace=True,
)
num_kpt = max(interp_df["kpt_idx"])
VB_idx = num_wann - 1

VB_derv = []
for i in range(1, num_kpt):
    curr_block = interp_df[interp_df["kpt_idx"] == i]
    VB_derv.append(
        [
            curr_block["Ex"].values[VB_idx],
            curr_block["Ey"].values[VB_idx],
            curr_block["Ez"].values[VB_idx],
        ]
    )
VB_derv = np.array(VB_derv)
VB_derv = abs(VB_derv)
np.savetxt(
    f"{id}_mean_VB_derivatives.csv",
    [np.mean(VB_derv, axis=0)],
)

# the output folder is huge, delete it if wanted
if delete_out:
    shutil.rmtree(os.path.join(wf_dir, "out"))

# delete the wannier files
delete_files = [".mmn", ".eig", ".chk", ".nnkp", ".amn"]
files = os.listdir()
for file in files:
    for ending in delete_files:
        if ending in file:
            os.remove(file)

# save calculation time
with open("timing.txt", "a+") as f:
    f.write(
        f"{os.path.basename(__file__):<25}  {(time.time() - start_time):7.2f} s  {ncores} cores\n"
    )
