# Predicting exciton binding energies from ground state properties

This repository was created to analyze how the exciton binding energies of various bulk
semiconducting and insulting materials are related to various ground state properties,
such as the band gap, dielectric constant, ionizity and others. More information can
be found in the associated paper. You may need to adjust the *start_calc* function in 
*src/utils/basic_utils* to suit the supercomputer that you are using or set *calc_setup = "local"* 
in the *main.py* file to perform single calculations on the current interactive node.

**REQUIREMENTS:**

- Python > v3.9 (miniconda with the following packages installed into base)
    - pymatgen
    - matplotlib
    - numpy
    - ase
- QuantumEspresso > v7.1 ("/bin" directory loaded into the path, compiled for MPI)
- Wannier90 > v3.1 ("/bin" and "/utility" directory loaded into the path, compiled for MPI)

**SETUP/USAGE**

There are two ways to supply structures to the workflows, either from the Materials Project or as a .cif file.
If the structure that you want to analyze is on the Materials Project, simply proceed as follows:
Create a file called *api_key.py* in the main directory that contains the following variable: api_key = *YOUR_API_KEY*, 
where *YOUR_API_KEY* can be obtained from https://next-gen.materialsproject.org/api.
Run the *main.py* script inside the main directory with the material-ids, workflows and settings you want.

Otherwise, if the structures you are interested in are only available as .cif files, simply copy the .cif files into the *cifs* directory
and run the *main_cif.py* script, indicating in it which cif-files you want to run, as well as settings and workflows, 
same as for the *main.py* script.

Possible workflows are found in /src/workflows.
**ACKNOWLEDGEMENT**

We want to thank Miguel A. L. Marques for the provision of the automated symmetry detection aiding the Quantum ESPRESSO workflows.

**LICENSE**

Copyright (c) 2023 Malte Grunert and Max Großmann 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
