The data set *exp_exciton_binding.csv* contains the average exciton binding energies
for each material that we found in our literature review, as well as the associated 
Materials Project identifiers and names from the Materials Project. All binding energies 
are given in meV.

All exciton binding found in the literature review including the references and
experimental method with comments from our side is gathered in the *exp_data_collection.xlsx* 
file in this folder.

We are aware that there are more binding energies available in the literature. Some were left 
out due to uncertainties in the experimental method/data and data analysis method used. Furthermore
we ignored all materials that contained transition metals, thanides and actinides due to the fact 
that strong correlation effects can heavily affect the electronic structure and exciton physics 
of these materials. These effects are difficult to analyse just using the employed standard density 
functional theory. Otherwise we also did not calculate materials like Ge, AlAs, InN as they are metals 
in DFT, which causes some problems in our workflow.

It is very probable that some materials were missed in our literature review. Therefore the
database should be seen as incomplete. If you find more suitable bindung energies in the 
literature please report them to us. Thank you in advance!