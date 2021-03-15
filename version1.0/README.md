# PremDPPI-1.0
## About
<font size=4> 
  
Before running PremDPPI, you need to create a folder including two input files. 

(see the example of 2021031508311331364218043)
  
</font>

## Two input files in the folder of 2021031508311331364218043
<font size=4> 

1. The 3D structure of a protein (1YAG.pdb2), which can be obtained from the Protein Data Bank (PDB) or created by the user.

2. The file containing mutation information (2021031508311331364218043.input), who's name must be consistent with the input folder name.

- PDBfile: coordinate file of a protein structure.
- Partner1: the target protein chains that will be taken into account during the calculation.
- Partner2: the interactor protein chains that will be taken into account during the calculation.
- MutChain: the protein chain where the mutation occurs.
- Mutation_PDB: the first character is one letter amino acid code for the wild-type residue, the second to penultimate characters indicate residue number in MutChain, and the final character indicates the mutant amino acid.
- Result_Id: a number defined by the user.
- isPI: equal to one or zero if mutant structure is produced or not for each mutation.

  The columns are separated by tabs.

</font>
