# phase-separationAKK


Readme file to submit on Git
Before you download/run any of the code, you should create a folder and call it ‘Phase_Spearation_AKK'
→ Within this directory, you should create two subdirectories called: 
"energy_plots" 
 "snapshots"
Download the code into the 'Phase_Spearation_AKK' folder you created 

The code will be separated into three parts summarized below. 

Generating two lattices 
-tetramers and solvent
-dimers and solvent
-filling lattice with random numbers

This code starts by defining functions to generate random numbers, which will be used to dictate whether a solvent or tetramer lie at a node and whether a solvent or dimer lies in a bond.
Then, the two lattices are initiated: one for the nodes and one for the bonds
After filling the nodes and bonds of the lattices, the overlaid lattices are displayed. These are the images that will be saved into your 'snapshots' folder and can be used to generate movies in ImageJ
####################################################
Assigning energies to lattices 
-interactions and associated energies

The #finding energy code, ranks the types of interactions that are allowed to occur in the system (tetramer-dimer interaction is most favorable and thus assigned the lowest energy)
In this system, we define the energies as the following interactions:
Molecule-solvent
Molecule-molecule
Solvent-solvent 

####################################################
Simulation (Metropolis monte carlo)
-initial energy states
-swapping of nodes 
-swapping bonds
-new energies

For the simulation, the initial states are determined by the neighboring molecules. 
Each lattice will be able to swap nodes based on the condition that the energy must be lower than the previous calculated energy. 
This is to ensure that nodes and/or bonds that are swapped contribute to decreasing the free energy of the system. 

####################################################
Data Analysis 
-energies of different probabilities
-how cluster size is determined
	(#tetramers - #solvent)/ total molecules in box

