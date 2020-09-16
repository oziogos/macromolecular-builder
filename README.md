# macromolecular-builder

A hybrid Monte Carlo/Molecular Mechanics code for the generation of macromolecular initial configurations

* This version currently utilizes a LAMMPS executable for energy minimization; will swap with internal optimization routine

Core functionalities:
- Bond-by-bond construction of aliphatic macromolecules (linear or branched)
- Building elements: -CH3, =CH2, -OH, =O; can reproduce the majority of aliphatic polymers such as vinyls, ethers, esters, acrylics, hydroxyalkanoates, and dienes
- Full conformational and tacticity control
- Introduction of rigid nano-inclusions (with or without flexible functionalization): nanocomposites modelling
- Applications: polymers, polymer matrix nanocomposites, polymer blends, polymer/substrate interfaces, core-shell morphologies, thermotropic/lyotropic liquid crystals
