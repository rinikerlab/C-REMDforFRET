# C-REMDforFRET

## Publications

J. Phys. Chem. A 2023, 127, 27, 5620–5628 DOI: [https://doi.org/10.1021/acs.jpca.3c01509](https://doi.org/10.1021/acs.jpca.3c01509)

## Abstract

Gas-phase Förster resonance energy transfer (FRET) combines the advantages of mass spectrometry and fluorescence spectroscopy for the conformational analysis of mass-selected biomolecules. While this implementation of FRET in the gas phase promises detailed insights for fundamental and applied studies, the gas-phase environment also poses great challenges. For FRET, fluorophore pairs are typically covalently attached to strategic binding sites in the backbone of a biomolecule, using short linkers. The linker further increases the mobility of the dye, contributing to rotational averaging of the relative orientation of the transition dipole moments of donor and acceptor. However, little is known about the fluorophore’s degrees of freedom in the gas phase and how it may be influenced by intramolecular interactions. In this study, we test the influence of a fluorophore’s linker length on the measured FRET efficiencies in the gas phase to probe the mobility of the fluorophore. An increased FRET efficiency was observed with increasing linker length, ranging from 5.3 % for a linker consisting of 2 atoms to 27.7 % for a linker length of 13 atoms. To rationalize this trend, we profiled the conformational landscape of each model system with MD simulations. Employing state-of-theart enhanced sampling techniques, we captured intramolecular interactions that promote a population shift towards smaller donoracceptor separation for longer linker lengths and induce a significant increase in their acceptor dipole. The presented methodology is a first step towards the explicit consideration of a fluoruophore’s range of motion in the interpretation of gas-phase FRET experiments.

## Install
```bash
mamba create --name lrep -f spec-file.txt
```
or 
```bash
mamba env create --file environment.yml
```

## Analysis

### Get data
The starting structures and topologies are available online in a curated data archive at ETH Zurich (DOI: [10.3929/ethz-b-000600522](https://doi.org/10.3929/ethz-b-000600522)). The simulation data is provided upon reasonable request by the corresponding authors.

### Perform analysis

The data analysis is performed in the [Analysis.ipynb](Analysis/Analysis.ipynb) notebook.
The comparison of QM vs classical dipoles is performed in the [ClassicalvsQCDipoles.ipynb](Analysis/ClassicalvsQCDipoles.ipynb) notebook.


## Simulation

### Get data

The topologies and starting coordinates are provided as a supplementary material with the publication.

How to run a C-REMD simulation is described in the [run_at_300K.py](Simulation/run_at_300K.py) runfile.

For example:
```bash
python run_at_300K.py -n 0.1 -d 0 -f topologies_and_starting_coordinates/ -e 10000 -r plus100ns_start_water_300K_1to025_id0_r0_q1_rep0 -i 0 -ran 0 -qa 1
```

### Post process

Distances, dipoles, and kappas can be extracted using the [get_dipole.py](Analysis/get_dipole.py) runfile.

```bash
python get_dipole.py -in trajectories/plus100ns_start_water_300K_1to025_id3_r0_q1_rep0_4_output.h5 -out test -l L13
```