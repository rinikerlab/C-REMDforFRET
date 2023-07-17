# C-REMDforFRET

## Publications

J. Phys. Chem. A 2023, 127, 27, 5620–5628 DOI: [https://doi.org/10.1021/acs.jpca.3c01509](https://doi.org/10.1021/acs.jpca.3c01509)

## Abstract

Gas-phase Förster resonance energy transfer (FRET) combines mass spectrometry and fluorescence spectroscopy for the conformational analysis of mass-selected biomolecular ions. In FRET, fluorophore pairs are typically covalently attached to a biomolecule using short linkers, which affect the mobility of the dye and the relative orientation of the transition dipole moments of the donor and acceptor. Intramolecular interactions may further influence the range of motion. Yet, little is known about this factor, despite the importance of intramolecular interactions in the absence of a solvent. In this study, we applied transition metal ion FRET (tmFRET) to probe the mobility of a single chromophore pair (Rhodamine 110 and Cu2+) as a function of linker lengths to assess the relevance of intramolecular interactions. Increasing FRET efficiencies were observed with increasing linker length, ranging from 5% (2 atoms) to 28% (13 atoms). To rationalize this trend, we profiled the conformational landscape of each model system using molecular dynamics (MD) simulations. We captured intramolecular interactions that promote a population shift toward smaller donor–acceptor separation for longer linker lengths and induce a significant increase in the acceptor’s transition dipole moment. The presented methodology is a first step toward the explicit consideration of a fluorophore’s range of motion in the interpretation of gas-phase FRET experiments.

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