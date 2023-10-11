# C-REMDforFRET

## Publications

[p1] J. Phys. Chem. A 2023, 127, 27, 5620–5628 DOI: [https://doi.org/10.1021/acs.jpca.3c01509](https://doi.org/10.1021/acs.jpca.3c01509)

[p2] Benzenberg LR, Katzberger P, et al. Probing the stability of a β-hairpin scaffold after desolvation. ChemRxiv. Cambridge: Cambridge Open Engage; 2023, DOI: TBD

## Abstracts

[p1] Gas-phase Förster resonance energy transfer (FRET) combines mass spectrometry and fluorescence spectroscopy for the conformational analysis of mass-selected biomolecular ions. In FRET, fluorophore pairs are typically covalently attached to a biomolecule using short linkers, which affect the mobility of the dye and the relative orientation of the transition dipole moments of the donor and acceptor. Intramolecular interactions may further influence the range of motion. Yet, little is known about this factor, despite the importance of intramolecular interactions in the absence of a solvent. In this study, we applied transition metal ion FRET (tmFRET) to probe the mobility of a single chromophore pair (Rhodamine 110 and Cu2+) as a function of linker lengths to assess the relevance of intramolecular interactions. Increasing FRET efficiencies were observed with increasing linker length, ranging from 5% (2 atoms) to 28% (13 atoms). To rationalize this trend, we profiled the conformational landscape of each model system using molecular dynamics (MD) simulations. We captured intramolecular interactions that promote a population shift toward smaller donor–acceptor separation for longer linker lengths and induce a significant increase in the acceptor’s transition dipole moment. The presented methodology is a first step toward the explicit consideration of a fluorophore’s range of motion in the interpretation of gas-phase FRET experiments.


[p2] Soft ionization techniques and native mass spectrometry (nMS) open up the possibility to investigate the folding of biomolecular ions in the gas phase, promising insights into the intrinsic conformation in the absence of solvents. Although non-covalent interactions and thus native fold features are believed to be largely retained upon desolvation, the conformation usually depends heavily on the charge state of the species investigated. Here, we utilize transition metal ion Förster Resonance Energy Transfer (tmFRET) and ion mobility-mass spectrometry (IM-MS) to probe the folding of the β-hairpin structure of GB1p in the gas phase. Fluorophore lifetime measurements and collisional cross sections suggest an unfolding of the β-hairpin motif for higher charge states, as expected for an increased Coulomb repulsion that overcompensates for residual intramolecular non-covalent interactions. Comparison of experimental results with computed values from molecular dynamics (MD) simulations allows for a suggestion of probable gas-phase structure candidates. MD simulations were in good agreement with experimental results and suggested that the β-hairpin motif could be retained in the gas phase. Intriguingly, a preservation of the β-hairpin was not only suggested for the 2+ state but also for the 4+ state.

## Install
```bash
mamba create --name lrep -f spec-file.txt
```
or 
```bash
mamba env create --file environment.yml
```

## Analysis and Simulation examples are provided for both publications [p1] and [p2]. The files used in publication [p2] contain a p2 in the filename.

## Analysis

### Get data
The simulation data is provided upon reasonable request by the corresponding authors.

### Perform analysis

The data analysis is performed in the [Analysis_p2.ipynb](Analysis/Analysis.ipynb) notebook.

## Simulation

### Get data

The topologies and starting coordinates are provided as a supplementary material with the publication.

How to run a C-REMD simulation is described in the [run_at_300K_p2.py](Simulation/run_at_300K_p2.py) runfile.

For example:
```bash
python run_at_300K_p2.py -n 10 -d 0 -f topologies_and_starting_coordinates/ -e 10000 -r beta_id0_r0 -i 0 -ran 0
```

### Post process

Distances, dipoles, and kappas can be extracted using the [get_dipole_p2.py](Analysis/get_dipole_p2.py) runfile.

```bash
python get_dipole_p2.py -in input_trajectory -out output_file
```

CCS values can be extracted using the [get_shape_p2.py](Analysis/get_shape_p2.py) runfile. *NOTE:* Impact program must be installed. 

```bash
python get_shape_p2.py -in input_trajectory -out output_file
```

DSSP values can be extracted using the [get_dssp_p2.py](Analysis/get_dssp_p2.py) runfile.

```bash
python get_dssp_p2.py -in input_trajectory -out output_file
```