# aace-rrce

A set of programs to calculate inter- and intra- protein contact energies using 
atom-atom and residue-residue statistical potentials.

Residue-residue and atom-atom contact energies were derived by maximizing the 
probability of observing native sequences in a non-redundant set of 6,319 
high-res protein structures. The optimization task was formulated as an inverse 
statistical mechanics problem applied to the Potts model. Its solution by 
pseudo-likelihood maximization provides consistent estimates of coupling 
constants at atomic and residue levels.

## Installation

### Compilation

        make
        make install

### Download database and set path


## Usage

### Atom-atom energies:

Intrachain:

        ./aace -r <structure.pdb> -t <AACE_TYPE> -d <dmax> -k <kmin>

Interchain:

        ./aace -r <receptor.pdb> -l <ligand.pdb> -t <AACE_TYPE> -d <dmax> -k <kmin>

### Residue-residue energies:

        ./rrce -r <structure.pdb> -t <RRCE_TYPE> -d <dmax> -k <kmin>
        ./rrce -r <receptor.pdb> -l <ligand.pdb> -t <RRCE_TYPE> -d <dmax> -k <kmin>

## References
[1] I Anishchenko, PJ Kundrotas, IA Vakser. Contact energies in proteins and 
    protein complexes inferred from the Potts model. (2017) In preparation
