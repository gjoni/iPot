# aace-rrce

A set of programs to calculate inter- and intra- protein contact energies using 
atom-atom and residue-residue statistical potentials.

### Description

Residue-residue and atom-atom contact energies were derived by maximizing the 
probability of observing native sequences in a non-redundant set of 6,319 
high-res protein structures. The optimization task was formulated as an inverse 
statistical mechanics problem applied to the Potts model. Its solution by 
pseudo-likelihood maximization provides consistent estimates of coupling 
constants at atomic and residue levels.


## Installation

### Program download and compilation

        git clone https://github.com/gjoni/aace-rrce
        cd aace-rrce
        make

### Download databases and set path

        export TMDOCKDAT=/path/to/database

## Usage

### Calculate atom-atom energies

intrachain:

        ./aace -r <structure.pdb> -t <AACE_TYPE> -d <dmax> -k <kmin>

interchain:

        ./aace -r <receptor.pdb> -l <ligand.pdb> -t <AACE_TYPE> -d <dmax> -k <kmin>

### Calculate residue-residue energies

        ./rrce -r <structure.pdb> -t <RRCE_TYPE> -d <dmax> -k <kmin>
        ./rrce -r <receptor.pdb> -l <ligand.pdb> -t <RRCE_TYPE> -d <dmax> -k <kmin>

## Acknowledgements

This package uses kdtree routine by John Tsiombikas <nuclear@member.fsf.org> available at https://github.com/jtsiomb/kdtree.

## Links

* [Vakser Lab](http://vakser.compbio.ku.edu/main/)

## References
[1] I Anishchenko, PJ Kundrotas, IA Vakser. Contact energies in proteins and 
    protein complexes inferred from the Potts model. (2017) In preparation
