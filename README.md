# aace-rrce

A set of programs to calculate inter- and intra- protein contact energies using 
atom-atom and residue-residue statistical potentials [1].

### Description

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

### Download databases and set path

        export TMDOCKDAT=/path/to/database

## Usage

### Calculate atom-atom energies

intrachain:

        ./aace -r <structure.pdb> -t <AACE_TYPE> -d <dmax> -k <kmin>

interchain:

        ./aace -r <receptor.pdb> -l <ligand.pdb> -t <AACE_TYPE> -d <dmax> -k <kmin>

### Calculate residue-residue energies

intrachain:

        ./rrce -r <structure.pdb> -t <RRCE_TYPE> -d <dmax> -k <kmin>

interchain:

        ./rrce -r <receptor.pdb> -l <ligand.pdb> -t <RRCE_TYPE> -d <dmax> -k <kmin>

## Links

* [Vakser Lab](http://vakser.compbio.ku.edu/main/)

## References
[1] I Anishchenko, PJ Kundrotas, IA Vakser. Contact potential for structure prediction 
of proteins and protein complexes from Potts model. (2018) In preparation
