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

        make
        make install

## Usage

### Atom-atom energies:
* for a single structure:
        ./aace -r <structure.pdb> -t <AACE_TYPE> -d <dmax> -k <kmin>
        
aace - atom-atom contact energies

rrce - residue-residue contact energies

[1] I Anishchenko, PJ Kundrotas, IA Vakser. Contact energies in proteins and 
    protein complexes inferred from the Potts model. (2017) In preparation
