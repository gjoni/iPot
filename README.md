# iPot

**iPot** (from **i**nverse **Pot**ts) - a set of programs to calculate inter- and intra- protein contact energies using 
atom-atom and residue-residue statistical contact potentials derived from the Potts model [1].

```
**rrce20** - contact energies between residue centroids
**aace18** - atomic contact energies (18 atom types from Ref. [2])
**aace20** - atomic contact energies (20 atom types, all heavy atoms in a residue belong to one type)
**aace167** - atomic contact energies (167 atom types, all heavy atoms in a residue belong to different types)
```


### Description

Residue-residue and atom-atom contact energies were derived by maximizing the 
probability of observing native sequences in a non-redundant set of 6,319 
high-res protein structures. The optimization task was formulated as an inverse 
statistical mechanics problem applied to the Potts model. Its solution by 
pseudo-likelihood maximization provides consistent estimates of coupling 
constants at atomic and residue levels.


## Installation

### Program download and compilation

```
git clone https://github.com/gjoni/aace-rrce
cd aace-rrce
make
```

### Download databases (optional)

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

## Acknowledgements

 - [kdtree library](https://github.com/jtsiomb/kdtree) by John Tsiombikas

## Links

 - [Vakser Lab](http://vakser.compbio.ku.edu/main/)

## References
[1] I Anishchenko, PJ Kundrotas, IA Vakser. Contact potential for structure prediction 
of proteins and protein complexes from Potts model. (2018) In preparation
[2] C Zhang, G Vasmatzis, JL Cornette, C DeLisi. Determination of atomic 
desolvation energies from the structures of crystallized proteins. (1997) 
J Mol Biol. 267(3): 707-6
