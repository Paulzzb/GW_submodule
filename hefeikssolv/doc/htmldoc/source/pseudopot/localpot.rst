Local Pseudo Potential
==================================

For a molecule, the local pseudo potential (PP) is formed by the summation
of the local PP of each atom. And it can usually be applied to wavefunctions
via direct multiplication in real space which can also be interpreted as
the multiplication of a diagonal matrix (local PP) with vectors
(wavefunctions). Due to the trans 


The calculation of the local pseudo potential (PP) is split into two parts.
The first part only depends on the atom type and is independent of molecule
settings whereas the second part relies on the location of each atom.

At the beginning of the local PP calculation, the local PP is the
:math:`V_{loc}(r)` function for each type of atom, and the
:math:`V_{loc}(r)` function given by the PP file is on a irregular grid
along the radius direction in the real space.  Meanwhile, by the end of the
local PP calculation, we want a :math:`V_{loc}(\vec{x})` of all atom on the
uniform grid in the real space as an output.


