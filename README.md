# cm_nbody
![Mposition](https://user-images.githubusercontent.com/47643825/223049218-2dfe164d-7362-4f1a-a468-46205e225f5a.gif)
This was an undergraduate project that numerically solved for the positions, velocities, and accelerations of an $n\times n \times n$ grid of Xenon atoms arranged in a 3D lattice.  The atoms are subjected to one of two possible forces, the Lennard-Jones and the Morse potential models.
The solutiion is computed via Verlet integration implemented in pure Fortran, compiled with `gfortran` with the included `Makefile`.

A detailed report of the results are in the file `nbodyverlet.pdf`, feel free to peruse there if you are curious!
