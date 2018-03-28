# palabos-lammps
This is the coupling interface between palabos and lammps, two widely used opensource code. The interface was built using the immersed boundary method. This approach can achieve high performance on thousands of computers. 
# How to use it:
make sure you have installed a mpi library.
1) download palabos source code
2) download lammps source code
3) download the coupling code
4) add the files in lammps-ext to lammps/src
5) compile lammps as a lib
6) modify the Makefile in the given example, make sure the coupling/src is in the includePaths, the libraryPaths should link to lammps/src
7) compile the embolism.cpp and run it using MPI.  

# Reference:
[Jifu Tan, Talid Sinno, Scott Diamond, A parallel fluid–solid coupling model using LAMMPS and Palabos based on the immersed boundary method, Journal of Computational Science, Vol 25, p89-100, 2018.](https://www.sciencedirect.com/science/article/pii/S187775031730981X)

