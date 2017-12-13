#ifndef IBM_LBM_3D_H
#define IBM_LBM_3D_H

#include "lammpsWrapper.h"

namespace plb {
  
 /* template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity3D(MultiBlockLattice3D<T,Descriptor> &lattice,
                           MultiTensorField3D<T,3> &velocity,
                           LammpsWrapper &wrapper);*/
  template<typename T>
  void interpolateVelocity3D(MultiTensorField3D<T,3> &velocity,
                           LammpsWrapper &wrapper);
  
  template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity3D(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LammpsWrapper &wrapper);

  template<typename T, template<typename U> class Descriptor>
  void spreadForce3D(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LammpsWrapper &wrapper);
  
  template<typename T, template<typename U> class Descriptor>
  void forceCoupling3D(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LammpsWrapper &wrapper);
  
}; /* namespace plb */

#include "ibm3D.hh"

#endif /* IBDATAEXCHANGEWRAPPERS_H_LBDEM */
