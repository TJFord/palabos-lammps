
#ifndef IBM_LBM_3D_HH
#define IBM_LBM_3D_HH

#include "atom.h"
#include "modify.h"
#include "fix.h"
#include "fix_fcm.h"
#include "update.h"

namespace plb {
  
  template<typename T>
  void weight(T r, std::vector<T> & w){
      T q = sqrt(1 + 4*r*(1-r));
      w[0] = (3 - 2*r - q)/8.0;
      w[1] = (3 - 2*r + q)/8.0;
      w[2] = (1 + 2*r + q)/8.0;
      w[3] = (1 + 2*r - q)/8.0;
  }

  template<typename T, template<typename U> class Descriptor>
  class Interpolation3D: public BoxProcessingFunctional3D_L<T,Descriptor>{
    public:
      Interpolation3D(LammpsWrapper &wrapper_):wrapper(wrapper_){
        plint i,ifix(0),nfix;
        nfix = wrapper.lmp->modify->nfix;
        for (i=0;i<nfix;i++)
          if (strcmp(wrapper.lmp->modify->fix[i]->style,"fcm")==0) ifix=i;
          
        f_fcm = static_cast<LAMMPS_NS::FixFCM *>(wrapper.lmp->modify->fix[ifix]);
        f_fcm->grow_arrays(wrapper.lmp->atom->nmax);
        f_fcm->init();
        groupbit = f_fcm->groupbit;//new code
        dt = wrapper.lmp->update->dt;
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        TensorField3D<T,Descriptor<T>::d> velocity(lattice.getNx(),lattice.getNy(),lattice.getNz());
        //computeVelocity(lattice, velocity);
        //std::cout<<"velocity size "<<lattice.getNx()<<" "<<lattice.getNy()<<" "<<lattice.getNz()<<std::endl;
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt;
        Array<T,3> us(0.,0.,0.);
        Array<T,3> uf;
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        int *mask = wrapper.lmp->atom->mask;//new code
        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        for(ix=domain.x0-2;ix<=domain.x1+2;ix++)
        for(iy=domain.y0-2;iy<=domain.y1+2;iy++)
        for(iz=domain.z0-2;iz<=domain.z1+2;iz++){
          lattice.get(ix,iy,iz).computeVelocity(velocity.get(ix,iy,iz));
        }
        for (plint iS=0; iS<nlocal; iS++){
          if (mask[iS] & groupbit){
            xl = floor(x[iS][0]); 
            yl = floor(x[iS][1]); 
            zl = floor(x[iS][2]);
            rx = x[iS][0] - xl;
            ry = x[iS][1] - yl;
            rz = x[iS][2] - zl;
            weight<T>(rx,wx);
            weight<T>(ry,wy);
            weight<T>(rz,wz);
            us[0] = us[1] = us[2]=0.0;
            for (ii=0;ii<4;ii++ )
              for (jj=0;jj<4;jj++ )
                for (kk=0;kk<4;kk++ ){
                  ix = xl-1 + ii - offset.x ;
                  iy = yl-1 + jj - offset.y ;
                  iz = zl-1 + kk - offset.z ;
          //        if (ix > domain.x1 || ix < domain.x0) continue; 
            //      if (iy > domain.y1 || iy < domain.y0) continue;
              //    if (iz > domain.z1 || iz < domain.z0) continue;
                 // vicinity = 1 can access x0-1 or x1+1 
             /*     if (ix > domain.x1 + 2 || ix < domain.x0 - 2) 
                      std::cout<<" out x "<<domain.x0<<" "<<domain.x1<<" x "<<x[iS][0]<<" nlocal, id "<<nlocal<<" "<<iS<<" ix "<<ix <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                  if (iy > domain.y1 + 2 || iy < domain.y0 - 2) 
                      std::cout<<" out y "<<domain.y0<<" "<<domain.y1<<" y "<<x[iS][1]<<" nlocal, id "<<nlocal<<" "<<iS<<" iy "<<iy <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                  if (iz > domain.z1 + 2 || iz < domain.z0 - 2) 
                      std::cout<<" out z "<<domain.z0<<" "<<domain.z1<<" z "<<x[iS][2]<<" nlocal, id "<<nlocal<<" "<<iS<<" iz "<<iz <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<"zl "<<zl<<" "<<" kk "<<kk<<std::endl;*/
                  uf = velocity.get(ix,iy,iz);
                  wgt = wx[ii]*wy[jj]*wz[kk];
                  us[0] += wgt*uf[0];
                  us[1] += wgt*uf[1];
                  us[2] += wgt*uf[2];
                }
            v[iS][0]=us[0];
            v[iS][1]=us[1];
            v[iS][2]=us[2];
            //Euler method to update position
            x[iS][0] += v[iS][0]*dt;
            x[iS][1] += v[iS][1]*dt;
            x[iS][2] += v[iS][2]*dt;
          }
        }
      }
      virtual Interpolation3D<T,Descriptor> * clone() const{
        return new Interpolation3D(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::nothing; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      class LAMMPS_NS::FixFCM *f_fcm;
      plint groupbit;
      T dt;
  };


  template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity3D(MultiBlockLattice3D<T,Descriptor> &lattice, LammpsWrapper &wrapper)
  {
    //plint envelopeWidth = 2;
    //applyProcessingFunctional(new Interpolation3D<T>(wrapper), velocity.getBoundingBox(),velocity, envelopeWidth); 
    applyProcessingFunctional(new Interpolation3D<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  }


  template<typename T>
  class Interpolation3D_T: public BoxProcessingFunctional3D_T<T,3>{
    public:
      Interpolation3D_T(LammpsWrapper &wrapper_):wrapper(wrapper_){
        plint i,ifix(0),nfix;
        nfix = wrapper.lmp->modify->nfix;
        for (i=0;i<nfix;i++)
          if (strcmp(wrapper.lmp->modify->fix[i]->style,"fcm")==0) ifix=i;
          
        f_fcm = static_cast<LAMMPS_NS::FixFCM *>(wrapper.lmp->modify->fix[ifix]);
        f_fcm->grow_arrays(wrapper.lmp->atom->nmax);
        f_fcm->init();
        groupbit = f_fcm->groupbit;//new code
        dt = wrapper.lmp->update->dt;
      }
      virtual void process(Box3D domain, TensorField3D<T,3> &velocity){
        Dot3D offset = velocity.getLocation();
        //TensorField3D<T,Descriptor<T>::d> velocity(lattice.getNx(),lattice.getNy(),lattice.getNz());
        //computeVelocity(lattice, velocity);
        //std::cout<<"velocity size "<<lattice.getNx()<<" "<<lattice.getNy()<<" "<<lattice.getNz()<<std::endl;
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt;
        Array<T,3> us(0.,0.,0.);
        Array<T,3> uf;
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        int *mask = wrapper.lmp->atom->mask;//new code
        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        /*for(ix=domain.x0-2;ix<=domain.x1+2;ix++)
        for(iy=domain.y0-2;iy<=domain.y1+2;iy++)
        for(iz=domain.z0-2;iz<=domain.z1+2;iz++){
          lattice.get(ix,iy,iz).computeVelocity(velocity.get(ix,iy,iz));
        }*/
        for (plint iS=0; iS<nlocal; iS++){
          if (mask[iS] & groupbit){
            xl = floor(x[iS][0]); 
            yl = floor(x[iS][1]); 
            zl = floor(x[iS][2]);
            rx = x[iS][0] - xl;
            ry = x[iS][1] - yl;
            rz = x[iS][2] - zl;
            weight<T>(rx,wx);
            weight<T>(ry,wy);
            weight<T>(rz,wz);
            us[0] = us[1] = us[2]=0.0;
            for (ii=0;ii<4;ii++ )
              for (jj=0;jj<4;jj++ )
                for (kk=0;kk<4;kk++ ){
                  ix = xl-1 + ii - offset.x ;
                  iy = yl-1 + jj - offset.y ;
                  iz = zl-1 + kk - offset.z ;
          //        if (ix > domain.x1 || ix < domain.x0) continue; 
            //      if (iy > domain.y1 || iy < domain.y0) continue;
              //    if (iz > domain.z1 || iz < domain.z0) continue;
                 // vicinity = 1 can access x0-1 or x1+1 
             /*     if (ix > domain.x1 + 2 || ix < domain.x0 - 2) 
                      std::cout<<" out x "<<domain.x0<<" "<<domain.x1<<" x "<<x[iS][0]<<" nlocal, id "<<nlocal<<" "<<iS<<" ix "<<ix <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                  if (iy > domain.y1 + 2 || iy < domain.y0 - 2) 
                      std::cout<<" out y "<<domain.y0<<" "<<domain.y1<<" y "<<x[iS][1]<<" nlocal, id "<<nlocal<<" "<<iS<<" iy "<<iy <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                  if (iz > domain.z1 + 2 || iz < domain.z0 - 2) 
                      std::cout<<" out z "<<domain.z0<<" "<<domain.z1<<" z "<<x[iS][2]<<" nlocal, id "<<nlocal<<" "<<iS<<" iz "<<iz <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<"zl "<<zl<<" "<<" kk "<<kk<<std::endl;*/
                  uf = velocity.get(ix,iy,iz);
                  wgt = wx[ii]*wy[jj]*wz[kk];
                  us[0] += wgt*uf[0];
                  us[1] += wgt*uf[1];
                  us[2] += wgt*uf[2];
                }
            v[iS][0]=us[0];
            v[iS][1]=us[1];
            v[iS][2]=us[2];
            //Euler method to update position
            x[iS][0] += v[iS][0]*dt;
            x[iS][1] += v[iS][1]*dt;
            x[iS][2] += v[iS][2]*dt;
          }
        }
      }

      /*virtual void process(Box3D domain, TensorField3D<T,3> &velocity){
        Dot3D offset = velocity.getLocation();
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt;
        Array<T,3> us(0.,0.,0.);
        Array<T,3> uf;
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        for (plint iS=0; iS<nlocal; iS++){
          xl = floor(x[iS][0]); 
          yl = floor(x[iS][1]); 
          zl = floor(x[iS][2]);
          rx = x[iS][0] - xl;
          ry = x[iS][1] - yl;
          rz = x[iS][2] - zl;
          weight<T>(rx,wx);
          weight<T>(ry,wy);
          weight<T>(rz,wz);
          us[0] = us[1] = us[2]=0.0;
          for (ii=0;ii<4;ii++ )
            for (jj=0;jj<4;jj++ )
              for (kk=0;kk<4;kk++ ){
                ix = xl-1 + ii - offset.x ;
                iy = yl-1 + jj - offset.y ;
                iz = zl-1 + kk - offset.z ;
        //        if (ix > domain.x1 || ix < domain.x0) continue; 
          //      if (iy > domain.y1 || iy < domain.y0) continue;
            //    if (iz > domain.z1 || iz < domain.z0) continue;
               // vicinity = 1 can access x0-1 or x1+1 
                if (ix > domain.x1 + 2 || ix < domain.x0 - 2) 
                    std::cout<<" out x "<<domain.x0<<" "<<domain.x1<<" x "<<x[iS][0]<<" nlocal, id "<<nlocal<<" "<<iS<<" ix "<<ix <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                if (iy > domain.y1 + 2 || iy < domain.y0 - 2) 
                    std::cout<<" out y "<<domain.y0<<" "<<domain.y1<<" y "<<x[iS][1]<<" nlocal, id "<<nlocal<<" "<<iS<<" iy "<<iy <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                if (iz > domain.z1 + 2 || iz < domain.z0 - 2) 
                    std::cout<<" out z "<<domain.z0<<" "<<domain.z1<<" z "<<x[iS][2]<<" nlocal, id "<<nlocal<<" "<<iS<<" iz "<<iz <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<"zl "<<zl<<" "<<" kk "<<kk<<std::endl;
                uf = velocity.get(ix,iy,iz);
                wgt = wx[ii]*wy[jj]*wz[kk];
                us[0] += wgt*uf[0];
                us[1] += wgt*uf[1];
                us[2] += wgt*uf[2];
              }
          v[iS][0]=us[0];
          v[iS][1]=us[1];
          v[iS][2]=us[2];
        }
      }*/
      virtual Interpolation3D_T<T> * clone() const{
        return new Interpolation3D_T(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::nothing; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      class LAMMPS_NS::FixFCM *f_fcm;
      plint groupbit;
      T dt;
  };


  template<typename T>
  void interpolateVelocity3D(MultiTensorField3D<T,3> &velocity, LammpsWrapper &wrapper)
  {
    //plint envelopeWidth = 2;
    //applyProcessingFunctional(new Interpolation3D<T>(wrapper), velocity.getBoundingBox(),velocity, envelopeWidth); 
    applyProcessingFunctional(new Interpolation3D_T<T>(wrapper), velocity.getBoundingBox(),velocity); 
  }

  template<typename T, template<typename U> class Descriptor>
  class Spreading3D: public BoxProcessingFunctional3D_L<T,Descriptor>{
      public:
      Spreading3D(LammpsWrapper &wrapper_):wrapper(wrapper_){
        plint i,ifix(0),nfix;
        nfix = wrapper.lmp->modify->nfix;
        for (i=0;i<nfix;i++)
          if (strcmp(wrapper.lmp->modify->fix[i]->style,"fcm")==0) ifix=i;
          
        f_fcm = static_cast<LAMMPS_NS::FixFCM *>(wrapper.lmp->modify->fix[ifix]);
        f_fcm->grow_arrays(wrapper.lmp->atom->nmax);
        f_fcm->init();
        groupbit = f_fcm->groupbit;//new code
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt;
        //Array<T,3> ff(0.,0.,0.);
        T **x = wrapper.lmp->atom->x;
        T **f = wrapper.lmp->atom->f;
        int *mask = wrapper.lmp->atom->mask;//new code
        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        for (plint iS=0; iS<nlocal; iS++){
          //f[iS][0]=0.;//debug
          //f[iS][1]=0.;
          //f[iS][2]=0.;
          if(mask[iS] & groupbit){
          xl = floor(x[iS][0]); 
          yl = floor(x[iS][1]); 
          zl = floor(x[iS][2]);
          rx = x[iS][0] - xl;
          ry = x[iS][1] - yl;
          rz = x[iS][2] - zl;
          weight<T>(rx,wx);
          weight<T>(ry,wy);
          weight<T>(rz,wz);
          for (ii=0;ii<4;ii++ )
            for (jj=0;jj<4;jj++ )
              for (kk=0;kk<4;kk++ ){
                ix = xl-1 + ii - offset.x ;
                iy = yl-1 + jj - offset.y ;
                iz = zl-1 + kk - offset.z ;
                //uf = velocity.get(ix,iy,iz);
                /*if (ix > domain.x1 + 2  || ix < domain.x0 - 2) 
                    std::cout<<" spread out x "<<domain.x0<<" "<<domain.x1<<" x "<<x[iS][0]<<" nlocal, id "<<nlocal<<" "<<iS<<" ix "<<ix <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                if (iy > domain.y1 + 2 || iy < domain.y0 - 2) 
                    std::cout<<" spread out y "<<domain.y0<<" "<<domain.y1<<" y "<<x[iS][1]<<" nlocal, id "<<nlocal<<" "<<iS<<" iy "<<iy <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                if (iz > domain.z1 + 2 || iz < domain.z0 - 2) 
                    std::cout<<" spread out z "<<domain.z0<<" "<<domain.z1<<" z "<<x[iS][2]<<" nlocal, id "<<nlocal<<" "<<iS<<" iz "<<iz <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<"zl "<<zl<<" "<<" kk "<<kk<<std::endl;*/
         
        //        if (ix > domain.x1 || ix < domain.x0) continue; 
        //        if (iy > domain.y1 || iy < domain.y0) continue;
        //        if (iz > domain.z1 || iz < domain.z0) continue;
        

                Cell<T,Descriptor>& cell  = lattice.get(ix,iy,iz);
                T *ff=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                wgt = wx[ii]*wy[jj]*wz[kk];
                ff[0] += wgt*f[iS][0]; 
                ff[1] += wgt*f[iS][1]; 
                ff[2] += wgt*f[iS][2]; 
                cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,ff );
              }
          }//mask[iS]
        }
        /*
        pcout<<"lbm x0 "<<domain.x0<<" "<<domain.x1<<" nx "<<lattice.getNx()<<std::endl;
        pcout<<"lbm y0 "<<domain.y0<<" "<<domain.y1<<" ny "<<lattice.getNy()<<std::endl;
        pcout<<"lbm z0 "<<domain.z0<<" "<<domain.z1<<" nz "<<lattice.getNz()<<std::endl;
        for (ix=domain.x0;ix<=domain.x1;ix++)
        for (iy=domain.y0;iy<=domain.y1;iy++)
        for (iz=domain.z0;iz<=domain.z1;iz++){
            Cell<T,Descriptor> &cell = lattice.get(ix+2,iy+2,iz+2); 
            Cell<T,Descriptor> &cell2 = lattice.get(ix-2,iy-2,iz-2); }*/
 
      }
      virtual Spreading3D<T,Descriptor> * clone() const{
        return new Spreading3D(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::staticVariables; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      class LAMMPS_NS::FixFCM *f_fcm;
      plint groupbit;

  };

  template<typename T, template<typename U> class Descriptor>
  void spreadForce3D(MultiBlockLattice3D<T,Descriptor> &lattice,
                   LammpsWrapper &wrapper ){
    //plint envelopeWidth = 2;
    applyProcessingFunctional(new Spreading3D<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  }

  template<typename T, template<typename U> class Descriptor>
  class ForceFSI3D: public BoxProcessingFunctional3D_L<T,Descriptor>{
    public:
      ForceFSI3D(LammpsWrapper &wrapper_):wrapper(wrapper_){
        plint i,ifix(0),nfix;
        nfix = wrapper.lmp->modify->nfix;
        for (i=0;i<nfix;i++)
          if (strcmp(wrapper.lmp->modify->fix[i]->style,"fcm")==0) ifix=i;
          
        f_fcm = static_cast<LAMMPS_NS::FixFCM *>(wrapper.lmp->modify->fix[ifix]);
        f_fcm->grow_arrays(wrapper.lmp->atom->nmax);
        f_fcm->init();
        groupbit = f_fcm->groupbit;//new code
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        TensorField3D<T,Descriptor<T>::d> velocity(lattice.getNx(),lattice.getNy(),lattice.getNz());
        //ScalarField3D<T> density(lattice.getNx(),lattice.getNy(),lattice.getNz());
        //computeVelocity(lattice, velocity);
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt;
        T rho;
        Array<T,3> us(0.,0.,0.);
        Array<T,3> fsi(0.,0.,0.);
        Array<T,3> uf(0.,0.,0.);
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        //T **f = wrapper.lmp->atom->f;
        T **fe = f_fcm->fexternal;
        int *mask = wrapper.lmp->atom->mask;//new code

        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        for(ix=domain.x0-2;ix<=domain.x1+2;ix++)
        for(iy=domain.y0-2;iy<=domain.y1+2;iy++)
        for(iz=domain.z0-2;iz<=domain.z1+2;iz++){
          lattice.get(ix,iy,iz).computeVelocity(velocity.get(ix,iy,iz));
          //density(ix,iy,iz)=lattice.get(ix,iy,iz).computeDensity();
        }
        for (plint iS=0; iS<nlocal; iS++){
          if (mask[iS] & groupbit ){
            xl = floor(x[iS][0]); 
            yl = floor(x[iS][1]); 
            zl = floor(x[iS][2]);
            rx = x[iS][0] - xl;
            ry = x[iS][1] - yl;
            rz = x[iS][2] - zl;
            weight<T>(rx,wx);
            weight<T>(ry,wy);
            weight<T>(rz,wz);
            us[0] = us[1] = us[2]=0.0;
            rho=0.0;
            for (ii=0;ii<4;ii++ )
              for (jj=0;jj<4;jj++ )
                for (kk=0;kk<4;kk++ ){
                  ix = xl-1 + ii - offset.x ;
                  iy = yl-1 + jj - offset.y ;
                  iz = zl-1 + kk - offset.z ;
          //        if (ix > domain.x1 || ix < domain.x0) continue; 
            //      if (iy > domain.y1 || iy < domain.y0) continue;
              //    if (iz > domain.z1 || iz < domain.z0) continue;
                 // vicinity = 1 can access x0-1 or x1+1 
             /*     if (ix > domain.x1 + 2 || ix < domain.x0 - 2) 
                      std::cout<<" out x "<<domain.x0<<" "<<domain.x1<<" x "<<x[iS][0]<<" nlocal, id "<<nlocal<<" "<<iS<<" ix "<<ix <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                  if (iy > domain.y1 + 2 || iy < domain.y0 - 2) 
                      std::cout<<" out y "<<domain.y0<<" "<<domain.y1<<" y "<<x[iS][1]<<" nlocal, id "<<nlocal<<" "<<iS<<" iy "<<iy <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<std::endl;
                  if (iz > domain.z1 + 2 || iz < domain.z0 - 2) 
                      std::cout<<" out z "<<domain.z0<<" "<<domain.z1<<" z "<<x[iS][2]<<" nlocal, id "<<nlocal<<" "<<iS<<" iz "<<iz <<" offset "<<offset.x<<" "<<offset.y<<" "<<offset.z<<" "<<"zl "<<zl<<" "<<" kk "<<kk<<std::endl;*/
                  uf = velocity.get(ix,iy,iz);
                  wgt = wx[ii]*wy[jj]*wz[kk];
                  us[0] += wgt*uf[0];
                  us[1] += wgt*uf[1];
                  us[2] += wgt*uf[2];
                  rho += wgt*lattice.get(ix,iy,iz).computeDensity();
                }
            // assume rho = 1.
            //fsi[0]=us[0]-v[iS][0];       //Eqn.(5) in K. Aidun. Int J. Numer. Meth. Fluids 2010:62:765-783   
            //fsi[1]=us[1]-v[iS][1];          
            //fsi[2]=us[2]-v[iS][2]; 
            // using rho value
            fsi[0]=rho*(us[0]-v[iS][0]);       //Eqn.(5) in K. Aidun. Int J. Numer. Meth. Fluids 2010:62:765-783   
            fsi[1]=rho*(us[1]-v[iS][1]);          
            fsi[2]=rho*(us[2]-v[iS][2]); 
            
            fe[iS][0] = fsi[0]; // no accumulation for external force
            fe[iS][1] = fsi[1];
            fe[iS][2] = fsi[2];
            //std::cout<<"force on particle "<<fsi[0]<<" "<<fsi[1]<<" "<<fsi[2]<<std::endl;

            for (ii=0;ii<4;ii++ )
              for (jj=0;jj<4;jj++ )
                for (kk=0;kk<4;kk++ ){
                  ix = xl-1 + ii - offset.x ;
                  iy = yl-1 + jj - offset.y ;
                  iz = zl-1 + kk - offset.z ;
                 
                  Cell<T,Descriptor>& cell  = lattice.get(ix,iy,iz);
                  T *ff=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                  wgt = wx[ii]*wy[jj]*wz[kk];
                  ff[0] -= wgt*fsi[0]; 
                  ff[1] -= wgt*fsi[1]; 
                  ff[2] -= wgt*fsi[2]; 
                  cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,ff );
                }
          }//mask[is]
        }
      }
      virtual ForceFSI3D<T,Descriptor> * clone() const{
        return new ForceFSI3D(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::staticVariables;
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      class LAMMPS_NS::FixFCM *f_fcm;
      plint groupbit;
  };


  template<typename T, template<typename U> class Descriptor>
  void forceCoupling3D(MultiBlockLattice3D<T,Descriptor> &lattice, LammpsWrapper &wrapper)
  {
    //plint envelopeWidth = 2;
    //applyProcessingFunctional(new Interpolation3D<T>(wrapper), velocity.getBoundingBox(),velocity, envelopeWidth); 
    applyProcessingFunctional(new ForceFSI3D<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  }


}; /* namespace plb */

#endif /* IBDATAEXCHANGEWRAPPERS_HH_LBDEM */
