//
// Writen by Xiao Yan Xu (wanderxu@gmail.com)
//
#ifndef __LATTICE_TRIANGULAR_MORE_H_
#define __LATTICE_TRIANGULAR_MORE_H_

#include "latticeplaque.h"

namespace itensor {

// lattice 4-plaque graph
Lattice4PlaqueGraph inline
triangularLattice4Plaque(int Nx, 
                  int Ny,
                  Args const& args = Args::global())
    {
    auto yperiodic = args.getBool("YPeriodic",true);
    // Periodicity on y is meaningless for one dimensional chain or a ladder
    yperiodic = yperiodic && (Ny > 2);
    auto N = Nx*Ny;
    auto Nplaque = (Ny<2 ? 0 : 3*N-4*Ny-Nx+1 + (yperiodic ? 0 : -3*Nx+4));
    Lattice4PlaqueGraph latt; 
    latt.reserve(Nplaque);

    for(int n = 1; n <= N; ++n)
        {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        if(Ny > 1)
            {
            //X-direction plaques
            if(x < Nx-1 && y < Ny) latt.emplace_back(n,n+Ny,n+2*Ny+1,n+Ny+1);
            if(x < Nx-1 && y == Ny && yperiodic) latt.emplace_back(n,n+Ny,n+Ny+1,n+1);

            //Diagonal plaques
            if(x < Nx && y < Ny-1) latt.emplace_back(n,n+Ny+1,n+Ny+2,n+1);
            if(x < Nx && y == Ny-1 && yperiodic) latt.emplace_back(n,n+Ny+1,n+2,n+1);

            //Vertical plaques
            if(x < Nx && y > 1) latt.emplace_back(n,n-1,n+Ny-1,n+Ny);
            if(x < Nx && y == 1 && yperiodic) latt.emplace_back(n,n+Ny-1,n+2*Ny-1,n+Ny);
            }
        }

    if(int(latt.size()) != Nplaque) Error("Wrong number of plaques");

    return latt;
    }

} //namespace itensor

#endif
