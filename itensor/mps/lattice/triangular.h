//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LATTICE_TRIANGULAR_H_
#define __ITENSOR_LATTICE_TRIANGULAR_H_

#include "itensor/mps/lattice/latticebond.h"
#include "itensor/mps/lattice/lattice4plaque.h"

namespace itensor {

LatticeGraph inline
triangularLattice(int Nx, 
                  int Ny,
                  Args const& args = Args::global())
    {
    auto yperiodic = args.getBool("YPeriodic",true);
    // Periodicity on y is meaningless for one dimensional chain or a ladder
    yperiodic = yperiodic && (Ny > 2);
    auto N = Nx*Ny;
    auto Nbond = 3*N-2*Ny + (yperiodic ? 0 : -2*Nx+1);
    LatticeGraph latt; 
    latt.reserve(Nbond);

    for(int n = 1; n <= N; ++n)
        {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        //X-direction bonds
        if(x < Nx) latt.emplace_back(n,n+Ny);

        if(Ny > 1) //2d bonds
            {
            //vertical bond / Y-periodic diagonal bond
            if((n+1 <= N) && ((y < Ny) || yperiodic)) 
                {
                latt.emplace_back(n,n+1);
                }

            //Periodic vertical bond
            if(yperiodic && y == 1) latt.emplace_back(n,n+Ny-1);

            //Diagonal bonds
            if(x < Nx && y < Ny) latt.emplace_back(n,n+Ny+1);
            }
        }

    if(int(latt.size()) != Nbond) Error("Wrong number of bonds");

    return latt;
    }

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
    auto NPlaque = (Ny<2 ? 0 : 3*N-4*Ny-Nx+1 + (yperiodic ? 0 : -3*Nx+4));
    Lattice4PlaqueGraph latt; 
    latt.reserve(NPlaque);

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

    if(int(latt.size()) != NPlaque) Error("Wrong number of bonds");

    return latt;
    }

} //namespace itensor

#endif
