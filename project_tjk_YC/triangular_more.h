//
// Writen by Xiao Yan Xu (wanderxu@gmail.com)
//
#ifndef __LATTICE_TRIANGULAR_MORE_H_
#define __LATTICE_TRIANGULAR_MORE_H_

#include "latticeplaque.h"

namespace itensor {
struct LatticeBondv2;

using LatticeGraphv2 = std::vector<LatticeBondv2>;

struct LatticeBondv2
    {
    int s1 = 0,
        s2 = 0;
    bool isbd = false;
    std::string type;
    Real x1 = NAN,
         y1 = NAN,
         x2 = NAN,
         y2 = NAN;

    LatticeBondv2() { }

    LatticeBondv2(int s1_, int s2_, bool isbd_)
      : s1{s1_}, 
        s2{s2_},
        isbd{isbd_}
        { }

    LatticeBondv2(int s1_, int s2_,
                Real x1_,  Real y1_,
                Real x2_,  Real y2_)
      : s1{s1_}, 
        s2{s2_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_}
        { }

    LatticeBondv2(int s1_, int s2_, std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        type{type_} 
        { }

    LatticeBondv2(int s1_, int s2_, 
                Real x1_,  Real y1_,
                Real x2_,  Real y2_,
                std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        type{type_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_}
        { }
    };

inline std::ostream& 
operator<<(std::ostream & s, LatticeBondv2 const& b) 
    { 
    //s << format("(%*d,%*d",3,b.s1,3,b.s2);
    s << format("(%d,%d",b.s1,b.s2);
    if(b.type.size()!=0) s << "," << b.type;
    s << ")";
    if(!std::isnan(b.x1) && !std::isnan(b.y1))
        {
        s << format("[%s,%s",b.x1,b.y1);
        if(!std::isnan(b.x2) && !std::isnan(b.y2))
            {
            s << format(";%s,%s]",b.x2,b.y2);
            }
        else
            {
            s << "]";
            }
        }
    return s;
    }

inline std::ostream& 
operator<<(std::ostream& s, LatticeGraphv2 const& G) 
    { 
    for(auto& b : G)
        {
        s << b << "\n";
        }
    return s;
    }

LatticeGraphv2 inline
triangularLatticev2(int Nx, 
                    int Ny,
                  Args const& args = Args::global())
    {
    auto yperiodic = args.getBool("YPeriodic",true);
    // Periodicity on y is meaningless for one dimensional chain or a ladder
    yperiodic = yperiodic && (Ny > 2);
    auto N = Nx*Ny;
    auto Nbond = 3*N-2*Ny + (yperiodic ? 0 : -2*Nx+1);
    LatticeGraphv2 latt; 
    latt.reserve(Nbond);

    for(int n = 1; n <= N; ++n)
        {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        //X-direction bonds
        if(x < Nx) latt.emplace_back(n,n+Ny,false);

        if(Ny > 1) //2d bonds
            {
            // vertical bond 
            if((n+1 <= N) && ((y < Ny)))
                {
                latt.emplace_back(n,n+1,false);
                }
            // Y-periodic diagonal bond
            else if((n+1 <= N) && yperiodic )
                {
                latt.emplace_back(n,n+1,true);
                }

            //Periodic vertical bond
            if(yperiodic && y == 1) latt.emplace_back(n,n+Ny-1,true);

            //Diagonal bonds
            if(x < Nx && y < Ny) latt.emplace_back(n,n+Ny+1,false);
            }
        }

    if(int(latt.size()) != Nbond) Error("Wrong number of bonds");

    return latt;
    }

LatticeGraphv2 inline
triangularLatticeYC(int Nx, 
                    int Ny,
                  Args const& args = Args::global())
  {
    auto yperiodic = args.getBool("YPeriodic",true);
    // Periodicity on y is meaningless for one dimensional chain or a ladder
    yperiodic = yperiodic && (Ny > 2);
    auto N = Nx*Ny;
    auto Nbond = 3*N-2*Ny + (yperiodic ? 0 : -2*Nx+1);
    LatticeGraphv2 latt; 
    latt.reserve(Nbond);

    for(int n = 1; n <= N; ++n) {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        //X-direction bonds
        if(x < Nx) latt.emplace_back(n,n+Ny,false);

        if(Ny > 1){ //2d bonds
            // vertical bond 
            if((n+1 <= N) && ((y < Ny))) {
                latt.emplace_back(n,n+1,false);
            }
            // Y-periodic diagonal bond
            if((n+2*Ny-1 <= N) && y==1 && x%2==1 && yperiodic ) {
                latt.emplace_back(n,n+2*Ny-1,true);
            }
            if((n+1 <= N) && y==Ny && x%2==0 && yperiodic ) {
                latt.emplace_back(n,n+1,true);
            }

            //Periodic vertical bond
            if(yperiodic && y == 1) latt.emplace_back(n,n+Ny-1,true);

            //Diagonal bonds
            if(x < Nx && y > 1 && x%2==1) {
                latt.emplace_back(n,n+Ny-1,false);
            } else if (x < Nx && y < Ny && x%2==0) {
                latt.emplace_back(n,n+Ny+1,false);
            }
        }
    }

    if(int(latt.size()) != Nbond) Error("Wrong number of bonds");

    return latt;
  }

LatticeGraphv2 inline
triangularLatticeXC(int Nx, 
                    int Ny,
                  Args const& args = Args::global())
  {
    auto yperiodic = args.getBool("YPeriodic",true);
    // Periodicity on y is meaningless for one dimensional chain or a ladder
    yperiodic = yperiodic && (Ny > 2);
    auto N = Nx*Ny;
    auto Nbond = 3*N-2*Ny + (yperiodic ? 0 : -2*Nx+1);
    LatticeGraphv2 latt; 
    latt.reserve(Nbond);

    for(int n = 1; n <= N; ++n) {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        //X-direction bonds
        if(x < Nx) latt.emplace_back(n,n+Ny,false);

        if(Ny > 1) { //2d bonds
            // vertical bond 
            if((n+1 <= N) && y%2==1 && (y < Ny)) {
                latt.emplace_back(n,n+1,false);
            }
            if((n-Ny+1 <= N) && y%2==0 && (y < Ny) && x>1) {
                latt.emplace_back(n,n-Ny+1,false);
            }

            // Y-periodic diagonal bond
            if((n+1 <= N) && y%2==1 && y==Ny && yperiodic ) {
                latt.emplace_back(n,n+1,true);
            }
            if((n-Ny+1 <= N) && y%2==0 && y==Ny && yperiodic ) {
                latt.emplace_back(n,n-Ny+1,true);
            }

            //Periodic vertical bond
            if(yperiodic && y==Ny && x>1 && y%2==0) latt.emplace_back(n,n-2*Ny+1,true);
            if(yperiodic && y==Ny && y%2==1) latt.emplace_back(n,n-Ny+1,true);

            //Diagonal bonds
            if( n+1<=N && y%2==0 && y<Ny) latt.emplace_back(n,n+1,false);
            if( n+Ny+1<=N && y%2==1 && y<Ny) latt.emplace_back(n,n+Ny+1,false);
        }
    }

    println("latt.size() = ", latt.size());
    println("Nbond = ", Nbond);
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
    auto Nplaque = (Ny<2 ? 0 : 3*N-4*Ny-Nx+1 + (yperiodic ? 0 : -3*Nx+4));
    Lattice4PlaqueGraph latt; 
    latt.reserve(Nplaque);

    for(int n = 1; n <= N; ++n)
        {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        if(Ny > 1)
            {
            //X-direction plaques, P1
            if(x < Nx-1 && y < Ny) latt.emplace_back(n,n+Ny,n+2*Ny+1,n+Ny+1);
            if(x < Nx-1 && y == Ny && yperiodic) latt.emplace_back(n,n+Ny,n+Ny+1,n+1);

            //Diagonal plaques, P2
            if(x < Nx && y < Ny-1) latt.emplace_back(n,n+Ny+1,n+Ny+2,n+1);
            if(x < Nx && y == Ny-1 && yperiodic) latt.emplace_back(n,n+Ny+1,n+2,n+1);

            //Vertical plaques, P3
            if(x < Nx && y > 1) latt.emplace_back(n,n-1,n+Ny-1,n+Ny);
            if(x < Nx && y == 1 && yperiodic) latt.emplace_back(n,n+Ny-1,n+2*Ny-1,n+Ny);
            }
        }

    if(int(latt.size()) != Nplaque) Error("Wrong number of plaques");

    return latt;
    }

// lattice 4-plaque graph
Lattice4PlaqueGraph inline
triangularLatticeYC4Plaque(int Nx, 
                  int Ny,
                  Args const& args = Args::global())
  {
    auto yperiodic = args.getBool("YPeriodic",true);
    // Periodicity on y is meaningless for one dimensional chain or a ladder
    yperiodic = yperiodic && (Ny > 2);
    auto N = Nx*Ny;
    auto Nplaque = (Ny<2 ? 0 : 3*N-4*Ny + (yperiodic ? 0 : -4*Nx+5));
    Lattice4PlaqueGraph latt; 
    latt.reserve(Nplaque);

    for(int n = 1; n <= N; ++n) {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        if(Ny > 1) {
            if( x%2 == 1 ) {
                //right-bottom plaques
                if(x < Nx && y > 1) latt.emplace_back(n,n-1,n+Ny-1,n+Ny);
                if(x < Nx && y == 1 && yperiodic) latt.emplace_back(n,n+Ny-1,n+2*Ny-1,n+Ny);

                //right-top plaques
                if(x < Nx && y < Ny && y>1) latt.emplace_back(n,n+Ny-1,n+Ny,n+1);
                if(x < Nx && y == Ny && yperiodic) latt.emplace_back(n,n+Ny-1,n+Ny,n-Ny+1);
                if(x < Nx && y == 1 && yperiodic) latt.emplace_back(n,n+2*Ny-1,n+Ny,n+1);

                //right plaques
                if(x < Nx-1 && y > 1) latt.emplace_back(n,n+Ny-1,n+2*Ny,n+Ny);
                if(x < Nx-1 && y == 1 && yperiodic) latt.emplace_back(n,n+2*Ny-1,n+2*Ny,n+Ny);
            } else {
                //right-bottom plaques
                if(x < Nx && y < Ny && y>1) latt.emplace_back(n,n-1,n+Ny,n+Ny+1);
                if(x < Nx && y == Ny && yperiodic) latt.emplace_back(n,n-1,n+Ny,n+1);
                if(x < Nx && y == 1 && yperiodic) latt.emplace_back(n,n+Ny-1,n+Ny,n+Ny+1);

                //right-top plaques
                if(x < Nx && y < Ny) latt.emplace_back(n,n+Ny,n+Ny+1,n+1);
                if(x < Nx && y == Ny && yperiodic) latt.emplace_back(n,n+Ny,n+1,n-Ny+1);

                //right plaques
                if(x < Nx-1 && y < Ny) latt.emplace_back(n,n+Ny,n+2*Ny,n+Ny+1);
                if(x < Nx-1 && y == Ny && yperiodic) latt.emplace_back(n,n+Ny,n+2*Ny,n+1);
            }
        }
    }

    println("latt.size() = ", latt.size());
    println("Nplaque = ", Nplaque);
    if(int(latt.size()) != Nplaque) Error("Wrong number of plaques");

    return latt;
  }

// lattice 4-plaque graph
Lattice4PlaqueGraph inline
triangularLatticeXC4Plaque(int Nx, 
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

    if(Ny%2==1 && yperiodic ) Error("Ny should be even when yperiodic is True!");

    for(int n = 1; n <= N; ++n) {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        if(Ny > 1) {
            if(y%2==1) {
                //X-direction plaques, P1
                if(x < Nx-1 && y < Ny) latt.emplace_back(n,n+Ny,n+2*Ny+1,n+Ny+1);
                if(x < Nx-1 && y == Ny && yperiodic) latt.emplace_back(n,n+Ny,n+Ny+1,n+1);

                //Diagonal plaques, P2
                if(x < Nx && y < Ny-1) latt.emplace_back(n,n+Ny+1,n+2,n+1);
                if(x < Nx && y == Ny-1 && yperiodic) latt.emplace_back(n,n+Ny+1,n+2-Ny,n+1);

                //Vertical plaques, P3
                if(x < Nx-1 && y > 1) latt.emplace_back(n,n+Ny-1,n+2*Ny-1,n+Ny);
                if(x < Nx-1 && y == 1 && yperiodic) latt.emplace_back(n,n+2*Ny-1,n+3*Ny-1,n+Ny);
            } else {
                //X-direction plaques, P1
                if(x < Nx && y < Ny) latt.emplace_back(n,n+Ny,n+Ny+1,n+1);
                if(x < Nx && y == Ny && yperiodic) latt.emplace_back(n,n+Ny,n+1,n-Ny+1);

                //Diagonal plaques, P2
                if(x > 1 && y < Ny-1) latt.emplace_back(n,n+1,n+2,n-Ny+1);
                if(x > 1 && y == Ny-1 && yperiodic) latt.emplace_back(n,n+1,n+2-Ny,n+1-Ny);

                //Vertical plaques, P3
                if(x < Nx && y > 1) latt.emplace_back(n,n-1,n+Ny-1,n+Ny);
            }
        }
    }

    println("latt.size() = ", latt.size());
    println("Nplaque = ", Nplaque);
    if(int(latt.size()) != Nplaque) Error("Wrong number of plaques");

    return latt;
  }

} //namespace itensor

#endif
