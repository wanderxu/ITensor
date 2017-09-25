#include "itensor/all.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");
    //Read in individual parameters from the input file
    //second argument to getXXX methods is a default
    //in case parameter not provided in input file
    auto Nx = input.getInt("Nx");
    auto Ny = input.getInt("Ny");
    auto N = Nx*Ny;
    auto yperiodic = input.getYesNo("yperiodic",true);
    auto J1 = input.getInt("J1");
    auto J2 = input.getInt("J2");
    auto gamma1 = input.getInt("gamma1");
    auto gamma2 = input.getInt("gamma2");
    auto quiet = input.getYesNo("quiet",true);

    // Read the sweeps parameters
    auto nsweeps = input.getInt("nsweeps");
    auto table = InputGroup(input,"sweeps");

    //Create the sweeps class & print
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = SpinHalf(N);

    //
    // Use the AutoMPO feature to create the 
    // nearest-neighbor XXZ model with ring exchange.
    //

    auto ampo = AutoMPO(sites);
    auto lattice = triangularLattice(Nx,Ny,{"YPeriodic=",yperiodic});
    auto lattice4plaque = triangularLattice4Plaque(Nx,Ny,{"YPeriodic=",yperiodic});
    for(auto bnd : lattice)
    {
        println( bnd.s1, " ", bnd.s2 );
    }
    println("Total number of nn bound: ", lattice.size());

    for(auto bnd : lattice4plaque)
    {
        println( bnd.s1, " ", bnd.s2, " ", bnd.s3, " ", bnd.s4 );
    }
    println("Total number of plaques: ", lattice4plaque.size());

    // two-body term, nearest neighbor
    for(auto bnd : lattice)
        {
        ampo += J1/8.0,"S+",bnd.s1,"S-",bnd.s2;
        ampo += J1/8.0,"S-",bnd.s1,"S+",bnd.s2;
        ampo += J1*gamma1/4.0,"Sz",bnd.s1,"Sz",bnd.s2;
        }


    // ring-exchange term
    for(auto bnd : lattice4plaque)
        {
        ampo += J2/32.0,"S+",bnd.s1,"S-",bnd.s2,"S+",bnd.s3,"S-",bnd.s4;
        ampo += J2/32.0,"S-",bnd.s1,"S+",bnd.s2,"S-",bnd.s3,"S+",bnd.s4;
        ampo += J2*gamma2/32.0,"S+",bnd.s1,"S-",bnd.s2,"Sz",bnd.s3,"Sz",bnd.s4;
        ampo += J2*gamma2/32.0,"S-",bnd.s1,"S+",bnd.s2,"Sz",bnd.s3,"Sz",bnd.s4;
        ampo += J2*gamma2/32.0,"Sz",bnd.s1,"Sz",bnd.s2,"S+",bnd.s3,"S-",bnd.s4;
        ampo += J2*gamma2/32.0,"Sz",bnd.s1,"Sz",bnd.s2,"S-",bnd.s3,"S+",bnd.s4;
        ampo += J2*gamma2/32.0,"Sz",bnd.s1,"S+",bnd.s2,"S-",bnd.s3,"Sz",bnd.s4;
        ampo += J2*gamma2/32.0,"Sz",bnd.s1,"S-",bnd.s2,"S+",bnd.s3,"Sz",bnd.s4;
        ampo += J2*gamma2/32.0,"S+",bnd.s1,"Sz",bnd.s2,"Sz",bnd.s3,"S-",bnd.s4;
        ampo += J2*gamma2/32.0,"S-",bnd.s1,"Sz",bnd.s2,"Sz",bnd.s3,"S+",bnd.s4;
        ampo += -J2/32.0,"S+",bnd.s1,"Sz",bnd.s2,"S-",bnd.s3,"Sz",bnd.s4;
        ampo += -J2/32.0,"S-",bnd.s1,"Sz",bnd.s2,"S+",bnd.s3,"Sz",bnd.s4;
        ampo += -J2/32.0,"Sz",bnd.s1,"S+",bnd.s2,"Sz",bnd.s3,"S-",bnd.s4;
        ampo += -J2/32.0,"Sz",bnd.s1,"S-",bnd.s2,"Sz",bnd.s3,"S+",bnd.s4;
        ampo += J2*(2.0*gamma2*gamma2-1.0)/16.0,"Sz",bnd.s1,"Sz",bnd.s2,"Sz",bnd.s3,"Sz",bnd.s4;
        }

    auto H = IQMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    // This choice implicitly sets the global Sz quantum number
    // of the wavefunction to zero. Since it is an IQMPS
    // it will remain in this quantum number sector.
    //
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
        }

    auto psi = IQMPS(state);

    //
    // overlap calculates matrix elements of MPO's with respect to MPS's
    // overlap(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f", overlap(psi,H,psi) );


    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H,sweeps,{"Quiet",quiet});

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing overlap = %.10f", overlap(psi,H,psi) );

    println("\nTotal QN of Ground State = ",totalQN(psi));

    return 0;
    }

