#include "itensor/all.h"
#include "triangular_more.h"
#include "measure.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    println("//////////////////////////");
    println("Reading input file ......\n");
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
    double t1 = input.getReal("t1");
    auto idop = input.getInt("idop");
    double J1 = input.getReal("J1");
    double J2 = input.getReal("J2");
    double gamma1 = input.getReal("gamma1");
    double gamma2 = input.getReal("gamma2");
    auto readmps = input.getYesNo("readmps",false);
    auto eneropt = input.getYesNo("eneropt",true);
    auto twist_ybc = input.getYesNo("twist_ybc",false);
    double ytheta = input.getReal("ytheta",0.0);
    Cplx expitheta = std::exp( Cplx(0.0,std::acos(-1))* ytheta );
    auto domeas = input.getYesNo("domeas",false);
    auto meas_spincorr = input.getYesNo("meas_spincorr",false);
    auto meas_dimer = input.getYesNo("meas_dimer",false);
    auto meas_dxcorr = input.getYesNo("meas_dxcorr",false);
    auto meas_dycorr = input.getYesNo("meas_dycorr",false);
    auto meas_dxycorr = input.getYesNo("meas_dxycorr",false);
    auto meas_chiralcorr = input.getYesNo("meas_chiralcorr",false);
    auto meas_pair = input.getYesNo("meas_pair",false);
    auto meas_pairbubble = input.getYesNo("meas_pairbubble",false);
    auto meas_pairodp = input.getYesNo("meas_pairodp",false);
    auto meas_paircorr = input.getYesNo("meas_paircorr",false);
    auto lfixi0 = input.getYesNo("lfixi0",false);
    auto x0 = input.getInt("x0",0);
    auto quiet = input.getYesNo("quiet",true);
    // as the measurement usually very slow, we divide it into Npart, submit Npart jobs, each one deals with its ithpart
    auto ithpart = input.getInt("ithpart",1);
    auto Npart = input.getInt("Npart",1);

    // Read the sweeps parameters
    auto nsweeps = input.getInt("nsweeps");
    auto table = InputGroup(input,"sweeps");

    // suggested output file name from model parameters
    std::string runlogfile="runlog_";
    runlogfile += "Nx="; runlogfile += std::to_string(Nx);
    runlogfile += "_Ny="; runlogfile += std::to_string(Ny);
    if(yperiodic) runlogfile += "_YPB";
    else runlogfile += "_OBC";
    runlogfile += "_J1="; runlogfile += std::to_string(J1);
    runlogfile += "_J2="; runlogfile += std::to_string(J2);
    runlogfile += "_g1="; runlogfile += std::to_string(gamma1);
    runlogfile += "_g2="; runlogfile += std::to_string(gamma2);
    runlogfile += ".out";

    println("output file name suggested: ", runlogfile);
    println("expitheta = ", expitheta);

    //Create the sweeps class & print
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    tJ sites;
    IQMPS psi;
    IQMPO H;

    // check whether "psi_file" exists
    {
        std::ifstream ifile("psi_file");
        if(ifile) {
            readmps=true;
            println("psi_file is found, readmps is set to true"); }
        else{
            readmps=false;
            println("psi_file is NOT found, readmps is set to false"); }
    }

    if(readmps) {
        println("\n//////////////////////////////////////////////////");
        println("Reading basis, wavefunction and H from files ......");
        //tJ sites;
        readFromFile("sites_file", sites);
        //IQMPS psi(sites);
        psi=IQMPS(sites);
        readFromFile("psi_file", psi);
        //IQMPO H(sites);
        H=IQMPO(sites);
        readFromFile("H_file", H);
        //auto psiHpsi = overlap(psi,H,psi);
        ////auto psiHHpsi = overlap(psi,H,H,psi);
        //printfln(" Intial energy information from input file: ");
        //printfln("\n<psi|H|psi> = %.10f", psiHpsi );
        ////printfln("\n<psi|H^2|psi> = %.10f", psiHHpsi );
        ////printfln("\n<psi|H^2|psi> - <psi|H|psi>^2 = %.10f", psiHHpsi-psiHpsi*psiHpsi );
        //println("\nTotal QN of Ground State = ",totalQN(psi));
    }
    else {
        println("\n//////////////////////////////////////////////////");
        println("Building basis, wavefunction and H from scratch ......\n");
        //
        // Initialize the site degrees of freedom.
        //
        //auto sites = tJ(N);
        sites = tJ(N);

        //
        // Use the AutoMPO feature to create the 
        // nearest-neighbor XXZ model with ring exchange.
        //

        auto ampo = AutoMPO(sites);
        //auto lattice = triangularLatticev2(Nx,Ny,{"YPeriodic=",yperiodic});
        //auto lattice4plaque = triangularLattice4Plaque(Nx,Ny,{"YPeriodic=",yperiodic});
        auto lattice = triangularLatticeXC(Nx,Ny,{"YPeriodic=",yperiodic});
        auto lattice4plaque = triangularLatticeXC4Plaque(Nx,Ny,{"YPeriodic=",yperiodic});

        println("H is made up of ");
        println("\nBound:\n", lattice);
        println("Total number of nn bound: ", lattice.size());

        println("\nPlaque:\n", lattice4plaque);
        println("Total number of plaques: ", lattice4plaque.size());

        // hopping term
        for(auto bnd : lattice)
        {
            if( (not bnd.isbd) or (not twist_ybc) ) {
                ampo += -t1,"Cdagup",bnd.s1,"Cup",bnd.s2;
                ampo += -t1,"Cdagup",bnd.s2,"Cup",bnd.s1;
                ampo += -t1,"Cdagdn",bnd.s1,"Cdn",bnd.s2;
                ampo += -t1,"Cdagdn",bnd.s2,"Cdn",bnd.s1;
            }
            else if ( ytheta == 1.0 ) {
                println( " Impose twist boundary here ", bnd.s1, " ", bnd.s2);
                ampo += t1,"Cdagup",bnd.s1,"Cup",bnd.s2; // Cup1^+ CupL e^{i\theta}
                ampo += t1,"Cdagup",bnd.s2,"Cup",bnd.s1; // CupL^+ Cup1 e^{-i\theta}
                ampo += t1,"Cdagdn",bnd.s1,"Cdn",bnd.s2; // Cdn1^+ CdnL e^{i\theta}
                ampo += t1,"Cdagdn",bnd.s2,"Cdn",bnd.s1; // CdnL^+ Cdn1 e^{-i\theta}
            }
            else {
                println( " Impose twist boundary here ", bnd.s1, " ", bnd.s2);
                ampo += -t1*expitheta,"Cdagup",bnd.s1,"Cup",bnd.s2; // Cup1^+ CupL e^{i\theta}
                ampo += -t1/expitheta,"Cdagup",bnd.s2,"Cup",bnd.s1; // CupL^+ Cup1 e^{-i\theta}
                ampo += -t1*expitheta,"Cdagdn",bnd.s1,"Cdn",bnd.s2; // Cdn1^+ CdnL e^{i\theta}
                ampo += -t1/expitheta,"Cdagdn",bnd.s2,"Cdn",bnd.s1; // CdnL^+ Cdn1 e^{-i\theta}
            }
        }

        // two-body term, nearest neighbor
        for(auto bnd : lattice)
        {
            if( (not bnd.isbd) or (not twist_ybc) ) {
                ampo += J1*2.0,"S+",bnd.s1,"S-",bnd.s2;
                ampo += J1*2.0,"S-",bnd.s1,"S+",bnd.s2;
                ampo += J1*gamma1*4.0,"Sz",bnd.s1,"Sz",bnd.s2;
            }
            else if ( ytheta == 1.0 ) {
                println( " Impose twist boundary here ", bnd.s1, " ", bnd.s2);
                ampo += -J1*2.0,"S+",bnd.s1,"S-",bnd.s2; // S1^+ SL^- e^{i\theta}
                ampo += -J1*2.0,"S-",bnd.s1,"S+",bnd.s2; // S1^- SL^+ e^{-i\theta}
                ampo += J1*gamma1*4.0,"Sz",bnd.s1,"Sz",bnd.s2;
            }
            else {
                println( " Impose twist boundary here ", bnd.s1, " ", bnd.s2);
                ampo += J1*2.0*expitheta,"S+",bnd.s1,"S-",bnd.s2; // S1^+ SL^- e^{i\theta}
                ampo += J1*2.0/expitheta,"S-",bnd.s1,"S+",bnd.s2; // S1^- SL^+ e^{-i\theta}
                ampo += J1*gamma1*4.0,"Sz",bnd.s1,"Sz",bnd.s2;
            }
        }


        // ring-exchange term
        for(auto bnd : lattice4plaque)
            {
            ampo += J2*8.0,"S+",bnd.s1,"S-",bnd.s2,"S+",bnd.s3,"S-",bnd.s4;
            ampo += J2*8.0,"S-",bnd.s1,"S+",bnd.s2,"S-",bnd.s3,"S+",bnd.s4;
            ampo += J2*gamma2*8.0,"S+",bnd.s1,"S-",bnd.s2,"Sz",bnd.s3,"Sz",bnd.s4;
            ampo += J2*gamma2*8.0,"S-",bnd.s1,"S+",bnd.s2,"Sz",bnd.s3,"Sz",bnd.s4;
            ampo += J2*gamma2*8.0,"Sz",bnd.s1,"Sz",bnd.s2,"S+",bnd.s3,"S-",bnd.s4;
            ampo += J2*gamma2*8.0,"Sz",bnd.s1,"Sz",bnd.s2,"S-",bnd.s3,"S+",bnd.s4;
            ampo += J2*gamma2*8.0,"Sz",bnd.s1,"S+",bnd.s2,"S-",bnd.s3,"Sz",bnd.s4;
            ampo += J2*gamma2*8.0,"Sz",bnd.s1,"S-",bnd.s2,"S+",bnd.s3,"Sz",bnd.s4;
            ampo += J2*gamma2*8.0,"S+",bnd.s1,"Sz",bnd.s2,"Sz",bnd.s3,"S-",bnd.s4;
            ampo += J2*gamma2*8.0,"S-",bnd.s1,"Sz",bnd.s2,"Sz",bnd.s3,"S+",bnd.s4;
            ampo += -J2*8.0,"S+",bnd.s1,"Sz",bnd.s2,"S-",bnd.s3,"Sz",bnd.s4;
            ampo += -J2*8.0,"S-",bnd.s1,"Sz",bnd.s2,"S+",bnd.s3,"Sz",bnd.s4;
            ampo += -J2*8.0,"Sz",bnd.s1,"S+",bnd.s2,"Sz",bnd.s3,"S-",bnd.s4;
            ampo += -J2*8.0,"Sz",bnd.s1,"S-",bnd.s2,"Sz",bnd.s3,"S+",bnd.s4;
            ampo += J2*(2.0*gamma2*gamma2-1.0)*16.0,"Sz",bnd.s1,"Sz",bnd.s2,"Sz",bnd.s3,"Sz",bnd.s4;
            }

        //auto H = IQMPO(ampo);
        H = IQMPO(ampo);

        // Set the initial wavefunction matrix product state
        // to be a Neel state.
        //
        // This choice implicitly sets the global Sz quantum number
        // of the wavefunction to zero when idop is even, and to one when idop is odd.
        // Since it is an IQMPS, it will remain in this quantum number sector.
        //
        auto state = InitState(sites);
        int p = N-idop;
        //for(int i = N; i >= 1; --i) {
        //  if(p > i) {
        //    println("Doubly occupying site ",i);
        //    state.set(i,"UpDn");
        //    p -= 2;
        //  } else if(p > 0) {
        //    println("Singly occupying site ",i);
        //    state.set(i,(i%2==1 ? "Up" : "Dn"));
        //    p -= 1;
        //  } else {
        //    state.set(i,"Emp");
        //  }
        //}
        int iper = N/idop;
        for(int i = N; i >= 1; --i) {
          if(i%iper != 0) {
            println("Singly occupying site ",i);
            state.set(i,(p%2==1 ? "Up" : "Dn"));
            p -= 1;
          } else {
            state.set(i,"Emp");
          }
        }

        //auto psi = IQMPS(state);
        psi = IQMPS(state);

        //
        // overlap calculates matrix elements of MPO's with respect to MPS's
        // overlap(psi,H,psi) = <psi|H|psi>
        //
        printfln("\nInitial energy = %.5f\n", overlap(psi,H,psi));
    }


    if(eneropt){
        println("\n//////////////////////////////////////////////////////////////////");
        println("Beigin energy optimization, to get ground state wavefunction ......");
        //
        // Begin the DMRG calculation
        //
        auto energy = dmrg(psi,H,sweeps,{"Quiet",quiet,"WriteM",600});

        // after the MPS converged, write basis, psi, and H to disk
        writeToFile("sites_file", sites);
        writeToFile("psi_file", psi);
        writeToFile("H_file", H);

        //
        // Calculate entanglement entropy
        for ( int b = 1; b <= N-1; ++b ) {
            psi.position(b);
            auto wf = psi.A(b)*psi.A(b+1);
            auto U = psi.A(b);
            IQTensor S, V;
            auto spectrum = svd(wf,U,S,V);
            Real SvN = 0.;
            for ( auto p : spectrum.eigs() ) {
                if ( p > 1E-12 ) SvN += -p*log(p);
            }
            printfln("Across bond b=%d, SvN = %.10f", b, SvN);
        }

        //
        // Print the final energy reported by DMRG
        //
        println("\nTotal QN of Ground State = ",totalQN(psi));
        printfln("\nGround State Energy = %.10f", energy);
        printfln("\n<psi|H|psi> / N = %.10f", energy/N );

        // since calculate overlap(psi,H,H,psi) is extremely memory expensive, we only calculate it when maxm <= max(640, 2560*16/N)
        if( sweeps.maxm( sweeps.nsweep() ) <= std::max(160, 640*16/N) ) { 
            auto psiHHpsi = overlap(psi,H,H,psi);
            printfln("\n<psi|H^2|psi> = %.10f", psiHHpsi );
            printfln("\n<psi|H^2|psi> - <psi|H|psi>^2 = %.10f", psiHHpsi-energy*energy);
            printfln("\n<psi|H^2|psi> / N^2 = %.10f", psiHHpsi/(N*N) );
            printfln("\n( <psi|H^2|psi> - <psi|H|psi>^2 ) / N^2 = %.10f", (psiHHpsi-energy*energy)/(N*N) );
            printfln("\nsqrt( | <psi|H^2|psi> - <psi|H|psi>^2 | ) / N = %.10f", sqrt(abs(psiHHpsi-energy*energy))/N );
        }

    }

    if(domeas && meas_spincorr) {
        println("\n////////////////////////////");
        println("Start to perform measurement of spin correlation\n");
        //
        // Measure Si.Sj of every {i,j}, and total M
        //
        auto totalM = 0.0;
        auto Msquare = 0.0;
        auto Mzsquare = 0.0;
        std::vector<double> SiSj_meas={};
        std::vector<double> SiSjzz_meas={};
        std::vector<double> SiSjpm_meas={};
        std::vector<double> Sz_meas={};
        std::vector<double> Nup_meas={};
        std::vector<double> Ndn_meas={};
        std::vector<double> Ntot_meas={};
        std::vector<Cplx> Sp_meas={};
        std::vector<Cplx> Sm_meas={};
        for ( int i = 1; i <= N; ++i ) {
            //'gauge' the MPS to site i
            psi.position(i); 
            
            //psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation

            // i == j part

            // magnetization
            auto ket = psi.A(i);
            auto bra = dag(prime(ket,Site));
            auto sz_tmp = ((bra*sites.op("Sz",i)*ket).cplx()).real();
            Sz_meas.emplace_back(sz_tmp);
            totalM +=  sz_tmp;

            Nup_meas.emplace_back( ((bra*sites.op("Nup",i)*ket).cplx()).real() );
            Ndn_meas.emplace_back( ((bra*sites.op("Ndn",i)*ket).cplx()).real() );
            Ntot_meas.emplace_back( ((bra*sites.op("Ntot",i)*ket).cplx()).real() );

            auto sp_tmp = (bra*sites.op("S+",i)*ket).cplx();
            Sp_meas.emplace_back(sp_tmp);
            auto sm_tmp = (bra*sites.op("S-",i)*ket).cplx();
            Sm_meas.emplace_back(sm_tmp);

            auto ss_tmp = 0.0;
            ss_tmp += 0.75*(((dag(ket)*ket).cplx()).real());
            SiSj_meas.emplace_back(ss_tmp);
            Msquare += ss_tmp;
            SiSjzz_meas.emplace_back(0.25);
            Mzsquare += 0.25;
            SiSjpm_meas.emplace_back(ss_tmp-0.25);
            println( i, " ", i, " ", ss_tmp );
            
            if ( i < N ) {
                // i != j part
                //index linking i to i+1:
                auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);
   
                auto op_ip = sites.op("S+",i);
                auto op_im = sites.op("S-",i);
                auto op_iz = sites.op("Sz",i);
                auto Cpm = psi.A(i)*op_ip*dag(prime(psi.A(i),Site,ir));
                auto Cmp = psi.A(i)*op_im*dag(prime(psi.A(i),Site,ir));
                auto Czz = psi.A(i)*op_iz*dag(prime(psi.A(i),Site,ir));
                for(int j = i+1; j <= N; ++j) {
                    Cpm *= psi.A(j);
                    Cmp *= psi.A(j);
                    Czz *= psi.A(j);

                    auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);

                    auto op_jm = sites.op("S-",j);
                    auto ss_tmp = 0.0;
                    ss_tmp += 0.5*(( (Cpm*op_jm)*dag(prime(psi.A(j),jl,Site)) ).cplx()).real();

                    auto op_jp = sites.op("S+",j);
                    ss_tmp += 0.5*(( (Cmp*op_jp)*dag(prime(psi.A(j),jl,Site)) ).cplx()).real();

                    auto spm_tmp = ss_tmp;
                    SiSjpm_meas.emplace_back(ss_tmp);

                    auto op_jz = sites.op("Sz",j);
                    ss_tmp += (( (Czz*op_jz)*dag(prime(psi.A(j),jl,Site)) ).cplx()).real();

                    SiSjzz_meas.emplace_back(ss_tmp-spm_tmp);
                    Mzsquare += (ss_tmp - spm_tmp)*2.0;
                    
                    SiSj_meas.emplace_back(ss_tmp);
                    Msquare += ss_tmp*2.0;
                    println( i, " ", j, " ", ss_tmp ); 

                    if(j < N) {
                        Cpm *= dag(prime(psi.A(j),Link));
                        Cmp *= dag(prime(psi.A(j),Link));
                        Czz *= dag(prime(psi.A(j),Link));
                    }
                }
            }
        }
        printfln("Total M = %.10e", totalM );
        printfln("Msquare = %.10e", Msquare );
        printfln("Mzsquare = %.10e", Mzsquare );

        std::ofstream fSzout("Siz.out",std::ios::out);
        fSzout.precision(12);
        for (std::vector<double>::const_iterator i = Sz_meas.begin(); i != Sz_meas.end(); ++i)
                fSzout << *i << ' ';

        std::ofstream fNupout("Nup.out",std::ios::out);
        fNupout.precision(12);
        for (std::vector<double>::const_iterator i = Nup_meas.begin(); i != Nup_meas.end(); ++i)
                fNupout << *i << ' ';

        std::ofstream fNdnout("Ndn.out",std::ios::out);
        fNdnout.precision(12);
        for (std::vector<double>::const_iterator i = Ndn_meas.begin(); i != Ndn_meas.end(); ++i)
                fNdnout << *i << ' ';

        std::ofstream fNtotout("Ntot.out",std::ios::out);
        fNtotout.precision(12);
        for (std::vector<double>::const_iterator i = Ntot_meas.begin(); i != Ntot_meas.end(); ++i)
                fNtotout << *i << ' ';

        std::ofstream fSpout("Sip.out",std::ios::out);
        fSpout.precision(12);
        for (std::vector<Cplx>::const_iterator i = Sp_meas.begin(); i != Sp_meas.end(); ++i)
                fSpout << *i << ' ';

        std::ofstream fSmout("Sim.out",std::ios::out);
        fSmout.precision(12);
        for (std::vector<Cplx>::const_iterator i = Sm_meas.begin(); i != Sm_meas.end(); ++i)
                fSmout << *i << ' ';

        std::ofstream fSiSjout("SiSj.out",std::ios::out);
        fSiSjout.precision(12);
        for (std::vector<double>::const_iterator i = SiSj_meas.begin(); i != SiSj_meas.end(); ++i)
                fSiSjout << *i << ' ';

        std::ofstream fSiSjzzout("SiSjzz.out",std::ios::out);
        fSiSjzzout.precision(12);
        for (std::vector<double>::const_iterator i = SiSjzz_meas.begin(); i != SiSjzz_meas.end(); ++i)
                fSiSjzzout << *i << ' ';

        std::ofstream fSiSjpmout("SiSjpm.out",std::ios::out);
        fSiSjpmout.precision(12);
        for (std::vector<double>::const_iterator i = SiSjpm_meas.begin(); i != SiSjpm_meas.end(); ++i)
                fSiSjpmout << *i << ' ';
    }

    if(domeas && meas_dimer) {
        println("\n////////////////////////////");
        println("Start to perform measurement of dimer correlation");
        // 
        // measure dimer correlation
        //
        
        // make the dimer table
        // x-direction
        auto num_x_dimer = (Nx-1)*Ny;
        auto num_y_dimer = N-Ny/2-(yperiodic ? 0 : Nx-1);
        auto num_xy_dimer = N-(Ny+1)/2 - (yperiodic ? 0 : Nx);
        LatticeGraph x_dimer;
        LatticeGraph y_dimer;
        LatticeGraph xy_dimer;
        x_dimer.reserve(num_x_dimer);
        y_dimer.reserve(num_y_dimer);
        xy_dimer.reserve(num_xy_dimer);
        for(int n = 1; n <= N; ++n) {
            int x = (n-1)/Ny+1; 
            int y = (n-1)%Ny+1;

            //X-direction bonds
            if(x < Nx) x_dimer.emplace_back(n,n+Ny);

            if(Ny > 1) { //2d bonds
                // vertical bond 
                if((n+1 <= N) && y%2==1 && (y < Ny)) {
                    y_dimer.emplace_back(n,n+1);
                }
                if((n-Ny+1 <= N) && y%2==0 && (y < Ny) && x>1) {
                    y_dimer.emplace_back(n,n-Ny+1);
                }

                // Y-periodic diagonal bond
                if((n+1 <= N) && y%2==1 && y==Ny && yperiodic ) {
                    xy_dimer.emplace_back(n,n+1);
                }
                if((n-Ny+1 <= N) && y%2==0 && y==Ny && yperiodic ) {
                    xy_dimer.emplace_back(n,n-Ny+1);
                }

                //Periodic vertical bond
                if(yperiodic && y==Ny && x>1 && y%2==0) y_dimer.emplace_back(n,n-2*Ny+1);
                if(yperiodic && y==Ny && y%2==1) y_dimer.emplace_back(n,n-Ny+1);

                //Diagonal bonds
                if( n+1<=N && y%2==0 && y<Ny) xy_dimer.emplace_back(n,n+1);
                if( n+Ny+1<=N && y%2==1 && y<Ny) xy_dimer.emplace_back(n,n+Ny+1);
            }
        }
        println( "x_dimer.size() = ", x_dimer.size());
        println( "num_x_dimer = ", num_x_dimer);
        println( "y_dimer.size() = ", y_dimer.size());
        println( "num_y_dimer = ", num_y_dimer);
        println( "xy_dimer.size() = ", xy_dimer.size());
        println( "num_xy_dimer = ", num_xy_dimer);
        if(int(x_dimer.size()) != num_x_dimer) Error("Wrong number of x_dimer");
        if(int(y_dimer.size()) != num_y_dimer) Error("Wrong number of y_dimer");
        if(int(xy_dimer.size()) != num_xy_dimer) Error("Wrong number of xy_dimer");
        println( "\nx_dimer: \n", x_dimer );
        println( "y_dimer: \n", y_dimer );
        println( "xy_dimer: \n", xy_dimer );

        println("measure x_dimer order parameter");
        std::vector<double> dx_meas={}; // store <Di>
        for(int i = 0; i < int(x_dimer.size()); ++i) {
            auto dimer_meas = 0.0;
            // note conjugation codition is used, assume s1!=s2, and it is really the case here
            dimer_meas += mtwobody(psi,sites,{x_dimer[i].s1, x_dimer[i].s2}, "S+", "S-");
            dimer_meas += mtwobody(psi,sites,{x_dimer[i].s1, x_dimer[i].s2}, "Sz", "Sz");
            dx_meas.emplace_back(dimer_meas);
        }
        // output to file
        std::ofstream fdxout("Dxi.out",std::ios::out);
        fdxout.precision(12);
        for (std::vector<double>::const_iterator i = dx_meas.begin(); i != dx_meas.end(); ++i)
                fdxout << *i << ' ';

    if( meas_dxcorr) {
        // measure x_dimer correlation
        println("measure x_dimer correlation");
        int i_start = 0;
        int i_end = int(x_dimer.size());
        if( Npart > 1 ) {
            i_start = ( int(x_dimer.size())/Npart )*(ithpart-1);
            i_end =  ( int(x_dimer.size())/Npart )*ithpart;
            if( ithpart==Npart ) { i_end =  int(x_dimer.size()); }
        }
        println("tot_numer = ", int(x_dimer.size()), " is divided into ", Npart, " part");
        println("i_start = ", i_start, "  i_end = ", i_end);

        std::vector<double> dxdx_meas( (i_end-i_start)*(x_dimer.size()-i_start+x_dimer.size()-i_end+1)/2 ); // store <DiDj>
        //for(int i = 0; i < int(x_dimer.size()); ++i) {
        for(int i = i_start; i < i_end; ++i) {
            std::vector< std::pair<int,int> > op34pair_vec ={}; // store (opk,opl) pair
            std::vector<int> corr_ind = {};  // index in dxdx_meas
            for(int j = i; j < int(x_dimer.size()); ++j) {
                std::vector<int> sites_tmp = { x_dimer[i].s1, x_dimer[i].s2, x_dimer[j].s1, x_dimer[j].s2 };
                // i,j,k,l, only select i!=j < k!=l
                for (auto n : sites_tmp ) { std::cout << n <<" "; }
                std::cout << '\n';
                if( (sites_tmp[0] != sites_tmp[1]) &&
                    (sites_tmp[0] <  sites_tmp[2]) &&
                    (sites_tmp[0] <  sites_tmp[3]) &&
                    (sites_tmp[1] <  sites_tmp[2]) &&
                    (sites_tmp[1] <  sites_tmp[3]) &&
                    (sites_tmp[2] != sites_tmp[3]) ) {
                    op34pair_vec.emplace_back( std::make_pair( sites_tmp[2], sites_tmp[3] ) );
                    int ind_meas = ( x_dimer.size()-i_start + x_dimer.size() - i + 1)*(i-i_start)/2 + j-i;
                    dxdx_meas[ind_meas] = 0.0; // initial the obserator
                    corr_ind.emplace_back( ind_meas );
                    std::cout << "ddcorr_meas = \n";
                    std::cout << '\n';
                    //std::cout << " j = " << j <<std::endl;
                    //std::cout << " ind_meas = " << ind_meas <<std::endl;
                }
                else {
                    // use ordinary method
                    // calculate correlation, Si*Sj*Sk*Sl
                    auto ddcorr_meas = 0.0;
                    if( (sites_tmp[0] == sites_tmp[1]) ||
                        (sites_tmp[0] == sites_tmp[2]) ||
                        (sites_tmp[0] == sites_tmp[3]) ||
                        (sites_tmp[1] == sites_tmp[2]) ||
                        (sites_tmp[1] == sites_tmp[3]) ||
                        (sites_tmp[2] == sites_tmp[3]) ) {
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S+","S-","S+","S-");
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S+","S-","S-","S+");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S+","S-","Sz","Sz");
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S-","S+","S+","S-");
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S-","S+","S-","S+");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S-","S+","Sz","Sz");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"Sz","Sz","S+","S-");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"Sz","Sz","S-","S+");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"Sz","Sz","Sz","Sz"); 
                    }
                    // use conjugatation condition when i,j,k,l not equal
                    else{
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S+","S-","S+","S-");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S+","S-","S-","S+");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"S+","S-","Sz","Sz");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"Sz","Sz","S+","S-");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"Sz","Sz","Sz","Sz"); 
                    }
                    int ind_meas = ( x_dimer.size()-i_start + x_dimer.size() - i + 1)*(i-i_start)/2 + j-i;
                    //std::cout << " ind_meas = " << ind_meas <<std::endl;
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas);
                    dxdx_meas[ind_meas]=ddcorr_meas;
                }
            } // end for(int j = i; j < int(x_dimer.size()); ++j) {

            if( op34pair_vec.size() > 0 ) {
                //std::cout << " op34pair_vec = " <<std::endl;
                //for (auto rr : op34pair_vec ) { std::cout << rr.first <<" "<< rr.second << '\n'; }
                mfourbody_str(psi, sites, {x_dimer[i].s1,x_dimer[i].s2}, "S+", "S-", op34pair_vec, "S+", "S-", corr_ind, dxdx_meas, 0.50 );
                mfourbody_str(psi, sites, {x_dimer[i].s1,x_dimer[i].s2}, "S+", "S-", op34pair_vec, "S-", "S+", corr_ind, dxdx_meas, 0.50 );
                mfourbody_str(psi, sites, {x_dimer[i].s1,x_dimer[i].s2}, "S+", "S-", op34pair_vec, "Sz", "Sz", corr_ind, dxdx_meas, 1.00 );
                mfourbody_str(psi, sites, {x_dimer[i].s1,x_dimer[i].s2}, "Sz", "Sz", op34pair_vec, "S+", "S-", corr_ind, dxdx_meas, 1.00 );
                mfourbody_str(psi, sites, {x_dimer[i].s1,x_dimer[i].s2}, "Sz", "Sz", op34pair_vec, "Sz", "Sz", corr_ind, dxdx_meas, 1.00 );
            }
        }
        // output to file
        std::ofstream fdxdxout("DxiDxj.out",std::ios::out);
        fdxdxout.precision(12);
        for (std::vector<double>::const_iterator i = dxdx_meas.begin(); i != dxdx_meas.end(); ++i)
                fdxdxout << *i << ' ';
    } // end if( meas_dxcorr) {

        println("measure y_dimer order parameter");
        std::vector<double> dy_meas={}; // store <Di>
        for(int i = 0; i < int(y_dimer.size()); ++i) {
            auto dimer_meas = 0.0;
            // note conjugation codition is used, assume s1!=s2, and it is really the case here
            dimer_meas += mtwobody(psi,sites,{y_dimer[i].s1, y_dimer[i].s2}, "S+", "S-");
            dimer_meas += mtwobody(psi,sites,{y_dimer[i].s1, y_dimer[i].s2}, "Sz", "Sz");
            dy_meas.emplace_back(dimer_meas);
        }
        // output to file
        std::ofstream fdyout("Dyi.out",std::ios::out);
        fdyout.precision(12);
        for (std::vector<double>::const_iterator i = dy_meas.begin(); i != dy_meas.end(); ++i)
                fdyout << *i << ' ';

    if( meas_dycorr) {
        // measure y_dimer correlation
        println("measure y_dimer correlation");
        int i_start = 0;
        int i_end = int(y_dimer.size());
        if( Npart > 1 ) {
            i_start = ( int(y_dimer.size())/Npart )*(ithpart-1);
            i_end =  ( int(y_dimer.size())/Npart )*ithpart;
            if( ithpart==Npart ) { i_end =  int(y_dimer.size()); }
        }
        println("tot_numer = ", int(y_dimer.size()), " is divided into ", Npart, " part");
        println("i_start = ", i_start, "  i_end = ", i_end);

        std::vector<double> dydy_meas( (i_end-i_start)*(y_dimer.size()-i_start+y_dimer.size()-i_end+1)/2 ); // store <DiDj>
        for(int i = i_start; i < i_end; ++i) {
            std::vector< std::pair<int,int> > op34pair_vec ={}; // store (opk,opl) pair
            std::vector<int> corr_ind = {};  // index in dydy_meas
            for(int j = i; j < int(y_dimer.size()); ++j) {
                std::vector<int> sites_tmp = { y_dimer[i].s1, y_dimer[i].s2, y_dimer[j].s1, y_dimer[j].s2 };
                // i,j,k,l, only select i!=j < k!=l
                for (auto n : sites_tmp ) { std::cout << n <<" "; }
                std::cout << '\n';
                if( (sites_tmp[0] != sites_tmp[1]) &&
                    (sites_tmp[0] <  sites_tmp[2]) &&
                    (sites_tmp[0] <  sites_tmp[3]) &&
                    (sites_tmp[1] <  sites_tmp[2]) &&
                    (sites_tmp[1] <  sites_tmp[3]) &&
                    (sites_tmp[2] != sites_tmp[3]) ) {
                    op34pair_vec.emplace_back( std::make_pair( sites_tmp[2], sites_tmp[3] ) );
                    int ind_meas = ( y_dimer.size()-i_start + y_dimer.size() - i + 1)*(i-i_start)/2 + j-i;
                    dydy_meas[ind_meas] = 0.0; // initial the obserator
                    corr_ind.emplace_back( ind_meas );
                    std::cout << "ddcorr_meas = \n";
                    std::cout << '\n';
                    //std::cout << " j = " << j <<std::endl;
                    //std::cout << " ind_meas = " << ind_meas <<std::endl;
                }
                else {
                    // use ordinary method
                    // calculate correlation, Si*Sj*Sk*Sl
                    auto ddcorr_meas = 0.0;
                    if( (sites_tmp[0] == sites_tmp[1]) ||
                        (sites_tmp[0] == sites_tmp[2]) ||
                        (sites_tmp[0] == sites_tmp[3]) ||
                        (sites_tmp[1] == sites_tmp[2]) ||
                        (sites_tmp[1] == sites_tmp[3]) ||
                        (sites_tmp[2] == sites_tmp[3]) ) {
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S+","S-","S+","S-");
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S+","S-","S-","S+");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S+","S-","Sz","Sz");
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S-","S+","S+","S-");
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S-","S+","S-","S+");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S-","S+","Sz","Sz");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"Sz","Sz","S+","S-");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"Sz","Sz","S-","S+");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"Sz","Sz","Sz","Sz"); 
                    }
                    // use conjugatation condition when i,j,k,l not equal
                    else{
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S+","S-","S+","S-");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S+","S-","S-","S+");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"S+","S-","Sz","Sz");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"Sz","Sz","S+","S-");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"Sz","Sz","Sz","Sz"); 
                    }
                    int ind_meas = ( y_dimer.size()-i_start + y_dimer.size() - i + 1)*(i-i_start)/2 + j-i;
                    //std::cout << " ind_meas = " << ind_meas <<std::endl;
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas);
                    dydy_meas[ind_meas]=ddcorr_meas;
                }
            } // end for(int j = i; j < int(y_dimer.size()); ++j) {

            if( op34pair_vec.size() > 0 ) {
                //std::cout << " op34pair_vec = " <<std::endl;
                //for (auto rr : op34pair_vec ) { std::cout << rr.first <<" "<< rr.second << '\n'; }
                mfourbody_str(psi, sites, {y_dimer[i].s1,y_dimer[i].s2}, "S+", "S-", op34pair_vec, "S+", "S-", corr_ind, dydy_meas, 0.50 );
                mfourbody_str(psi, sites, {y_dimer[i].s1,y_dimer[i].s2}, "S+", "S-", op34pair_vec, "S-", "S+", corr_ind, dydy_meas, 0.50 );
                mfourbody_str(psi, sites, {y_dimer[i].s1,y_dimer[i].s2}, "S+", "S-", op34pair_vec, "Sz", "Sz", corr_ind, dydy_meas, 1.00 );
                mfourbody_str(psi, sites, {y_dimer[i].s1,y_dimer[i].s2}, "Sz", "Sz", op34pair_vec, "S+", "S-", corr_ind, dydy_meas, 1.00 );
                mfourbody_str(psi, sites, {y_dimer[i].s1,y_dimer[i].s2}, "Sz", "Sz", op34pair_vec, "Sz", "Sz", corr_ind, dydy_meas, 1.00 );
            }
        }
        // output to file
        std::ofstream fdydyout("DyiDyj.out",std::ios::out);
        fdydyout.precision(12);
        for (std::vector<double>::const_iterator i = dydy_meas.begin(); i != dydy_meas.end(); ++i)
                fdydyout << *i << ' ';
    } // if( meas_dycorr) {

        println("measure xy_dimer order parameter");
        std::vector<double> dxy_meas={}; // store <Di>
        for(int i = 0; i < int(xy_dimer.size()); ++i) {
            auto dimer_meas = 0.0;
            // note conjugation codition is used, assume s1!=s2, and it is really the case here
            dimer_meas += mtwobody(psi,sites,{xy_dimer[i].s1, xy_dimer[i].s2}, "S+", "S-");
            dimer_meas += mtwobody(psi,sites,{xy_dimer[i].s1, xy_dimer[i].s2}, "Sz", "Sz");
            dxy_meas.emplace_back(dimer_meas);
        }
        // output to file
        std::ofstream fdxyout("Dxyi.out",std::ios::out);
        fdxyout.precision(12);
        for (std::vector<double>::const_iterator i = dxy_meas.begin(); i != dxy_meas.end(); ++i)
                fdxyout << *i << ' ';

    if( meas_dxycorr) {
        // measure xy_dimer correlation
        println("measure xy_dimer correlation");
        int i_start = 0;
        int i_end = int(xy_dimer.size());
        if( Npart > 1 ) {
            i_start = ( int(xy_dimer.size())/Npart )*(ithpart-1);
            i_end =  ( int(xy_dimer.size())/Npart )*ithpart;
            if( ithpart==Npart ) { i_end =  int(xy_dimer.size()); }
        }
        println("tot_numer = ", int(xy_dimer.size()), " is divided into ", Npart, " part");
        println("i_start = ", i_start, "  i_end = ", i_end);

        std::vector<double> dxydxy_meas( (i_end-i_start)*(xy_dimer.size()-i_start+xy_dimer.size()-i_end+1)/2 ); // store <DiDj>
        for(int i = i_start; i < i_end; ++i) {
            std::vector< std::pair<int,int> > op34pair_vec ={}; // store (opk,opl) pair
            std::vector<int> corr_ind = {};  // index in dxydxy_meas
            for(int j = i; j < int(xy_dimer.size()); ++j) {
                std::vector<int> sites_tmp = { xy_dimer[i].s1, xy_dimer[i].s2, xy_dimer[j].s1, xy_dimer[j].s2 };
                // i,j,k,l, only select i!=j < k!=l
                for (auto n : sites_tmp ) { std::cout << n <<" "; }
                std::cout << '\n';
                if( (sites_tmp[0] != sites_tmp[1]) &&
                    (sites_tmp[0] <  sites_tmp[2]) &&
                    (sites_tmp[0] <  sites_tmp[3]) &&
                    (sites_tmp[1] <  sites_tmp[2]) &&
                    (sites_tmp[1] <  sites_tmp[3]) &&
                    (sites_tmp[2] != sites_tmp[3]) ) {
                    op34pair_vec.emplace_back( std::make_pair( sites_tmp[2], sites_tmp[3] ) );
                    int ind_meas = ( xy_dimer.size()-i_start + xy_dimer.size() - i + 1)*(i-i_start)/2 + j-i;
                    dxydxy_meas[ind_meas] = 0.0; // initial the obserator
                    corr_ind.emplace_back( ind_meas );
                    std::cout << "ddcorr_meas = \n";
                    std::cout << '\n';
                    //std::cout << " j = " << j <<std::endl;
                    //std::cout << " ind_meas = " << ind_meas <<std::endl;
                }
                else {
                    // use ordinary method
                    // calculate correlation, Si*Sj*Sk*Sl
                    auto ddcorr_meas = 0.0;
                    if( (sites_tmp[0] == sites_tmp[1]) ||
                        (sites_tmp[0] == sites_tmp[2]) ||
                        (sites_tmp[0] == sites_tmp[3]) ||
                        (sites_tmp[1] == sites_tmp[2]) ||
                        (sites_tmp[1] == sites_tmp[3]) ||
                        (sites_tmp[2] == sites_tmp[3]) ) {
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S+","S-","S+","S-");
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S+","S-","S-","S+");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S+","S-","Sz","Sz");
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S-","S+","S+","S-");
                        ddcorr_meas += 0.25*mfourbody(psi,sites,sites_tmp,"S-","S+","S-","S+");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S-","S+","Sz","Sz");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"Sz","Sz","S+","S-");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"Sz","Sz","S-","S+");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"Sz","Sz","Sz","Sz"); 
                    }
                    // use conjugatation condition when i,j,k,l not equal
                    else{
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S+","S-","S+","S-");
                        ddcorr_meas += 0.50*mfourbody(psi,sites,sites_tmp,"S+","S-","S-","S+");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"S+","S-","Sz","Sz");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"Sz","Sz","S+","S-");
                        ddcorr_meas += 1.00*mfourbody(psi,sites,sites_tmp,"Sz","Sz","Sz","Sz"); 
                    }
                    int ind_meas = ( xy_dimer.size()-i_start + xy_dimer.size() - i + 1)*(i-i_start)/2 + j-i;
                    //std::cout << " ind_meas = " << ind_meas <<std::endl;
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas);
                    dxydxy_meas[ind_meas]=ddcorr_meas;
                }
            } // end for(int j = i; j < int(xy_dimer.size()); ++j) {

            if( op34pair_vec.size() > 0 ) {
                //std::cout << " op34pair_vec = " <<std::endl;
                //for (auto rr : op34pair_vec ) { std::cout << rr.first <<" "<< rr.second << '\n'; }
                mfourbody_str(psi, sites, {xy_dimer[i].s1,xy_dimer[i].s2}, "S+", "S-", op34pair_vec, "S+", "S-", corr_ind, dxydxy_meas, 0.50 );
                mfourbody_str(psi, sites, {xy_dimer[i].s1,xy_dimer[i].s2}, "S+", "S-", op34pair_vec, "S-", "S+", corr_ind, dxydxy_meas, 0.50 );
                mfourbody_str(psi, sites, {xy_dimer[i].s1,xy_dimer[i].s2}, "S+", "S-", op34pair_vec, "Sz", "Sz", corr_ind, dxydxy_meas, 1.00 );
                mfourbody_str(psi, sites, {xy_dimer[i].s1,xy_dimer[i].s2}, "Sz", "Sz", op34pair_vec, "S+", "S-", corr_ind, dxydxy_meas, 1.00 );
                mfourbody_str(psi, sites, {xy_dimer[i].s1,xy_dimer[i].s2}, "Sz", "Sz", op34pair_vec, "Sz", "Sz", corr_ind, dxydxy_meas, 1.00 );
            }
        }
        // output to file
        std::ofstream fdxydxyout("DxyiDxyj.out",std::ios::out);
        fdxydxyout.precision(12);
        for (std::vector<double>::const_iterator i = dxydxy_meas.begin(); i != dxydxy_meas.end(); ++i)
                fdxydxyout << *i << ' ';
    }  // if( meas_dxycorr) {
    }

    if(domeas && meas_chiralcorr) {
        println("\n////////////////////////////");
        println("Start to perform measurement of chiral correlation\n");
        // 
        // measure chiral correlation
        //

        // make the chiral table
        auto num_tri_plaq = 2*(Nx-1)*(yperiodic ? Ny : Ny-1);
        Lattice3PlaqueGraph tri_plaq;
        tri_plaq.reserve(num_tri_plaq);
        for(int n = 1; n <= N; ++n) {
            int x = (n-1)/Ny+1;
            int y = (n-1)%Ny+1;

            if( y%2 == 1 ) {
                if((x < Nx) && (y > 1))  tri_plaq.emplace_back(n, n+Ny-1, n+Ny); // down tri plaq
                if((x < Nx) && (y < Ny)) tri_plaq.emplace_back(n, n+Ny, n+Ny+1); // up tri plaq
                if((x < Nx) && (y == 1) && yperiodic)  tri_plaq.emplace_back(n, n+2*Ny-1, n+Ny);
                if((x < Nx) && (y == Ny) && yperiodic) tri_plaq.emplace_back(n, n+Ny, n+1);
            } else {
                if((x < Nx) && (y < Ny)) tri_plaq.emplace_back(n, n-1, n+Ny); // down tri plaq
                if((x < Nx) && (y < Ny)) tri_plaq.emplace_back(n, n+Ny, n+1); //  up tri plaq
                if((x < Nx) && (y == Ny)) tri_plaq.emplace_back(n, n-1, n+Ny); // down tri plaq
                if((x < Nx) && (y == Ny) && yperiodic) tri_plaq.emplace_back(n, n+Ny, n-Ny+1); // up tri plaq
            }
        }
        println( "tri_plaq.size() = ", tri_plaq.size() );
        println( "num_tri_plaq = ", num_tri_plaq );
        if(int(tri_plaq.size()) != num_tri_plaq) Error("Wrong number of tri_plaq");
        println( "tri_plaq: \n", tri_plaq );

        // measure chiral correlation
        println("measure chiral correlation");
        int i_start = 0;
        int i_end = int(tri_plaq.size());
        if( Npart > 1 ) {
            i_start = ( int(tri_plaq.size())/Npart )*(ithpart-1);
            i_end =  ( int(tri_plaq.size())/Npart )*ithpart;
            if( ithpart==Npart ) { i_end =  int(tri_plaq.size()); }
        }
        println("tot_numer = ", int(tri_plaq.size()), " is divided into ", Npart, " part");
        println("i_start = ", i_start, "  i_end = ", i_end);

        std::vector<double> XiXj_meas( (i_end-i_start)*(tri_plaq.size()-i_start+tri_plaq.size()-i_end+1)/2 ); // store <DiDj>
        std::vector<double> Xi_meas={}; // store <Xi>
        for(int i = i_start; i < i_end; ++i) {
            std::vector< std::vector<int> > op456vec_vec ={}; // store (opl,opm,opn) pair
            std::vector<int> corr_ind = {};  // index in XiXj_meas
            for(int j = i; j < int(tri_plaq.size()); ++j) {
                std::vector<int> sites_tmp = { tri_plaq[i].s1, tri_plaq[i].s2, tri_plaq[i].s3, tri_plaq[j].s1, tri_plaq[j].s2, tri_plaq[j].s3 };
                // i,j,k,l,m,n, only select i!=j!=k < l!=m!=n
                for (auto n : sites_tmp ) { std::cout << n <<" "; }
                std::cout << '\n';
                if( (sites_tmp[0] != sites_tmp[1]) &&
                    (sites_tmp[1] != sites_tmp[2]) &&
                    (sites_tmp[2] != sites_tmp[0]) &&
                    (sites_tmp[0] <  sites_tmp[3]) &&
                    (sites_tmp[0] <  sites_tmp[4]) &&
                    (sites_tmp[0] <  sites_tmp[5]) &&
                    (sites_tmp[1] <  sites_tmp[3]) &&
                    (sites_tmp[1] <  sites_tmp[4]) &&
                    (sites_tmp[1] <  sites_tmp[5]) &&
                    (sites_tmp[2] <  sites_tmp[3]) &&
                    (sites_tmp[2] <  sites_tmp[4]) &&
                    (sites_tmp[2] <  sites_tmp[5]) &&
                    (sites_tmp[3] != sites_tmp[4]) &&
                    (sites_tmp[4] != sites_tmp[5]) &&
                    (sites_tmp[5] != sites_tmp[3]) ) {
                    op456vec_vec.emplace_back( std::initializer_list<int>{sites_tmp[3], sites_tmp[4], sites_tmp[5] } );
                    int ind_meas = ( tri_plaq.size()-i_start + tri_plaq.size() - i + 1)*(i-i_start)/2 + j-i;
                    XiXj_meas[ind_meas] = 0.0; // initial the obserator
                    corr_ind.emplace_back( ind_meas );
                    std::cout << "XXcorr_meas = \n";
                    std::cout << '\n';
                    //std::cout << " j = " << j <<std::endl;
                    //std::cout << " ind_meas = " << ind_meas <<std::endl;
                }
                else {
                    // use ordinary method
                    // calculate correlation, Si*Sj*Sk*Sl*Sm*Sn
                    if(i == j){
                        auto chiral_meas = 0.0;
                        // note conjugation codition is used, assume s1!=s2!=s3, and it is really the case here
                        // using the conjugate codition
                        chiral_meas +=  (Cplx_i*mthreebodyC(psi,sites,{tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3},"S+","S-","Sz")).real();
                        chiral_meas +=  (Cplx_i*mthreebodyC(psi,sites,{tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3},"S-","Sz","S+")).real();
                        chiral_meas +=  (Cplx_i*mthreebodyC(psi,sites,{tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3},"Sz","S+","S-")).real();
                        Xi_meas.emplace_back(chiral_meas);
                    }
                    auto XXcorr_meas = 0.0;
                    if( (sites_tmp[0] == sites_tmp[1]) ||
                        (sites_tmp[0] == sites_tmp[2]) ||
                        (sites_tmp[0] == sites_tmp[3]) ||
                        (sites_tmp[0] == sites_tmp[4]) ||
                        (sites_tmp[0] == sites_tmp[5]) ||
                        (sites_tmp[1] == sites_tmp[2]) ||
                        (sites_tmp[1] == sites_tmp[3]) ||
                        (sites_tmp[1] == sites_tmp[4]) ||
                        (sites_tmp[1] == sites_tmp[5]) ||
                        (sites_tmp[2] == sites_tmp[3]) ||
                        (sites_tmp[2] == sites_tmp[4]) ||
                        (sites_tmp[2] == sites_tmp[5]) ||
                        (sites_tmp[3] == sites_tmp[4]) ||
                        (sites_tmp[3] == sites_tmp[5]) ||
                        (sites_tmp[4] == sites_tmp[5]) ) {
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "S+", "S-", "Sz" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "S-", "S+", "Sz" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "S-", "Sz", "S+" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "S+", "Sz", "S-" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "Sz", "S+", "S-" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "Sz", "S-", "S+" );

XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S-", "S+", "Sz", "S+", "S-", "Sz" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S-", "S+", "Sz", "S-", "S+", "Sz" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S-", "S+", "Sz", "S-", "Sz", "S+" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S-", "S+", "Sz", "S+", "Sz", "S-" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S-", "S+", "Sz", "Sz", "S+", "S-" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S-", "S+", "Sz", "Sz", "S-", "S+" );

XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "S+", "S-", "Sz" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "S-", "S+", "Sz" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "S-", "Sz", "S+" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "S+", "Sz", "S-" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "Sz", "S+", "S-" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "Sz", "S-", "S+" );

XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S+", "Sz", "S-", "S+", "S-", "Sz" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S+", "Sz", "S-", "S-", "S+", "Sz" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S+", "Sz", "S-", "S-", "Sz", "S+" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S+", "Sz", "S-", "S+", "Sz", "S-" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "S+", "Sz", "S-", "Sz", "S+", "S-" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "S+", "Sz", "S-", "Sz", "S-", "S+" );

XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "S+", "S-", "Sz" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "S-", "S+", "Sz" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "S-", "Sz", "S+" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "S+", "Sz", "S-" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "Sz", "S+", "S-" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "Sz", "S-", "S+" );

XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "Sz", "S-", "S+", "S+", "S-", "Sz" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "Sz", "S-", "S+", "S-", "S+", "Sz" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "Sz", "S-", "S+", "S-", "Sz", "S+" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "Sz", "S-", "S+", "S+", "Sz", "S-" );
XXcorr_meas += -0.25*msixbody(psi, sites, sites_tmp, "Sz", "S-", "S+", "Sz", "S+", "S-" );
XXcorr_meas +=  0.25*msixbody(psi, sites, sites_tmp, "Sz", "S-", "S+", "Sz", "S-", "S+" );
                    }
                    // use conjugatation condition when i,j,k,l,m,n not equal
                    else{
XXcorr_meas +=  0.50*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "S+", "S-", "Sz" );
XXcorr_meas += -0.50*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "S-", "S+", "Sz" );
XXcorr_meas +=  0.50*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "S-", "Sz", "S+" );
XXcorr_meas += -0.50*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "S+", "Sz", "S-" );
XXcorr_meas +=  0.50*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "Sz", "S+", "S-" );
XXcorr_meas += -0.50*msixbody(psi, sites, sites_tmp, "S+", "S-", "Sz", "Sz", "S-", "S+" );

XXcorr_meas +=  0.50*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "S+", "S-", "Sz" );
XXcorr_meas += -0.50*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "S-", "S+", "Sz" );
XXcorr_meas +=  0.50*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "S-", "Sz", "S+" );
XXcorr_meas += -0.50*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "S+", "Sz", "S-" );
XXcorr_meas +=  0.50*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "Sz", "S+", "S-" );
XXcorr_meas += -0.50*msixbody(psi, sites, sites_tmp, "S-", "Sz", "S+", "Sz", "S-", "S+" );

XXcorr_meas +=  0.50*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "S+", "S-", "Sz" );
XXcorr_meas += -0.50*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "S-", "S+", "Sz" );
XXcorr_meas +=  0.50*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "S-", "Sz", "S+" );
XXcorr_meas += -0.50*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "S+", "Sz", "S-" );
XXcorr_meas +=  0.50*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "Sz", "S+", "S-" );
XXcorr_meas += -0.50*msixbody(psi, sites, sites_tmp, "Sz", "S+", "S-", "Sz", "S-", "S+" );
                    }
                    int ind_meas = ( tri_plaq.size()-i_start + tri_plaq.size() - i + 1)*(i-i_start)/2 + j-i;
                    //std::cout << " ind_meas = " << ind_meas <<std::endl;
                    XXcorr_meas = -XXcorr_meas;   // Note that "-" comes from i^2
                    printfln("XXcorr_meas = %.8f\n", XXcorr_meas);
                    XiXj_meas[ind_meas]=XXcorr_meas;
                }
            } // end for(int j = i; j < int(tri_plaq.size()); ++j) {

            if( op456vec_vec.size() > 0 ) {
                //std::cout << " op456vec_vec = " <<std::endl;
                //for (auto rr : op456vec_vec ) { std::cout << rr[0] <<" "<< rr[1] <<" "<< rr[2] << '\n'; }
// use conjugate condition, also note the (-1.0) comes from i^2
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S+", "S-", "Sz", op456vec_vec, "S+", "S-", "Sz", corr_ind, XiXj_meas, 0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S+", "S-", "Sz", op456vec_vec, "S-", "S+", "Sz", corr_ind, XiXj_meas,-0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S+", "S-", "Sz", op456vec_vec, "S-", "Sz", "S+", corr_ind, XiXj_meas, 0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S+", "S-", "Sz", op456vec_vec, "S+", "Sz", "S-", corr_ind, XiXj_meas,-0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S+", "S-", "Sz", op456vec_vec, "Sz", "S+", "S-", corr_ind, XiXj_meas, 0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S+", "S-", "Sz", op456vec_vec, "Sz", "S-", "S+", corr_ind, XiXj_meas,-0.50*(-1.0) );

msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S-", "Sz", "S+", op456vec_vec, "S+", "S-", "Sz", corr_ind, XiXj_meas, 0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S-", "Sz", "S+", op456vec_vec, "S-", "S+", "Sz", corr_ind, XiXj_meas,-0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S-", "Sz", "S+", op456vec_vec, "S-", "Sz", "S+", corr_ind, XiXj_meas, 0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S-", "Sz", "S+", op456vec_vec, "S+", "Sz", "S-", corr_ind, XiXj_meas,-0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S-", "Sz", "S+", op456vec_vec, "Sz", "S+", "S-", corr_ind, XiXj_meas, 0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "S-", "Sz", "S+", op456vec_vec, "Sz", "S-", "S+", corr_ind, XiXj_meas,-0.50*(-1.0) );

msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "Sz", "S+", "S-", op456vec_vec, "S+", "S-", "Sz", corr_ind, XiXj_meas, 0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "Sz", "S+", "S-", op456vec_vec, "S-", "S+", "Sz", corr_ind, XiXj_meas,-0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "Sz", "S+", "S-", op456vec_vec, "S-", "Sz", "S+", corr_ind, XiXj_meas, 0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "Sz", "S+", "S-", op456vec_vec, "S+", "Sz", "S-", corr_ind, XiXj_meas,-0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "Sz", "S+", "S-", op456vec_vec, "Sz", "S+", "S-", corr_ind, XiXj_meas, 0.50*(-1.0) );
msixbody_str(psi, sites, {tri_plaq[i].s1,tri_plaq[i].s2,tri_plaq[i].s3}, "Sz", "S+", "S-", op456vec_vec, "Sz", "S-", "S+", corr_ind, XiXj_meas,-0.50*(-1.0) );
            }
        } // end for(int i = 0; i < int(tri_plaq.size()); ++i) {
        std::ofstream fXiXjout("XiXj.out",std::ios::out);
        fXiXjout.precision(12);
        for (std::vector<double>::const_iterator i = XiXj_meas.begin(); i != XiXj_meas.end(); ++i)
                fXiXjout << *i << ' ';

        std::ofstream fXiout("Xi.out",std::ios::out);
        fXiout.precision(12);
        for (std::vector<double>::const_iterator i = Xi_meas.begin(); i != Xi_meas.end(); ++i)
                fXiout << *i << ' ';
    }

    if(domeas && meas_pair) {
        println("\n////////////////////////////");
        println("Start to perform measurement of dimer correlation");
        // 
        // measure pairing correlation
        //

        // set the nnlist
        LatticeNeighborGraph nnlist;
        for(int n = 1; n<=N; ++n) {
            int x = (n-1)/Ny+1;
            int y = (n-1)%Ny+1;
            int id0 = n;
            int id1 = 0, id2 = 0, id3 = 0, id4 = 0, id5 = 0, id6 = 0;
            if(y%2==1) {
                id1 = (x<Nx ? n+Ny: n+Ny-N);    // (1,0) dir
                id2 = (y>1 ? (x<Nx ? n+Ny-1 : N-n+1): (x<Nx ? n+2*Ny-1 : N-n+1) ); // (0,-1) dir
                id3 = (y>1 ? n-1 : n+Ny-1); // (-1,-1) dir
                id4 = (x>1 ? n-Ny : n-Ny+N); // (-1,0) dir
                id5 = (y<Ny ? n+1: n+1-Ny);  // (0,1) dir
                id6 = (x<Nx ? (y<Ny ? n+Ny+1 : n-Ny+1) : ( y<Ny ? n+Ny+1-N : n+Ny+1-N-Ny) ); // (1,1) dir
            } else {
                id1 = (x<Nx ? n+Ny: n+Ny-N); // (1,0) dir
                id2 = n-1;   // (0,-1) dir
                id3 = (x>1 ? n-Ny-1 : N-Ny+n-1);  // (-1,-1) dir
                id4 = (x>1 ? n-Ny: n-Ny+N); // (-1,0) dir
                id5 = (x>1 ? (y<Ny ? n-Ny+1 : n-2*Ny+1) : ( y<Ny ? n-Ny+1+N : n-2*Ny+1+N) ); // (0,1) dir
                id6 = (y<Ny ? n+1 : n-Ny+1);  // (1,1) dir
            }
            std::vector<int> nntmp = {id1, id2, id3, id4, id5, id6};
            nnlist.emplace_back( id0, nntmp );
        }
        println( "\nnnlist: \n", nnlist );

        // measure single-particle Green's function, and get pairbubble from it
        if( meas_pairbubble ) {
            std::vector<Cplx> grup={};
            std::vector<Cplx> grdn={};
            for(int i = 1; i<=N; ++i) {
                psi.position(i);
                auto ket = psi.A(i);
                auto bra = dag(prime(ket,Site));
                grup.emplace_back( ((bra*sites.op("Nup",i)*ket).cplx()) );
                grdn.emplace_back( ((bra*sites.op("Ndn",i)*ket).cplx()) );
                if ( i < N ) {
                    // i < j case
                    // c^+_iup c_jup = a^+_iup F_i ... F_{j-1} a_jup
                    // c^+_idn c_jdn = a^+_idn F_{i+1} ... F_j a_jdn
                    auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);
                    auto gup = noprime( psi.A(i)*sites.op("Adagup",i), Site );
                    auto gdn = psi.A(i)*sites.op("Adagdn",i);
                    gup *= sites.op("F",i);
                    gup *= dag(prime(psi.A(i),Site,ir));
                    gdn *= dag(prime(psi.A(i),Site,ir));
                    for(int j = i+1; j <= N; ++j) {
                        gup *= psi.A(j);
                        gdn *= psi.A(j);
                        gdn *= sites.op("F",j);

                        auto gup_tmp = gup*sites.op("Aup",j);
                        auto jr = commonIndex(psi.A(j),psi.A(j-1),Link);
                        gup_tmp *= dag( prime(psi.A(j),Site,jr) );
                        grup.emplace_back( -gup_tmp.cplx() ); // minus as we apply Adag first

                        auto gdn_tmp = noprime(gdn,Site)*sites.op("Adn",j);
                        gdn_tmp *= dag( prime(psi.A(j),Site,jr) );
                        grdn.emplace_back( -gdn_tmp.cplx() ); // minus as we apply Adag first

                        if( j < N ) {
                            gup *= sites.op("F",j);
                            gup *= dag( prime(psi.A(j), Site, Link) );
                            gdn *= dag( prime(psi.A(j), Site, Link) );
                        }

                    }
                } // if ( i < N ) {
            }  // for(int i = 1; i<=N; ++i) {
            std::ofstream fgrupout("grup.out",std::ios::out);
            fgrupout.precision(12);
            for (std::vector<Cplx>::const_iterator i = grup.begin(); i != grup.end(); ++i)
                    fgrupout << *i << ' ';
            std::ofstream fgrdnout("grdn.out",std::ios::out);
            fgrdnout.precision(12);
            for (std::vector<Cplx>::const_iterator i = grdn.begin(); i != grdn.end(); ++i)
                    fgrdnout << *i << ' ';

            // cal. pair_bubble term
            std::vector<Cplx> pair_bubble={};
            for(int n1 = 1; n1 <= N; ++n1) {
                //int n1 = 1;
                int i = nnlist[n1-1].s0;
                for (int id1 = 0; id1<6; id1++) {
                    //int id1=0;
                    int j = nnlist[n1-1].snn[id1];
                    for(int n2 = n1; n2 <= N ; ++n2) {
                        //int n2 = 8;
                        int k = nnlist[n2-1].s0;
                        for (int id2 = 0; id2<6; id2++) {
                            //int id2=0;
                            int l = nnlist[n2-1].snn[id2];

                            // pair1a4_bubble = < c^+_jdn c_ldn> <c^+_iup c_kup> + <c^+_jup c_lup> <c^+_idn c_kdn>
                            Cplx gjlup;
                            Cplx gjldn;
                            if( j <= l ) {
                                gjlup = grup[(2*N-j+2)*(j-1)/2+l-j];
                                gjldn = grdn[(2*N-j+2)*(j-1)/2+l-j];
                            } else {
                                gjlup = std::conj(grup[(2*N-l+2)*(l-1)/2+j-l]);
                                gjldn = std::conj(grdn[(2*N-l+2)*(l-1)/2+j-l]);
                            }
                            Cplx gikup;
                            Cplx gikdn;
                            if ( i <= k ){
                                gikup = grup[(2*N-i+2)*(i-1)/2+k-i];
                                gikdn = grdn[(2*N-i+2)*(i-1)/2+k-i];
                            } else {
                                gikup = std::conj(grup[(2*N-k+2)*(k-1)/2+i-k]);
                                gikdn = std::conj(grdn[(2*N-k+2)*(k-1)/2+i-k]);
                            }
                            Cplx pair1a4_bubble = gjldn*gikup + gjlup*gikdn;

                            // pair2a3_bubble = < c^+_jup c_kup> <c^+_idn c_ldn> + <c^+_jdn c_kdn> <c^+_iup c_lup>
                            Cplx gilup;
                            Cplx gildn;
                            if( i <= l ) {
                                gilup = grup[(2*N-i+2)*(i-1)/2+l-i];
                                gildn = grdn[(2*N-i+2)*(i-1)/2+l-i];
                            } else {
                                gilup = std::conj(grup[(2*N-l+2)*(l-1)/2+i-l]);
                                gildn = std::conj(grdn[(2*N-l+2)*(l-1)/2+i-l]);
                            }
                            Cplx gjkup;
                            Cplx gjkdn;
                            if ( j <= k ){
                                gjkup = grup[(2*N-j+2)*(j-1)/2+k-j];
                                gjkdn = grdn[(2*N-j+2)*(j-1)/2+k-j];
                            } else {
                                gjkup = std::conj(grup[(2*N-k+2)*(k-1)/2+j-k]);
                                gjkdn = std::conj(grdn[(2*N-k+2)*(k-1)/2+j-k]);
                            }
                            Cplx pair2a3_bubble = gjkup*gildn + gjkdn*gilup;

                            pair_bubble.emplace_back( pair1a4_bubble );
                            pair_bubble.emplace_back( pair2a3_bubble );
                            ////printfln(" %d, %d, %d, %d, pair1a4_bubble = %.12f", j, i, k, l, pair1a4);
                            ////printfln(" %d, %d, %d, %d, pair2a3_bubble = %.12f", j, i, k, l, pair2a3);
                        } // for (int id2 = 0; id2<6; id2++) {
                    }
                    printfln(" pair_bubble of j,i = %d, %d, k,l=?,? done ", j, i);
                } // for (int id1 = 0; id1<6; id1++) {
            } // for(int n1 = 1; n1 <= N ; ++n1) {
            std::ofstream fpair_bubbleout("pair_bubble.out",std::ios::out);
            fpair_bubbleout.precision(12);
            for (std::vector<Cplx>::const_iterator i = pair_bubble.begin(); i != pair_bubble.end(); ++i)
                    fpair_bubbleout << *i << ' ';
        }

        // measure pairing order parameter
        if( meas_pairodp) {
            std::vector<Cplx> pairodp={};
            for(int n = 1; n <= N ; ++n) {
                int i = nnlist[n-1].s0;
                for (int id = 0; id<6; id++) {
                    int j = nnlist[n-1].snn[id];

                    // i < j case:
                    if(i < j ){
                        // println( " i<j case ");
                        psi.position(i); 
                        // c_iup c_jdn = a_iup F_i ... F_j a_jdn
                        auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);
                        auto cpair = noprime( psi.A(i)*sites.op("Aup",i), Site );
                        cpair = cpair * sites.op("F",i) * dag(prime(psi.A(i),Site,ir));
                        for (int k=i+1; k<j; ++k) {
                            cpair *= psi.A(k);
                            cpair = cpair*sites.op("F",k)*dag(prime(psi.A(k),Site,Link));
                        }
                        cpair *= psi.A(j);
                        cpair = noprime(cpair*sites.op("F",j), Site);
                        auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);
                        cpair = cpair*sites.op("Adn",j)*dag(prime(psi.A(j),Site,jl));
                        pairodp.emplace_back(cpair.cplx()); // store it
                        printfln(" %d, %d, pair = %.12f", i,j, cpair.cplx());

                        // c_idn c_jup = - c_jup c_idn = - a_jup F_{j-1} ... F_{i+1} a_idn
                        psi.position(j);
                        ir = commonIndex(psi.A(j),psi.A(j-1),Link);
                        cpair = psi.A(j)*sites.op("Aup",j) * dag(prime(psi.A(j),Site,ir));
                        for (int k=j-1; k>i; --k) {
                            cpair *= psi.A(k);
                            cpair = cpair*sites.op("F",k)*dag(prime(psi.A(k),Site,Link));
                        }
                        cpair *= psi.A(i);
                        jl = commonIndex(psi.A(i+1),psi.A(i),Link);
                        cpair = cpair*sites.op("Adn",i)*dag(prime(psi.A(i),Site,jl));
                        pairodp.emplace_back(-cpair.cplx()); // store it
                        printfln(" %d, %d, pair = %.12f", i,j, -cpair.cplx());
                    // i > j case:
                    } else if ( i > j ) {
                        //println( " i>j case ");
                        // c_iup c_jdn = -c_jdn c_iup = a_iup F_{i-1} ... F_{j+1} a_jdn
                        psi.position(i);
                        auto ir = commonIndex(psi.A(i),psi.A(i-1),Link);
                        auto cpair = psi.A(i)*sites.op("Aup",i) * dag(prime(psi.A(i),Site,ir));
                        for (int k=i-1; k>j; --k) {
                            cpair *= psi.A(k);
                            cpair = cpair*sites.op("F",k)*dag(prime(psi.A(k),Site,Link));
                        }
                        cpair *= psi.A(j);
                        auto jl = commonIndex(psi.A(j+1),psi.A(j),Link);
                        cpair = cpair*sites.op("Adn",j)*dag(prime(psi.A(j),Site,jl));
                        pairodp.emplace_back(cpair.cplx()); // store it
                        printfln(" %d, %d, pair = %.12f", i,j, cpair.cplx());

                        // c_idn c_jup > = - c_jup c_idn = - a_jup F_j ... F_i a_idn
                        ir = commonIndex(psi.A(j),psi.A(j+1),Link);
                        cpair = noprime( psi.A(j)*sites.op("Aup",j), Site );
                        cpair = cpair * sites.op("F",j) * dag(prime(psi.A(j),Site,ir));
                        for (int k=j+1; k<i; ++k) {
                            cpair *= psi.A(k);
                            cpair = cpair*sites.op("F",k)*dag(prime(psi.A(k),Site,Link));
                        }
                        cpair *= psi.A(i);
                        cpair = noprime(cpair*sites.op("F",i), Site);
                        jl = commonIndex(psi.A(i),psi.A(i-1),Link);
                        cpair = cpair*sites.op("Adn",i)*dag(prime(psi.A(i),Site,jl));
                        pairodp.emplace_back(-cpair.cplx()); // store it
                        printfln(" %d, %d, pair = %.12f", i,j, -cpair.cplx());

                    } else {
                        Error("Error: i and j should be different!");
                    }
                }
            }

            // output
            std::ofstream pairodpout("pairodp.out",std::ios::out);
            pairodpout.precision(12);
            for (std::vector<Cplx>::const_iterator i = pairodp.begin(); i != pairodp.end(); ++i)
                    pairodpout << *i << ' ';
        }

        // measure pairing correlation
        if(meas_paircorr) {
            int i_start = 0;
            int i_end = N;
            if( Npart > 1 ) {
                i_start = ( N/Npart )*(ithpart-1);
                i_end =  ( N/Npart )*ithpart;
                if( ithpart==Npart ) { i_end =  N; }
            }
            println("tot_numer = ", N, " is divided into ", Npart, " part");
            println("i_start = ", i_start, "  i_end = ", i_end);

            std::vector<Cplx> pair1a4( ((i_end*6-i_start*6)*(nnlist.size()*6-i_start*6+nnlist.size()*6-i_end*6+6)/2) );
            std::fill(pair1a4.begin(), pair1a4.end(), Cplx(0.0,0.0));
            std::vector<Cplx> pair2a3( ((i_end*6-i_start*6)*(nnlist.size()*6-i_start*6+nnlist.size()*6-i_end*6+6)/2) );
            std::fill(pair2a3.begin(), pair2a3.end(), Cplx(0.0,0.0));
            int numofcalcorr = 0;
            //println("pair1a4.size = ", pair1a4.size());
            //println("pair2a3.size = ", pair2a3.size());
            //for(int i = 0; i < nnlist.size(); ++i) {
            for(int n1 = i_start; n1 < i_end; ++n1) {
                int i = nnlist[n1].s0;
                int n1y = i%Ny;
                for (int id1 = 0; id1<6; id1++) {
                    int j = nnlist[n1].snn[id1];
                    std::vector< std::pair<int,int> > op34pair_vec ={}; // store (opk,opl) pair
                    std::vector<int> corr_ind = {};  // index in paircorr
                    for(int n2 = n1; n2 < nnlist.size(); ++n2) {
                        int k = nnlist[n2].s0;
                        int n2y = k%Ny;
                        for (int id2 = 0; id2<6; id2++) {
                            int l = nnlist[n2].snn[id2];

                            std::vector<int> sites_tmp = { j, i, k, l };
                            // i,j,k,l, only select i!=j < k!=l
                            // printfln("n1, id1, n2, id2 = %d, %d, %d, %d", n1, id1, n2, id2 );
                            // printfln("j, i, k, l = %d, %d, %d, %d", j, i, k, l );
                            if( (sites_tmp[0] != sites_tmp[1]) &&
                                (sites_tmp[0] <  sites_tmp[2]) &&
                                (sites_tmp[0] <  sites_tmp[3]) &&
                                (sites_tmp[1] <  sites_tmp[2]) &&
                                (sites_tmp[1] <  sites_tmp[3]) &&
                                (sites_tmp[2] != sites_tmp[3]) ) {
                                if(  (!lfixi0 || ((n1 == x0*Ny + n2%Ny) || (n2 == x0*Ny + n1%Ny)) ) && n1y==n2y
                                    && (id1==0 || id1 == 4 || id1 ==5)
                                    && (id2==0 || id2 == 4 || id2 ==5) ) { // fix n1
                                numofcalcorr += 1;
                                op34pair_vec.emplace_back( std::make_pair( sites_tmp[2], sites_tmp[3] ) );
                                int ind_meas = ( nnlist.size()-i_start + nnlist.size()-n1+1)*(n1-i_start)/2*36 + (nnlist.size()-n1)*6*id1 + (n2-n1)*6 + id2;
                                //std::cout << " ind_meas = " << ind_meas <<std::endl;
                                pair1a4[ind_meas] = 0.0; // initial the obserator
                                pair2a3[ind_meas] = 0.0; // initial the obserator
                                corr_ind.emplace_back( ind_meas );
                                // println("j,i<k,l case");
                                //std::cout << " k, l = " << k << " " << l <<std::endl;
                                //std::cout << " ind_meas = " << ind_meas <<std::endl;
                                }
                            } else {
                                ///if( !lfixi0 || ( ((n1 == x0*Ny + n2%Ny) || (n2 == x0*Ny + n1%Ny) ) && 
                                ///                                   (id1==0 || id1 == 4 || id1 ==5) &&
                                ///                                   (id2==0 || id2 == 4 || id2 ==5)) ) { // fix n1
                                if(  (!lfixi0 || ((n1 == x0*Ny + n2%Ny) || (n2 == x0*Ny + n1%Ny)) ) && n1y==n2y
                                    && (id1==0 || id1 == 4 || id1 ==5)
                                    && (id2==0 || id2 == 4 || id2 ==5) ) { // fix n1
                                numofcalcorr += 1;
                                // use mfourbodyf
                                auto pair1a4_tmp = mfourbodyf(psi,sites,sites_tmp,"Adagdn","Adagup","Aup","Adn") + 
                                                   mfourbodyf(psi,sites,sites_tmp,"Adagup","Adagdn","Adn","Aup");
                                auto pair2a3_tmp = mfourbodyf(psi,sites,sites_tmp,"Adagup","Adagdn","Aup","Adn") + 
                                                   mfourbodyf(psi,sites,sites_tmp,"Adagdn","Adagup","Adn","Aup");
                                int ind_meas = ( nnlist.size()-i_start + nnlist.size()-n1+1)*(n1-i_start)/2*36 + (nnlist.size()-n1)*6*id1 + (n2-n1)*6 + id2;
                                // std::cout << " ind_meas = " << ind_meas <<std::endl;
                                pair1a4[ind_meas] = pair1a4_tmp ;
                                pair2a3[ind_meas] = pair2a3_tmp ;
                                // printfln(" %d, %d, %d, %d, paircorr1a4 = %.12f", j, i, k, l, pair1a4_tmp);
                                // printfln(" %d, %d, %d, %d, paircorr2a3 = %.12f", j, i, k, l, pair2a3_tmp);
                                }
                            }
                        }
                    } // for(int n2 = n1; n2 < int(nnlist.size()); ++n2) {
                    if( op34pair_vec.size() > 0 ) {
                        //std::cout << " op34pair_vec = " <<std::endl;
                        //for (auto rr : op34pair_vec ) { std::cout << rr.first <<" "<< rr.second << '\n'; }
                        mfourbodyf_str(psi, sites, {j, i}, "Adagdn", "Adagup", op34pair_vec, "Aup", "Adn", corr_ind, pair1a4, Cplx(1.0,0.0) );
                        mfourbodyf_str(psi, sites, {j, i}, "Adagup", "Adagdn", op34pair_vec, "Adn", "Aup", corr_ind, pair1a4, Cplx(1.0,0.0) );
                        mfourbodyf_str(psi, sites, {j, i}, "Adagup", "Adagdn", op34pair_vec, "Aup", "Adn", corr_ind, pair2a3, Cplx(1.0,0.0) );
                        mfourbodyf_str(psi, sites, {j, i}, "Adagdn", "Adagup", op34pair_vec, "Adn", "Aup", corr_ind, pair2a3, Cplx(1.0,0.0) );
                    }
                    printfln(" paircorr of j,i = %d, %d, k,l=?,? done ", j, i);
                } // for (int id1 = 0; id1<6; id1++) {
            } // for(int n1 = i_start; n1 < i_end; ++n1) {
            println("numofcalcorr = ", numofcalcorr);
            std::vector<Cplx> paircorr;
            for ( int i = 0; i < pair1a4.size(); ++i) {
                paircorr.emplace_back( pair1a4[i] );
                paircorr.emplace_back( pair2a3[i] );
            }
            // output
            std::ofstream paircorrout;
            if( lfixi0 ) {
            paircorrout.open("paircorrfixi0.out",std::ios::out);
            } else {
            paircorrout.open("paircorr.out",std::ios::out);
            }
            paircorrout.precision(12);
            for (std::vector<Cplx>::const_iterator i = paircorr.begin(); i != paircorr.end(); ++i)
                    paircorrout << *i << ' ';
        } // if(meas_paircorr) {

    }  // if(domeas && meas_pair) {

    println( "\nRUNNING FINISHED ^_^ !!! " );
    return 0;
}
