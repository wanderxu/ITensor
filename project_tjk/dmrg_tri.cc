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
    auto meas_dimercorr = input.getYesNo("meas_dimercorr",false);
    auto meas_dxcorr = input.getYesNo("meas_dxcorr",false);
    auto meas_dycorr = input.getYesNo("meas_dycorr",false);
    auto meas_dxycorr = input.getYesNo("meas_dxycorr",false);
    auto meas_chiralcorr = input.getYesNo("meas_chiralcorr",false);
    auto meas_pair = input.getYesNo("meas_pair",false);
    auto meas_pairodp = input.getYesNo("meas_pairodp",false);
    auto meas_paircorr = input.getYesNo("meas_paircorr",false);
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
        auto lattice = triangularLatticev2(Nx,Ny,{"YPeriodic=",yperiodic});
        auto lattice4plaque = triangularLattice4Plaque(Nx,Ny,{"YPeriodic=",yperiodic});

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
        for(int i = N; i >= 1; --i) {
          if(p > i) {
            println("Doubly occupying site ",i);
            state.set(i,"UpDn");
            p -= 2;
          } else if(p > 0) {
            println("Singly occupying site ",i);
            state.set(i,(i%2==1 ? "Up" : "Dn"));
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
        for (std::vector<double>::const_iterator i = Sz_meas.begin(); i != Sz_meas.end(); ++i)
                fSzout << *i << ' ';

        std::ofstream fSpout("Sip.out",std::ios::out);
        for (std::vector<Cplx>::const_iterator i = Sp_meas.begin(); i != Sp_meas.end(); ++i)
                fSpout << *i << ' ';

        std::ofstream fSmout("Sim.out",std::ios::out);
        for (std::vector<Cplx>::const_iterator i = Sm_meas.begin(); i != Sm_meas.end(); ++i)
                fSmout << *i << ' ';

        std::ofstream fSiSjout("SiSj.out",std::ios::out);
        for (std::vector<double>::const_iterator i = SiSj_meas.begin(); i != SiSj_meas.end(); ++i)
                fSiSjout << *i << ' ';

        std::ofstream fSiSjzzout("SiSjzz.out",std::ios::out);
        for (std::vector<double>::const_iterator i = SiSjzz_meas.begin(); i != SiSjzz_meas.end(); ++i)
                fSiSjzzout << *i << ' ';

        std::ofstream fSiSjpmout("SiSjpm.out",std::ios::out);
        for (std::vector<double>::const_iterator i = SiSjpm_meas.begin(); i != SiSjpm_meas.end(); ++i)
                fSiSjpmout << *i << ' ';
    }

    if(domeas && meas_dimercorr) {
        println("\n////////////////////////////");
        println("Start to perform measurement of dimer correlation");
        // 
        // measure dimer correlation
        //
        
        // make the dimer table
        // x-direction
        auto num_x_dimer = (Nx-1)*Ny;
        auto num_y_dimer = Nx*(yperiodic ? Ny : Ny-1);
        auto num_xy_dimer = (Nx-1)*(yperiodic ? Ny : Ny-1);
        LatticeGraph x_dimer;
        LatticeGraph y_dimer;
        LatticeGraph xy_dimer;
        x_dimer.reserve(num_x_dimer);
        y_dimer.reserve(num_y_dimer);
        xy_dimer.reserve(num_xy_dimer);
        for(int n = 1; n <= N; ++n)
            {
            int x = (n-1)/Ny+1;
            int y = (n-1)%Ny+1;

            //X-direction bonds
            if(x < Nx) x_dimer.emplace_back(n,n+Ny);

            if(Ny > 1) //2d bonds
                {
                //vertical bond
                if(y < Ny) y_dimer.emplace_back(n,n+1);
                if((y == Ny) && yperiodic) y_dimer.emplace_back(n,n-Ny+1);

                //Diagonal bonds
                if((x < Nx) && (y < Ny)) xy_dimer.emplace_back(n,n+Ny+1);
                if((x < Nx) && (y == Ny) && yperiodic) xy_dimer.emplace_back(n,n+1);
                }
            }
        if(int(x_dimer.size()) != num_x_dimer) Error("Wrong number of x_dimer");
        if(int(y_dimer.size()) != num_y_dimer) Error("Wrong number of y_dimer");
        if(int(xy_dimer.size()) != num_xy_dimer) Error("Wrong number of xy_dimer");
        println( "\nx_dimer: \n", x_dimer );
        println( "y_dimer: \n", y_dimer );
        println( "xy_dimer: \n", xy_dimer );

        // measure x_dimer correlation
    if( meas_dxcorr) {
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
        std::vector<double> dx_meas={}; // store <Di>
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
                    if(i == j){
                        auto dimer_meas = 0.0;
                        // note conjugation codition is used, assume s1!=s2, and it is really the case here
                        dimer_meas += mtwobody(psi,sites,{x_dimer[i].s1, x_dimer[i].s2}, "S+", "S-");
                        dimer_meas += mtwobody(psi,sites,{x_dimer[i].s1, x_dimer[i].s2}, "Sz", "Sz");
                        dx_meas.emplace_back(dimer_meas);
                    }
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
        for (std::vector<double>::const_iterator i = dxdx_meas.begin(); i != dxdx_meas.end(); ++i)
                fdxdxout << *i << ' ';

        std::ofstream fdxout("Dxi.out",std::ios::out);
        for (std::vector<double>::const_iterator i = dx_meas.begin(); i != dx_meas.end(); ++i)
                fdxout << *i << ' ';
    } // end if( meas_dxcorr) {

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
        std::vector<double> dy_meas={}; // store <Di>
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
                    if(i == j){
                        auto dimer_meas = 0.0;
                        // note conjugation codition is used, assume s1!=s2, and it is really the case here
                        dimer_meas += mtwobody(psi,sites,{y_dimer[i].s1, y_dimer[i].s2}, "S+", "S-");
                        dimer_meas += mtwobody(psi,sites,{y_dimer[i].s1, y_dimer[i].s2}, "Sz", "Sz");
                        dy_meas.emplace_back(dimer_meas);
                    }
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
        for (std::vector<double>::const_iterator i = dydy_meas.begin(); i != dydy_meas.end(); ++i)
                fdydyout << *i << ' ';

        std::ofstream fdyout("Dyi.out",std::ios::out);
        for (std::vector<double>::const_iterator i = dy_meas.begin(); i != dy_meas.end(); ++i)
                fdyout << *i << ' ';
    } // if( meas_dycorr) {

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
        std::vector<double> dxy_meas={}; // store <Di>
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
                    if(i == j){
                        auto dimer_meas = 0.0;
                        // note conjugation codition is used, assume s1!=s2, and it is really the case here
                        dimer_meas += mtwobody(psi,sites,{xy_dimer[i].s1, xy_dimer[i].s2}, "S+", "S-");
                        dimer_meas += mtwobody(psi,sites,{xy_dimer[i].s1, xy_dimer[i].s2}, "Sz", "Sz");
                        dxy_meas.emplace_back(dimer_meas);
                    }
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
        for (std::vector<double>::const_iterator i = dxydxy_meas.begin(); i != dxydxy_meas.end(); ++i)
                fdxydxyout << *i << ' ';

        std::ofstream fdxyout("Dxyi.out",std::ios::out);
        for (std::vector<double>::const_iterator i = dxy_meas.begin(); i != dxy_meas.end(); ++i)
                fdxyout << *i << ' ';
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

            if((x < Nx) && (y < Ny)) {
                tri_plaq.emplace_back(n, n+Ny+1, n+Ny); // x-direction plaq
                tri_plaq.emplace_back(n, n+1, n+Ny+1); // y-direction plaq
            }
            if((x < Nx) && (y == Ny) && yperiodic) {
                tri_plaq.emplace_back(n, n+1, n+Ny);
                tri_plaq.emplace_back(n, n-Ny+1, n+1);
            }

        }
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
        for (std::vector<double>::const_iterator i = XiXj_meas.begin(); i != XiXj_meas.end(); ++i)
                fXiXjout << *i << ' ';

        std::ofstream fXiout("Xi.out",std::ios::out);
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
            int id1 = (x<Nx ? n+Ny: n+Ny-N); // (1,0) dir
            int id2 = (y>1 ? n-1: n-1+Ny); // (0,-1) dir
            int id3 = (x>1 ? (y>1 ? n-Ny-1 : n-Ny-1+Ny) : ( y>1 ? n-Ny-1+N : n-Ny-1+Ny+N ) ); // (-1,-1) dir
            int id4 = (x>1 ? n-Ny : n-Ny+N); // (-1,0) dir
            int id5 = (y<Ny ? n+1: n+1-Ny); // (0,1) dir
            int id6 = (x<Nx ? (y<Ny ? n+Ny+1 : n+Ny+1-Ny) : ( y<Ny ? n+Ny+1-N : n+Ny+1-N-Ny) ); // (1,1) dir
            std::vector<int> nntmp = {id1, id2, id3, id4, id5, id6};
            nnlist.emplace_back( id0, nntmp );
        }
        println( "\nnnlist: \n", nnlist );

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
                        printfln(" %d, %d, pair = %.12f, %.12f\n", i,j, cpair.cplx());

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
                        printfln(" %d, %d, pair = %.12f, %.12f\n", i,j, -cpair.cplx());
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
                        printfln(" %d, %d, pair = %.12f, %.12f\n", i,j, cpair.cplx());

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
                        printfln(" %d, %d, pair = %.12f, %.12f\n", i,j, -cpair.cplx());

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
            std::vector<Cplx> paircorr={};
            for(int n1 = 1; n1 <= N ; ++n1) {
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
                            // i,j < k,l
                            // printfln(" j,i,k,l = %d, %d, %d, %d ", j, i, k, l);
                            if( i<k && i<l && j<k && j<l) {
                                IQTensor cpairdu; 
                                IQTensor cpairud; 
                                if(i<j) {
                                    // i<j
                                    //   c^+_jdn c^+_iup =-c^+_iup c^+_jdn
                                    // =-a^+_iup F_i ... F_{j-1} a^+_jdn F_j
                                    psi.position(i);
                                    cpairdu = noprime(psi.A(i)*sites.op("Adagup",i), Site);
                                    cpairdu *= sites.op("F",i);
                                    auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);
                                    cpairdu *= dag(prime(psi.A(i),Site,ir));
                                    for( int i1 = i+1; i1<j; ++i1 ) {
                                        cpairdu *= psi.A(i1);
                                        cpairdu *= sites.op("F",i1);
                                        cpairdu *= dag(prime(psi.A(i1),Site,Link));
                                    }
                                    cpairdu *= psi.A(j);
                                    cpairdu = noprime(cpairdu*sites.op("Adagdn",j), Site);
                                    cpairdu *= sites.op("F",j);
                                    cpairdu *= dag(prime(psi.A(j),Site,Link));

                                    //   c^+_jup c^+_idn =-c^+_idn c^+_jup
                                    // =-a^+_idn F_{i+1} ... F_{j-1} a^+_jup
                                    // psi.position(i);
                                    cpairud = psi.A(i)*sites.op("Adagdn",i);
                                    ir = commonIndex(psi.A(i),psi.A(i+1),Link);
                                    cpairud *= dag(prime(psi.A(i),Site,ir));
                                    for( int i1 = i+1; i1<j; ++i1 ) {
                                        cpairud *= psi.A(i1);
                                        cpairud *= sites.op("F",i1);
                                        cpairud *= dag(prime(psi.A(i1),Site,Link));
                                    }
                                    cpairud *= psi.A(j);
                                    cpairud *= sites.op("Adagup",j);
                                    cpairud *= dag(prime(psi.A(j),Site,Link));
                                } else if(i>j) {
                                    // j<i
                                    //   c^+_jdn c^+_iup
                                    // = a^+_jdn F_{j+1} ... F_{i-1} a^+_iup
                                    psi.position(j);
                                    cpairdu = psi.A(j)*sites.op("Adagdn",j);
                                    auto ir = commonIndex(psi.A(j),psi.A(j+1),Link);
                                    cpairdu *= dag(prime(psi.A(j),Site,ir));
                                    for( int j1 = j+1; j1<i; ++j1 ) {
                                        cpairdu *= psi.A(j1);
                                        cpairdu *= sites.op("F",j1);
                                        cpairdu *= dag(prime(psi.A(j1),Site,Link));
                                    }
                                    cpairdu *= psi.A(i);
                                    cpairdu *= sites.op("Adagup",i);
                                    cpairdu *= dag(prime(psi.A(i),Site,Link));

                                    //   c^+_jup c^+_idn
                                    // = a^+_jup F_j ... F_{i-1} a^+_idn F_i
                                    // psi.position(j);
                                    cpairud = noprime(psi.A(j)*sites.op("Adagup",j), Site);
                                    cpairud *= sites.op("F",j);
                                    ir = commonIndex(psi.A(j),psi.A(j+1),Link);
                                    cpairud *= dag(prime(psi.A(j),Site,ir));
                                    for( int j1 = j+1; j1<i; ++j1 ) {
                                        cpairud *= psi.A(j1);
                                        cpairud *= sites.op("F",j1);
                                        cpairud *= dag(prime(psi.A(j1),Site,Link));
                                    }
                                    cpairud *= psi.A(i);
                                    cpairud = noprime(cpairud*sites.op("Adagdn",i), Site);
                                    cpairud *= sites.op("F",i);
                                    cpairud *= dag(prime(psi.A(i),Site,Link));
                                } else { Error("i should not equal to j !"); }

                                auto cpairduud = cpairdu; // 1
                                auto cpairdudu = cpairdu; // 3
                                auto cpairudud = cpairud; // 2
                                auto cpairuddu = cpairud; // 4
                                for (int j1 = std::max(i,j)+1; j1<std::min(k,l); ++j1) {
                                    cpairduud *= psi.A(j1);
                                    cpairdudu *= psi.A(j1);
                                    cpairudud *= psi.A(j1);
                                    cpairuddu *= psi.A(j1);
                                    cpairduud *= dag(prime(psi.A(j1),Link));
                                    cpairdudu *= dag(prime(psi.A(j1),Link));
                                    cpairudud *= dag(prime(psi.A(j1),Link));
                                    cpairuddu *= dag(prime(psi.A(j1),Link));
                                }
                                if(k<l) {
                                    // k<l
                                    //   c_kup c_ldn
                                    // = a_kup F_k ... F_l a_ldn
                                    cpairduud *= psi.A(k); // 1
                                    cpairudud *= psi.A(k); // 2
                                    cpairduud = noprime(cpairduud*sites.op("Aup",k), Site);
                                    cpairudud = noprime(cpairudud*sites.op("Aup",k), Site);
                                    cpairduud *= sites.op("F",k);
                                    cpairudud *= sites.op("F",k);
                                    cpairduud *= dag(prime(psi.A(k),Site,Link));
                                    cpairudud *= dag(prime(psi.A(k),Site,Link));
                                    for (int k1 = k+1; k1<l; ++k1) {
                                        cpairduud *= psi.A(k1);
                                        cpairudud *= psi.A(k1);
                                        cpairduud *= sites.op("F",k1);
                                        cpairudud *= sites.op("F",k1);
                                        cpairduud *= dag(prime(psi.A(k1),Site,Link));
                                        cpairudud *= dag(prime(psi.A(k1),Site,Link));
                                    }
                                    cpairduud *= psi.A(l);
                                    cpairudud *= psi.A(l);
                                    cpairduud = noprime(cpairduud*sites.op("F",l),Site);
                                    cpairudud = noprime(cpairudud*sites.op("F",l),Site);
                                    cpairduud *= sites.op("Adn",l);
                                    cpairudud *= sites.op("Adn",l);
                                    auto ir = commonIndex(psi.A(l),psi.A(l-1),Link);
                                    cpairduud *= dag(prime(psi.A(l),Site,ir));
                                    cpairudud *= dag(prime(psi.A(l),Site,ir));

                                    //   c_kdn c_lup
                                    // = -a_kdn F_{k+1} ... F_{l-1} a_lup
                                    cpairdudu *= psi.A(k); // 3
                                    cpairuddu *= psi.A(k); // 4
                                    cpairdudu *= sites.op("Adn",k);
                                    cpairuddu *= sites.op("Adn",k);
                                    cpairdudu *= dag(prime(psi.A(k),Site,Link));
                                    cpairuddu *= dag(prime(psi.A(k),Site,Link));
                                    for (int k1 = k+1; k1<l; ++k1) {
                                        cpairdudu *= psi.A(k1);
                                        cpairuddu *= psi.A(k1);
                                        cpairdudu *= sites.op("F",k1);
                                        cpairuddu *= sites.op("F",k1);
                                        cpairdudu *= dag(prime(psi.A(k1),Site,Link));
                                        cpairuddu *= dag(prime(psi.A(k1),Site,Link));
                                    }
                                    cpairdudu *= psi.A(l);
                                    cpairuddu *= psi.A(l);
                                    cpairdudu *= sites.op("Aup",l);
                                    cpairuddu *= sites.op("Aup",l);
                                    ir = commonIndex(psi.A(l),psi.A(l-1),Link);
                                    cpairdudu *= dag(prime(psi.A(l),Site,ir));
                                    cpairuddu *= dag(prime(psi.A(l),Site,ir));
                                } else if (k>l ) {
                                    // l<k
                                    //   c_kup c_ldn = - c_ldn c_kup
                                    // = a_ldn F_{l+1} ... F_{k-1} a_kup // sign is cancled
                                    cpairduud *= psi.A(l); // 1
                                    cpairudud *= psi.A(l); // 2
                                    cpairduud *= sites.op("Adn",l);
                                    cpairudud *= sites.op("Adn",l);
                                    cpairduud *= dag(prime(psi.A(l),Site,Link));
                                    cpairudud *= dag(prime(psi.A(l),Site,Link));
                                    for (int l1 = l+1; l1<k; ++l1) {
                                        cpairduud *= psi.A(l1);
                                        cpairudud *= psi.A(l1);
                                        cpairduud *= sites.op("F",l1);
                                        cpairudud *= sites.op("F",l1);
                                        cpairduud *= dag(prime(psi.A(l1),Site,Link));
                                        cpairudud *= dag(prime(psi.A(l1),Site,Link));
                                    }
                                    cpairduud *= psi.A(k);
                                    cpairudud *= psi.A(k);
                                    cpairduud *= sites.op("Aup",k);
                                    cpairudud *= sites.op("Aup",k);
                                    auto ir = commonIndex(psi.A(k),psi.A(k-1),Link);
                                    cpairduud *= dag(prime(psi.A(k),Site,ir));
                                    cpairudud *= dag(prime(psi.A(k),Site,ir));

                                    //   c_kdn c_lup = - c_lup c_kdn
                                    // = -a_lup F_l ... F_k a_kdn
                                    cpairdudu *= psi.A(l); // 3
                                    cpairuddu *= psi.A(l); // 4
                                    cpairdudu = noprime(cpairdudu*sites.op("Aup",l), Site);
                                    cpairuddu = noprime(cpairuddu*sites.op("Aup",l), Site);
                                    cpairdudu *= sites.op("F",l);
                                    cpairuddu *= sites.op("F",l);
                                    cpairdudu *= dag(prime(psi.A(l),Site,Link));
                                    cpairuddu *= dag(prime(psi.A(l),Site,Link));
                                    for (int l1 = l+1; l1<k; ++l1) {
                                        cpairdudu *= psi.A(l1);
                                        cpairuddu *= psi.A(l1);
                                        cpairdudu *= sites.op("F",l1);
                                        cpairuddu *= sites.op("F",l1);
                                        cpairdudu *= dag(prime(psi.A(l1),Site,Link));
                                        cpairuddu *= dag(prime(psi.A(l1),Site,Link));
                                    }
                                    cpairdudu *= psi.A(k);
                                    cpairuddu *= psi.A(k);
                                    cpairdudu = noprime(cpairdudu*sites.op("F",k), Site);
                                    cpairuddu = noprime(cpairuddu*sites.op("F",k), Site);
                                    cpairdudu *= sites.op("Adn",k);
                                    cpairuddu *= sites.op("Adn",k);
                                    ir = commonIndex(psi.A(k),psi.A(k-1),Link);
                                    cpairdudu *= dag(prime(psi.A(k),Site,ir));
                                    cpairuddu *= dag(prime(psi.A(k),Site,ir));
                                } else { Error("k should not equal to l !"); }
                                // store and output
                                Cplx pair1a4;
                                Cplx pair2a3;
                                // consider the sign
                                if ( i>j ) {
                                    pair1a4 = cpairduud.cplx() - cpairuddu.cplx(); // 1+4
                                    pair2a3 = cpairudud.cplx() - cpairdudu.cplx(); // 2+3
                                } else {
                                    pair1a4 = -cpairduud.cplx() + cpairuddu.cplx(); // 1+4
                                    pair2a3 = -cpairudud.cplx() + cpairdudu.cplx(); // 2+3
                                }
                                paircorr.emplace_back( pair1a4  );
                                paircorr.emplace_back( pair2a3  );

                                std::vector<int> sites_tmp = { j, i, k, l };
                                auto pair1a4_tmp = mfourbodyf(psi,sites,sites_tmp,"Adagdn","Adagup","Aup","Adn") + 
                                                   mfourbodyf(psi,sites,sites_tmp,"Adagup","Adagdn","Adn","Aup");
                                auto pair2a3_tmp = mfourbodyf(psi,sites,sites_tmp,"Adagup","Adagdn","Aup","Adn") + 
                                                   mfourbodyf(psi,sites,sites_tmp,"Adagdn","Adagup","Adn","Aup");

                                printfln(" %d, %d, %d, %d, paircorr1a4 = %.12f, %.12f\n", j, i, k, l, pair1a4);
                                printfln(" %d, %d, %d, %d, pair1a4_tmp = %.12f, %.12f\n", j, i, k, l, pair1a4_tmp);
                                printfln(" %d, %d, %d, %d, paircorr2a3 = %.12f, %.12f\n", j, i, k, l, pair2a3);
                                printfln(" %d, %d, %d, %d, pair2a3_tmp = %.12f, %.12f\n", j, i, k, l, pair2a3_tmp);
                            } // if( i<k && i<l && j<k && j<l) {
                            else {
                                // use AutoMPO
                                auto ampo = AutoMPO(sites);
                            }
                        }
                    }
                } // for (int id1 = 0; id1<6; id1++) {
            } // for(int n1 = 1; n1 <= N ; ++n1) {
            // output
            std::ofstream paircorrout("paircorr.out",std::ios::out);
            paircorrout.precision(12);
            for (std::vector<Cplx>::const_iterator i = paircorr.begin(); i != paircorr.end(); ++i)
                    paircorrout << *i << ' ';
        } // if(meas_paircorr) {
    }

    //// test fourbody
    //srand (time(NULL));
    //std::vector<int> sites_tmp={rand()%N+1,rand()%N+1,rand()%N+1,rand()%N+1};
    //sites_tmp = {1,7,8,14};
    //for (auto n : sites_tmp ) { std::cout << n <<" "; }
    //std::cout << '\n';
    //println( " <sijkl> = ", mfourbody(psi,sites,sites_tmp,"S+","S-","S+","S-") );

    //std::vector< std::pair<int,int> > op34pair_vec ={};
    //std::vector<int> corr_ind = {};
    //std::vector<double> dxdx_meas( 1 );
    //dxdx_meas[0]=0.0;
    //op34pair_vec.emplace_back( std::make_pair( sites_tmp[2], sites_tmp[3] ) );
    //corr_ind.emplace_back(0);
    //mfourbody_str(psi, sites, {sites_tmp[0],sites_tmp[1]}, "S+", "S-", op34pair_vec, "S+", "S-", corr_ind, dxdx_meas,1.0);
    //println( " mfourbody_str <sijkl> = ");
    //for (auto rr : dxdx_meas ) { std::cout << rr <<" "; }
    //std::cout << '\n';

    ////// test sixbody
    ////std::vector<int> sites_tmp = {1,2,3,4,5,6};
    ////for (auto n : sites_tmp ) { std::cout << n <<" "; }
    ////std::cout << '\n';
    ////println( "with msixbody <sijklmn> = ", msixbody(psi,sites,sites_tmp,"S+","S-","S+","S-","Sz","Sz") );

    ////auto tmp_mpo = AutoMPO(sites);
    ////tmp_mpo += 1.0,"S+",sites_tmp[0], "S-", sites_tmp[1],"S+",sites_tmp[2], "S-", sites_tmp[3],"Sz",sites_tmp[4], "Sz", sites_tmp[5];
    ////auto tmp_corr = IQMPO(tmp_mpo);
    ////println( "with overlap <sijklmn> =", overlap(psi,tmp_corr,psi));

    ////////sites_tmp = {3,3,3,3};
    ////////for (auto n : sites_tmp ) { std::cout << n <<" "; }
    ////////std::cout << '\n';
    ////////println( "with mfourbody <sijkl> = ", mfourbody(psi,sites,sites_tmp,"S+","S-","S+","S-") );

    ////////tmp_mpo = AutoMPO(sites);
    ////////tmp_mpo += 1.0,"S+",sites_tmp[0], "S-", sites_tmp[1],"S+",sites_tmp[2], "S-", sites_tmp[3];
    ////////tmp_corr = IQMPO(tmp_mpo);
    ////////println( "with overlap <sijkl> =", overlap(psi,tmp_corr,psi));

    ////std::vector<double> tmp_meas(1);
    ////std::vector<int> corr_ind = {};
    ////tmp_meas[0]=0.0;
    ////corr_ind.emplace_back(0);
    ////msixbody_str(psi, sites, {sites_tmp[0],sites_tmp[1],sites_tmp[2]}, "S+", "S-", "S+",
    ////                           { {sites_tmp[3],sites_tmp[4],sites_tmp[5]} }, "S-", "Sz", "Sz", corr_ind, tmp_meas,1.0);
    ////println( "with msixbody_str <sijklmn> =", tmp_meas[0]);

    println( "\nRUNNING FINISHED ^_^ !!! " );
    return 0;
    }
