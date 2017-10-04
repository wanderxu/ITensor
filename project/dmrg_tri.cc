#include "itensor/all.h"
#include "triangular_more.h"

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
    double J1 = input.getReal("J1");
    double J2 = input.getReal("J2");
    double gamma1 = input.getReal("gamma1");
    double gamma2 = input.getReal("gamma2");
    auto readmps = input.getYesNo("readmps",false);
    auto eneropt = input.getYesNo("eneropt",true);
    auto domeas = input.getYesNo("domeas",false);
    auto meas_spincorr = input.getYesNo("meas_spincorr",false);
    auto meas_dimercorr = input.getYesNo("meas_dimercorr",false);
    auto meas_chiralcorr = input.getYesNo("meas_chiralcorr",false);
    auto quiet = input.getYesNo("quiet",true);

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

    //Create the sweeps class & print
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    SpinHalf sites;
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
        //SpinHalf sites;
        readFromFile("sites_file", sites);
        //IQMPS psi(sites);
        psi=IQMPS(sites);
        readFromFile("psi_file", psi);
        //IQMPO H(sites);
        H=IQMPO(sites);
        readFromFile("H_file", H);
        auto psiHpsi = overlap(psi,H,psi);
        //auto psiHHpsi = overlap(psi,H,H,psi);
        printfln(" Intial energy information from input file: ");
        printfln("\n<psi|H|psi> = %.10f", psiHpsi );
        //printfln("\n<psi|H^2|psi> = %.10f", psiHHpsi );
        //printfln("\n<psi|H^2|psi> - <psi|H|psi>^2 = %.10f", psiHHpsi-psiHpsi*psiHpsi );
        println("\nTotal QN of Ground State = ",totalQN(psi));
    }
    else {
        println("\n//////////////////////////////////////////////////");
        println("Building basis, wavefunction and H from scratch ......\n");
        //
        // Initialize the site degrees of freedom.
        //
        //auto sites = SpinHalf(N);
        sites = SpinHalf(N);

        //
        // Use the AutoMPO feature to create the 
        // nearest-neighbor XXZ model with ring exchange.
        //

        auto ampo = AutoMPO(sites);
        auto lattice = triangularLattice(Nx,Ny,{"YPeriodic=",yperiodic});
        auto lattice4plaque = triangularLattice4Plaque(Nx,Ny,{"YPeriodic=",yperiodic});

        println("H is made up of ");
        println("\nBound:\n", lattice);
        println("Total number of nn bound: ", lattice.size());

        println("\nPlaque:\n", lattice4plaque);
        println("Total number of plaques: ", lattice4plaque.size());

        // two-body term, nearest neighbor
        for(auto bnd : lattice)
            {
            ampo += J1*2.0,"S+",bnd.s1,"S-",bnd.s2;
            ampo += J1*2.0,"S-",bnd.s1,"S+",bnd.s2;
            ampo += J1*gamma1*4.0,"Sz",bnd.s1,"Sz",bnd.s2;
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
        auto energy = dmrg(psi,H,sweeps,{"Quiet",quiet});

        //
        // Print the final energy reported by DMRG
        //
        printfln("\nGround State Energy = %.10f", energy);

        auto psiHpsi = overlap(psi,H,psi);
        auto psiHHpsi = overlap(psi,H,H,psi);
        printfln("\n<psi|H|psi> = %.10f", psiHpsi );
        printfln("\n<psi|H^2|psi> = %.10f", psiHHpsi );
        printfln("\n<psi|H^2|psi> - <psi|H|psi>^2 = %.10f", psiHHpsi-psiHpsi*psiHpsi );

        printfln("\n<psi|H|psi> / N = %.10f", psiHpsi/N );
        printfln("\n<psi|H^2|psi> / N^2 = %.10f", psiHHpsi/(N*N) );
        printfln("\nsqrt( <psi|H^2|psi> - <psi|H|psi>^2 ) / N = %.10f", sqrt(psiHHpsi-psiHpsi*psiHpsi)/N );

        println("\nTotal QN of Ground State = ",totalQN(psi));

        // after the MPS converged, write basis, psi, and H to disk
        writeToFile("sites_file", sites);
        writeToFile("psi_file", psi);
        writeToFile("H_file", H);
    }

    if(domeas && meas_spincorr) {
        println("\n////////////////////////////");
        println("Start to perform measurement of spin correlation\n");
        //
        // Measure Si.Sj of every {i,j}, and total M
        //
        auto totalM = 0.0;
        std::vector<double> SiSj_meas={};
        std::vector<double> Sz_meas={};
        std::vector<double> Sp_meas={};
        std::vector<double> Sm_meas={};
        for ( int i = 1; i <= N; ++i ) {
            //'gauge' the MPS to site i
            psi.position(i); 
            
            //psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation

            // i == j part

            // magnetization
            auto ket = psi.A(i);
            auto bra = dag(prime(ket,Site));
            auto sz_tmp = (bra*sites.op("Sz",i)*ket).real();
            Sz_meas.emplace_back(sz_tmp);
            totalM +=  sz_tmp;

            auto sp_tmp = (bra*sites.op("S+",i)*ket).real();
            Sp_meas.emplace_back(sp_tmp);
            auto sm_tmp = (bra*sites.op("S-",i)*ket).real();
            Sm_meas.emplace_back(sm_tmp);

            auto ss_tmp = 0.0;
            ss_tmp += 0.75*((dag(ket)*ket).real());
            SiSj_meas.emplace_back(ss_tmp);
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
                    ss_tmp += 0.5*( (Cpm*op_jm)*dag(prime(psi.A(j),jl,Site)) ).real();

                    auto op_jp = sites.op("S+",j);
                    ss_tmp += 0.5*( (Cmp*op_jp)*dag(prime(psi.A(j),jl,Site)) ).real();

                    auto op_jz = sites.op("Sz",j);
                    ss_tmp += ( (Czz*op_jz)*dag(prime(psi.A(j),jl,Site)) ).real();
                    
                    SiSj_meas.emplace_back(ss_tmp);
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

        std::ofstream fSzout("Siz.out",std::ios::out);
        for (std::vector<double>::const_iterator i = Sz_meas.begin(); i != Sz_meas.end(); ++i)
                fSzout << *i << ' ';

        std::ofstream fSpout("Sip.out",std::ios::out);
        for (std::vector<double>::const_iterator i = Sp_meas.begin(); i != Sp_meas.end(); ++i)
                fSpout << *i << ' ';

        std::ofstream fSmout("Sim.out",std::ios::out);
        for (std::vector<double>::const_iterator i = Sm_meas.begin(); i != Sm_meas.end(); ++i)
                fSmout << *i << ' ';

        std::ofstream fSiSjout("SiSj.out",std::ios::out);
        for (std::vector<double>::const_iterator i = SiSj_meas.begin(); i != SiSj_meas.end(); ++i)
                fSiSjout << *i << ' ';
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


        auto DDmpo = AutoMPO(sites);
        auto DDcorr = IQMPO(DDmpo);
        auto ddcorr_meas = overlap(psi,DDcorr,psi);
        // measure x_dimer correlation
        println("measure x_dimer correlation");
        std::vector<double> dxdx_meas={};
        std::vector<double> dx_meas={};
        for(int i = 0; i < int(x_dimer.size()); ++i) {
            for(int j = i; j < int(x_dimer.size()); ++j) {
                std::vector<int> sites_tmp = { x_dimer[i].s1, x_dimer[i].s2, x_dimer[j].s1, x_dimer[j].s2 };
                //std::sort( sites_tmp.begin(), sites_tmp.end() ); // sort the pair
                //println( sites_tmp );
                for (auto n : sites_tmp ) { std::cout << n <<" "; }
                std::cout << '\n';

                // calculate correlation, Si*Sj*Sk*Sl
                if(i == j){
                    DDmpo = AutoMPO(sites);
                    DDmpo += 0.5,"S+",x_dimer[i].s1,"S-",x_dimer[i].s2;
                    DDmpo += 0.5,"S-",x_dimer[i].s1,"S+",x_dimer[i].s2;
                    DDmpo += 1.0,"Sz",x_dimer[i].s1,"Sz",x_dimer[i].s2;
                    DDcorr = IQMPO(DDmpo);
                    ddcorr_meas = overlap(psi,DDcorr,DDcorr,psi);
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas); 
                    dxdx_meas.emplace_back(ddcorr_meas);
                    dx_meas.emplace_back(overlap(psi,DDcorr,psi));
                }
                else {
                    DDmpo = AutoMPO(sites);
                    DDmpo += 0.25,"S+",x_dimer[i].s1,"S-",x_dimer[i].s2,"S+",x_dimer[j].s1,"S-",x_dimer[j].s2;
                    DDmpo += 0.25,"S+",x_dimer[i].s1,"S-",x_dimer[i].s2,"S-",x_dimer[j].s1,"S+",x_dimer[j].s2;
                    DDmpo += 0.50,"S+",x_dimer[i].s1,"S-",x_dimer[i].s2,"Sz",x_dimer[j].s1,"Sz",x_dimer[j].s2;
                    DDmpo += 0.25,"S-",x_dimer[i].s1,"S+",x_dimer[i].s2,"S+",x_dimer[j].s1,"S-",x_dimer[j].s2;
                    DDmpo += 0.25,"S-",x_dimer[i].s1,"S+",x_dimer[i].s2,"S-",x_dimer[j].s1,"S+",x_dimer[j].s2;
                    DDmpo += 0.50,"S-",x_dimer[i].s1,"S+",x_dimer[i].s2,"Sz",x_dimer[j].s1,"Sz",x_dimer[j].s2;
                    DDmpo += 0.50,"Sz",x_dimer[i].s1,"Sz",x_dimer[i].s2,"S+",x_dimer[j].s1,"S-",x_dimer[j].s2;
                    DDmpo += 0.50,"Sz",x_dimer[i].s1,"Sz",x_dimer[i].s2,"S-",x_dimer[j].s1,"S+",x_dimer[j].s2;
                    DDmpo += 1.00,"Sz",x_dimer[i].s1,"Sz",x_dimer[i].s2,"Sz",x_dimer[j].s1,"Sz",x_dimer[j].s2;
                    DDcorr = IQMPO(DDmpo);
                    ddcorr_meas = overlap(psi,DDcorr,psi);
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas);
                    dxdx_meas.emplace_back(ddcorr_meas);
                }
            }
        }
        std::ofstream fdxdxout("DxiDxj.out",std::ios::out);
        for (std::vector<double>::const_iterator i = dxdx_meas.begin(); i != dxdx_meas.end(); ++i)
                fdxdxout << *i << ' ';

        std::ofstream fdxout("Dxi.out",std::ios::out);
        for (std::vector<double>::const_iterator i = dx_meas.begin(); i != dx_meas.end(); ++i)
                fdxout << *i << ' ';

        // measure y_dimer correlation
        // note when we have yperiodic boundary condition, care should be take for y_dimer correlation
        auto Dmpoi = AutoMPO(sites);
        auto Dmpoj = AutoMPO(sites);
        auto Diop = IQMPO(Dmpoi);
        auto Djop = IQMPO(Dmpoj);
        println("measure y_dimer correlation");
        std::vector<double> dydy_meas={};
        std::vector<double> dy_meas={};
        for(int i = 0; i < int(y_dimer.size()); ++i) {
            for(int j = i; j < int(y_dimer.size()); ++j) {
                std::vector<int> sites_tmp = { y_dimer[i].s1, y_dimer[i].s2, y_dimer[j].s1, y_dimer[j].s2 };
                //std::sort( sites_tmp.begin(), sites_tmp.end() ); // sort the pair
                //println( sites_tmp );
                for (auto n : sites_tmp ) { std::cout << n <<" "; }
                std::cout << '\n';

                // calculate correlation, Si*Sj*Sk*Sl
                if(i == j){
                    DDmpo = AutoMPO(sites);
                    DDmpo += 0.5,"S+",y_dimer[i].s1,"S-",y_dimer[i].s2;
                    DDmpo += 0.5,"S-",y_dimer[i].s1,"S+",y_dimer[i].s2;
                    DDmpo += 1.0,"Sz",y_dimer[i].s1,"Sz",y_dimer[i].s2;
                    DDcorr = IQMPO(DDmpo);
                    ddcorr_meas = overlap(psi,DDcorr,DDcorr,psi);
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas); 
                    dydy_meas.emplace_back(ddcorr_meas);
                    dy_meas.emplace_back(overlap(psi,DDcorr,psi));
                }
                else if( yperiodic && ((i+2)%Ny==0) && ((j+1)%Ny==0)) { //this is a specical case accross the y-boundary
                    Dmpoi = AutoMPO(sites);
                    Dmpoi += 0.5,"S+",y_dimer[i].s1,"S-",y_dimer[i].s2;
                    Dmpoi += 0.5,"S-",y_dimer[i].s1,"S+",y_dimer[i].s2;
                    Dmpoi += 1.0,"Sz",y_dimer[i].s1,"Sz",y_dimer[i].s2;
                    Dmpoj = AutoMPO(sites);
                    Dmpoj += 0.5,"S+",y_dimer[j].s1,"S-",y_dimer[j].s2;
                    Dmpoj += 0.5,"S-",y_dimer[j].s1,"S+",y_dimer[j].s2;
                    Dmpoj += 1.0,"Sz",y_dimer[j].s1,"Sz",y_dimer[j].s2;
                    Diop = IQMPO(Dmpoi);
                    Djop = IQMPO(Dmpoj);
                    ddcorr_meas = overlap(psi,Diop,Djop,psi);
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas); 
                    dydy_meas.emplace_back(ddcorr_meas);
                }
                else {
                    DDmpo = AutoMPO(sites);
                    DDmpo += 0.25,"S+",y_dimer[i].s1,"S-",y_dimer[i].s2,"S+",y_dimer[j].s1,"S-",y_dimer[j].s2;
                    DDmpo += 0.25,"S+",y_dimer[i].s1,"S-",y_dimer[i].s2,"S-",y_dimer[j].s1,"S+",y_dimer[j].s2;
                    DDmpo += 0.50,"S+",y_dimer[i].s1,"S-",y_dimer[i].s2,"Sz",y_dimer[j].s1,"Sz",y_dimer[j].s2;
                    DDmpo += 0.25,"S-",y_dimer[i].s1,"S+",y_dimer[i].s2,"S+",y_dimer[j].s1,"S-",y_dimer[j].s2;
                    DDmpo += 0.25,"S-",y_dimer[i].s1,"S+",y_dimer[i].s2,"S-",y_dimer[j].s1,"S+",y_dimer[j].s2;
                    DDmpo += 0.50,"S-",y_dimer[i].s1,"S+",y_dimer[i].s2,"Sz",y_dimer[j].s1,"Sz",y_dimer[j].s2;
                    DDmpo += 0.50,"Sz",y_dimer[i].s1,"Sz",y_dimer[i].s2,"S+",y_dimer[j].s1,"S-",y_dimer[j].s2;
                    DDmpo += 0.50,"Sz",y_dimer[i].s1,"Sz",y_dimer[i].s2,"S-",y_dimer[j].s1,"S+",y_dimer[j].s2;
                    DDmpo += 1.00,"Sz",y_dimer[i].s1,"Sz",y_dimer[i].s2,"Sz",y_dimer[j].s1,"Sz",y_dimer[j].s2;
                    DDcorr = IQMPO(DDmpo);
                    ddcorr_meas = overlap(psi,DDcorr,psi);
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas);
                    dydy_meas.emplace_back(ddcorr_meas);
                }
            }
        }
        std::ofstream fdydyout("DyiDyj.out",std::ios::out);
        for (std::vector<double>::const_iterator i = dydy_meas.begin(); i != dydy_meas.end(); ++i)
                fdydyout << *i << ' ';

        std::ofstream fdyout("Dyi.out",std::ios::out);
        for (std::vector<double>::const_iterator i = dy_meas.begin(); i != dy_meas.end(); ++i)
                fdyout << *i << ' ';

        // measure xy_dimer correlation
        println("measure xy_dimer correlation");
        std::vector<double> dxydxy_meas={};
        std::vector<double> dxy_meas={};
        for(int i = 0; i < int(xy_dimer.size()); ++i) {
            for(int j = i; j < int(xy_dimer.size()); ++j) {
                std::vector<int> sites_tmp = { xy_dimer[i].s1, xy_dimer[i].s2, xy_dimer[j].s1, xy_dimer[j].s2 };
                //std::sort( sites_tmp.begin(), sites_tmp.end() ); // sort the pair
                //println( sites_tmp );
                for (auto n : sites_tmp ) { std::cout << n <<" "; }
                std::cout << '\n';

                // calculate correlation, Si*Sj*Sk*Sl
                if(i == j){
                    DDmpo = AutoMPO(sites);
                    DDmpo += 0.5,"S+",xy_dimer[i].s1,"S-",xy_dimer[i].s2;
                    DDmpo += 0.5,"S-",xy_dimer[i].s1,"S+",xy_dimer[i].s2;
                    DDmpo += 1.0,"Sz",xy_dimer[i].s1,"Sz",xy_dimer[i].s2;
                    DDcorr = IQMPO(DDmpo);
                    ddcorr_meas = overlap(psi,DDcorr,DDcorr,psi);
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas); 
                    dxydxy_meas.emplace_back(ddcorr_meas);
                    dxy_meas.emplace_back(overlap(psi,DDcorr,psi));
                }
                else {
                    DDmpo = AutoMPO(sites);
                    DDmpo += 0.25,"S+",xy_dimer[i].s1,"S-",xy_dimer[i].s2,"S+",xy_dimer[j].s1,"S-",xy_dimer[j].s2;
                    DDmpo += 0.25,"S+",xy_dimer[i].s1,"S-",xy_dimer[i].s2,"S-",xy_dimer[j].s1,"S+",xy_dimer[j].s2;
                    DDmpo += 0.50,"S+",xy_dimer[i].s1,"S-",xy_dimer[i].s2,"Sz",xy_dimer[j].s1,"Sz",xy_dimer[j].s2;
                    DDmpo += 0.25,"S-",xy_dimer[i].s1,"S+",xy_dimer[i].s2,"S+",xy_dimer[j].s1,"S-",xy_dimer[j].s2;
                    DDmpo += 0.25,"S-",xy_dimer[i].s1,"S+",xy_dimer[i].s2,"S-",xy_dimer[j].s1,"S+",xy_dimer[j].s2;
                    DDmpo += 0.50,"S-",xy_dimer[i].s1,"S+",xy_dimer[i].s2,"Sz",xy_dimer[j].s1,"Sz",xy_dimer[j].s2;
                    DDmpo += 0.50,"Sz",xy_dimer[i].s1,"Sz",xy_dimer[i].s2,"S+",xy_dimer[j].s1,"S-",xy_dimer[j].s2;
                    DDmpo += 0.50,"Sz",xy_dimer[i].s1,"Sz",xy_dimer[i].s2,"S-",xy_dimer[j].s1,"S+",xy_dimer[j].s2;
                    DDmpo += 1.00,"Sz",xy_dimer[i].s1,"Sz",xy_dimer[i].s2,"Sz",xy_dimer[j].s1,"Sz",xy_dimer[j].s2;
                    DDcorr = IQMPO(DDmpo);
                    ddcorr_meas = overlap(psi,DDcorr,psi);
                    printfln("ddcorr_meas = %.8f\n", ddcorr_meas);
                    dxydxy_meas.emplace_back(ddcorr_meas);
                }
            }
        }
        std::ofstream fdxydxyout("DxyiDxyj.out",std::ios::out);
        for (std::vector<double>::const_iterator i = dxydxy_meas.begin(); i != dxydxy_meas.end(); ++i)
                fdxydxyout << *i << ' ';

        std::ofstream fdxyout("Dxyi.out",std::ios::out);
        for (std::vector<double>::const_iterator i = dxy_meas.begin(); i != dxy_meas.end(); ++i)
                fdxyout << *i << ' ';
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

        auto Xmpoi = AutoMPO(sites);
        auto Xmpoj = AutoMPO(sites);
        auto Xopi = IQMPO(Xmpoi);
        auto Xopj = IQMPO(Xmpoj);
        auto XXcorr_meas = overlap(psi,Xopi,psi);
        // measure chiral correlation
        std::vector<double> XiXj_meas={};
        std::vector<double> Xi_meas={};
        for(int i = 0; i < int(tri_plaq.size()); ++i) {
            Xmpoi = AutoMPO(sites);
            Xmpoi +=  0.5,"S+",tri_plaq[i].s1,"S-",tri_plaq[i].s2,"Sz",tri_plaq[i].s3;
            Xmpoi += -0.5,"S-",tri_plaq[i].s1,"S+",tri_plaq[i].s2,"Sz",tri_plaq[i].s3;
            Xmpoi +=  0.5,"S+",tri_plaq[i].s3,"S-",tri_plaq[i].s1,"Sz",tri_plaq[i].s2;
            Xmpoi += -0.5,"S-",tri_plaq[i].s3,"S+",tri_plaq[i].s1,"Sz",tri_plaq[i].s2;
            Xmpoi +=  0.5,"S+",tri_plaq[i].s2,"S-",tri_plaq[i].s3,"Sz",tri_plaq[i].s1;
            Xmpoi += -0.5,"S-",tri_plaq[i].s2,"S+",tri_plaq[i].s3,"Sz",tri_plaq[i].s1;
            Xopi = IQMPO(Xmpoi);
            for(int j = i; j < int(tri_plaq.size()); ++j) {
                std::vector<int> sites_tmp = { tri_plaq[i].s1, tri_plaq[i].s2, tri_plaq[i].s3, tri_plaq[j].s1, tri_plaq[j].s2, tri_plaq[j].s3 };
                //std::sort( sites_tmp.begin(), sites_tmp.end() ); // sort the pair
                //println( sites_tmp );
                for (auto n : sites_tmp ) { std::cout << n <<" "; }
                std::cout << '\n';

                // calculate correlation, Si*Sj*Sk*Sl
                Xmpoj = AutoMPO(sites);
                Xmpoj +=  0.5,"S+",tri_plaq[j].s1,"S-",tri_plaq[j].s2,"Sz",tri_plaq[j].s3;
                Xmpoj += -0.5,"S-",tri_plaq[j].s1,"S+",tri_plaq[j].s2,"Sz",tri_plaq[j].s3;
                Xmpoj +=  0.5,"S+",tri_plaq[j].s3,"S-",tri_plaq[j].s1,"Sz",tri_plaq[j].s2;
                Xmpoj += -0.5,"S-",tri_plaq[j].s3,"S+",tri_plaq[j].s1,"Sz",tri_plaq[j].s2;
                Xmpoj +=  0.5,"S+",tri_plaq[j].s2,"S-",tri_plaq[j].s3,"Sz",tri_plaq[j].s1;
                Xmpoj += -0.5,"S-",tri_plaq[j].s2,"S+",tri_plaq[j].s3,"Sz",tri_plaq[j].s1;
                Xopj = IQMPO(Xmpoj);
                XXcorr_meas = -overlap(psi,Xopi,Xopj,psi);  // Note that "-" comes of i^2
                printfln("XXcorr_meas = %.8f\n", XXcorr_meas); 
                XiXj_meas.emplace_back(XXcorr_meas);
                if(j==i) Xi_meas.emplace_back(overlap(psi,Xopi,psi));
            }
        }
        std::ofstream fXiXjout("XiXj.out",std::ios::out);
        for (std::vector<double>::const_iterator i = XiXj_meas.begin(); i != XiXj_meas.end(); ++i)
                fXiXjout << *i << ' ';

        std::ofstream fXiout("Xi.out",std::ios::out);
        for (std::vector<double>::const_iterator i = Xi_meas.begin(); i != Xi_meas.end(); ++i)
                fXiout << *i << ' ';
    }

    return 0; 
    }
