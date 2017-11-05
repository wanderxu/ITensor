//
// Writen by Xiao Yan Xu (wanderxu@gmail.com)
//
#ifndef __MEASURE_H_
#define __MEASURE_H_

//#include "latticeplaque.h"

namespace itensor {

typedef std::pair<std::string, int> opair;

bool cmp_by_value(const opair& lhs, const opair& rhs) {  
      return lhs.second < rhs.second;  
}

template <class Tensor>
Real
mfourbody(MPSt<Tensor>& psi, 
     SpinHalf const& sites,
     std::vector<int> const& sites_tmp,
     std::string const& op1_label,
     std::string const& op2_label,
     std::string const& op3_label,
     std::string const& op4_label )
    {
    std::vector<opair> opstr = { std::make_pair(op1_label,sites_tmp[0]),
                                 std::make_pair(op2_label,sites_tmp[1]),
                                 std::make_pair(op3_label,sites_tmp[2]),
                                 std::make_pair(op4_label,sites_tmp[3])};
    // sort by site index
    std::sort(opstr.begin(), opstr.end(), cmp_by_value);

    auto opi = sites.op(opstr[0].first,opstr[0].second);
    auto opj = sites.op(opstr[1].first,opstr[1].second);
    auto opk = sites.op(opstr[2].first,opstr[2].second);
    auto opl = sites.op(opstr[3].first,opstr[3].second);

    psi.position(opstr[0].second);
    IQTensor SStmp=psi.A(opstr[0].second);
    if( opstr[1].second != opstr[0].second ) {
        //println("i!=j");
        SStmp *= opi;
        auto ir1 = commonIndex(psi.A(opstr[0].second), psi.A(opstr[0].second+1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[0].second),Site),ir1));
        // propogate until opj
        for ( int i1 = opstr[0].second+1; i1<opstr[1].second; ++i1) { 
            SStmp *= psi.A(i1);
            SStmp *= dag(prime(psi.A(i1),Link));
        }
        if( opstr[2].second != opstr[1].second ){
            //println("j!=k");
            SStmp *= psi.A( opstr[1].second );
            SStmp *= opj;
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),Link));
            // propogate until opk
            for ( int i2 = opstr[1].second+1; i2<opstr[2].second; ++i2) { 
                SStmp *= psi.A(i2);
                SStmp *= dag(prime(psi.A(i2),Link));
            }
            if( opstr[3].second != opstr[2].second ){
                //println("k!=l");
                SStmp *= psi.A( opstr[2].second );
                SStmp *= opk;
                SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),Link));
                // propogate until opl
                for ( int i3 = opstr[2].second+1; i3<opstr[3].second; ++i3) { 
                    SStmp *= psi.A(i3);
                    SStmp *= dag(prime(psi.A(i3),Link));
                }
                SStmp *= psi.A( opstr[3].second );
                SStmp *= opl;
                auto ir3 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second-1),Link);
                SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir3)); // return i!=j!=k!=l
            }
            else { // opstr[3].second == opstr[2].second
                //println("k=l");
                SStmp *= psi.A( opstr[2].second );
                SStmp *= opk;
                //auto SStmp_tmp = SStmp*opl;
                //SStmp = prime(SStmp_tmp, Site);
                SStmp = noprime(SStmp,Site)*opl;
                auto ir3 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second-1),Link);
                SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir3)); // return i!=j!=k=l
            }
        }
        else { // opstr[2].second == opstr[1].second
            //println("j=k");
            if( opstr[3].second != opstr[2].second ){  // opstr[3].second != opstr[2].second == opstr[1].second
                //println("k!=l");
                SStmp *= psi.A( opstr[1].second );
                SStmp *= opj;
                //auto SStmp_tmp = SStmp*opk;
                //SStmp = prime(SStmp_tmp, Site);
                SStmp = noprime(SStmp,Site)*opk;
                SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),Link));
                // propogate until opl
                for ( int i3 = opstr[2].second+1; i3<opstr[3].second; ++i3) { 
                    SStmp *= psi.A(i3);
                    SStmp *= dag(prime(psi.A(i3),Link));
                }
                SStmp *= psi.A( opstr[3].second );
                SStmp *= opl;
                auto ir3 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second-1),Link);
                SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir3)); // return i!=j=k!=l
            }
            else { // opstr[3].second == opstr[2].second == opstr[1].second
                //println("k=l");
                SStmp *= psi.A( opstr[1].second );
                SStmp *= opj;
                //SStmp *= opk;
                //SStmp *= opl;
                SStmp = noprime(SStmp,Site)*opk;
                SStmp = noprime(SStmp,Site)*opl;
                auto ir3 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second-1),Link);
                SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir3));  // return i!=j=k=l
            }
        }
    }
    else{ // opstr[1].second == opstr[0].second
        //println("i=j");
        if( opstr[2].second != opstr[1].second ){ // opstr[2].second != opstr[1].second == sites[0]
            //println("j!=k");
            SStmp *= opi;
            //auto SStmp_tmp = SStmp*opj;
            //SStmp = prime(SStmp_tmp,Site);
            SStmp = noprime(SStmp,Site)*opj;
            auto ir1 = commonIndex(psi.A(opstr[1].second), psi.A(opstr[1].second+1),Link);
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),ir1));
            // propogate until opk
            for ( int i2 = opstr[1].second+1; i2<opstr[2].second; ++i2) { 
                SStmp *= psi.A(i2);
                SStmp *= dag(prime(psi.A(i2),Link));
            }
            if( opstr[3].second != opstr[2].second ){
                //println("k!=l");
                SStmp *= psi.A( opstr[2].second );
                SStmp *= opk;
                SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),Link));
                // propogate until opl
                for ( int i3 = opstr[2].second+1; i3<opstr[3].second; ++i3) { 
                    SStmp *= psi.A(i3);
                    SStmp *= dag(prime(psi.A(i3),Link));
                }
                SStmp *= psi.A( opstr[3].second );
                SStmp *= opl;
                auto ir3 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second-1),Link);
                SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir3)); // return i=j!=k!=l
            }
            else { // opstr[3].second == opstr[2].second
                //println("k=l");
                SStmp *= psi.A( opstr[2].second );
                SStmp *= opk;
                //auto SStmp_tmp = SStmp*opl;
                //SStmp = prime(SStmp_tmp, Site);
                SStmp = noprime(SStmp,Site)*opl;
                auto ir3 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second-1),Link);
                SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir3)); // return i=j!=k=l
            }
        }
        else{ // opstr[2].second == opstr[1].second == opstr[0].second
            //println("j=k");
            if( opstr[3].second != opstr[2].second ){ // opstr[3].second != opstr[2].second == opstr[1].second == opstr[0].second
                //println("k!=l");
                SStmp *= opi;
                //SStmp *= opj;
                //SStmp *= opk;
                SStmp = noprime(SStmp,Site)*opj;
                SStmp = noprime(SStmp,Site)*opk;
                auto ir2 = commonIndex(psi.A(opstr[2].second), psi.A(opstr[2].second+1),Link);
                SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),ir2));
                // propogate until opl
                for ( int i3 = opstr[2].second+1; i3<opstr[3].second; ++i3) { 
                    SStmp *= psi.A(i3);
                    SStmp *= dag(prime(psi.A(i3),Link));
                }
                SStmp *= psi.A( opstr[3].second );
                SStmp *= opl;
                auto ir3 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second-1),Link);
                SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir3)); // return i=j=k!=l
            }
            else{ // opstr[3].second = opstr[2].second == opstr[1].second == opstr[0].second
                //println("k=l");
                SStmp *= opi;
                //SStmp *= opj;
                //SStmp *= opk;
                //auto SStmp_tmp = SStmp*opl;
                //SStmp = prime(SStmp_tmp, Site);
                SStmp = noprime(SStmp,Site)*opj;
                SStmp = noprime(SStmp,Site)*opk;
                SStmp = noprime(SStmp,Site)*opl;
                //auto ir3 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second-1),Link);
                //SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir3)); // return i=j=k=l
                SStmp *= dag(prime(psi.A(opstr[3].second),Site)); // return i=j=k=l
            }
        }
    }
    return SStmp.real();
    }
} //namespace itensor

#endif
