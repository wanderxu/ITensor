#!/bin/bash

source cal_para.sh
#istep=9
istep=inf

WORKDIR="$PWD"
echo $WORKDIR
##EXE=$HOME/mywork/itensor/project/dmrg_tri
#exe=$HOME/mycode/itensor/project_tjk_XC/analysis-moresweeps/chi-square.py
#exe=$HOME/mycode/itensor/project_tjk_XC/analysis-moresweeps/chi-square-v2.py
exe=$HOME/mycode/itensor/project_tjk_XC/analysis-moresweeps/chi-square-v3.py
exelog=$HOME/mycode/itensor/project_tjk_XC/analysis-moresweeps/chi-square-exp.py
cd $WORKDIR
for Ny in ${Nyarray}; do
  for Nx  in ${Nxarray}; do
    for J2 in ${J2array}; do
      # prepare work
      maindir=Nx${Nx}_Ny${Ny}_Jone${J1}_Jtwo${J2}_gone${gamma1}_gtwo${gamma2}_t${t1}_idop${idop}
      pretag=$HOME/mycode/itensor/
      kfac=1
      let nn=$Nx*$Ny

      ##  rename several files for automatically processing

      #for ((istep=$firststep;istep<=$maxstep;istep++)); do
        cd $WORKDIR
        #datafile=$pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Ntot.out
        datafile=$pretag/project_tjk_XC/analysis-moresweeps/${maindir}/step${istep}/dd11vsx.dat

        if [ -f $datafile ];  then
            echo " processing $maindir step$istep ..."

            if [ ! -d $maindir ]; then
                mkdir $maindir
            fi
            cd $maindir
            if [ ! -d step$istep ]; then
                mkdir step$istep
            fi
            cd step$istep
            tagarray=$( echo "dd11 dd12 dd13 dd22 dd23 dd33")
            #tagarray=""
            for tag in $tagarray; do
                # a powrlaw fit
                python -W ignore $exe ${tag}vsx.dat > tmp.plot
                mv tmp.pdf ${tag}vsx_ft.pdf
                pdf2eps ${tag}vsx_ft.pdf
                eps2png ${tag}vsx_ft.eps
                sed -i 's/^/+/g' tmp.plot
                #cat tmp.plot |tr -d '\n' > ${tag}vsx_fitfunc
                func=$( cat tmp.plot |tr -d '\n' )
                rm tmp.plot

                # a exponential fit
                python -W ignore $exelog ${tag}vsx.dat > tmp.plot
                funclog=$( cat tmp.plot |tr -d '\n' )
                rm tmp.plot

                cd ..
                rm plot_${tag}_fit.gnu
                cp plot_${tag}.gnu plot_${tag}_fit.gnu
                sed -i 's/\.eps/_fit.eps/g'  plot_${tag}_fit.gnu
                echo "\"stepinf/${tag}vsx.dat\" u 1:2     w p lc 7 pt 11  ps 3 title \" inf\", \\" >> plot_${tag}_fit.gnu
                echo "\"stepinf/${tag}vsx.dat\" u 1:(-\$2) w p lc 7 pt 10  ps 3 notitle , \\" >> plot_${tag}_fit.gnu
                #echo "exp($func) w l lc 7 notitle" >> plot_${tag}_fit.gnu
                #echo "exp($func) w l lc 7" >> plot_${tag}_fit.gnu
                echo "abs($func) w l lc 7 title \"\|$func\|\", \\" >> plot_${tag}_fit.gnu
                echo "exp($funclog) w l lc -1 title \"exp\($funclog\)\"" >> plot_${tag}_fit.gnu
                gnuplot plot_${tag}_fit.gnu
                cd step$istep

                ###ytagarray=$( echo "y0 y1 y2 y3")
                for ytag in $ytagarray; do
                    # a powrlaw fit
                    python -W ignore $exe ${tag}vsx_${ytag}.dat > tmp.plot
                    mv tmp.pdf ${tag}vsx_${ytag}_ft.pdf
                    pdf2eps ${tag}vsx_${ytag}_ft.pdf
                    eps2png ${tag}vsx_${ytag}_ft.eps
                    sed -i 's/^/+/g' tmp.plot
                    #cat tmp.plot |tr -d '\n' > ${tag}vsx_${ytag}_fitfunc
                    func=$( cat tmp.plot |tr -d '\n' )
                    rm tmp.plot

                    # a exponential fit
                    python -W ignore $exelog ${tag}vsx_${ytag}.dat > tmp.plot
                    funclog=$( cat tmp.plot |tr -d '\n' )
                    rm tmp.plot

                    cd ..
                    rm plot_${tag}_${ytag}_fit.gnu
                    # construct plot_${tag}_${ytag}.gnu
                    cp plot_${tag}.gnu plot_${tag}_${ytag}.gnu
                    sed -i "s/vsx/vsx_${ytag}/g" plot_${tag}_$ytag.gnu
                    sed -i "s/${tag}\.eps/${tag}_${ytag}.eps/g" plot_${tag}_$ytag.gnu

                    cp plot_${tag}_${ytag}.gnu plot_${tag}_${ytag}_fit.gnu
                    sed -i 's/\.eps/_fit.eps/g'  plot_${tag}_${ytag}_fit.gnu
                    echo "\"stepinf/${tag}_${ytag}vsx_${ytag}.dat\" u 1:2     w p lc 7 pt 11  ps 3 title \" inf\", \\" >> plot_${tag}_${ytag}_fit.gnu
                    echo "\"stepinf/${tag}_${ytag}vsx_${ytag}.dat\" u 1:(-\$2) w p lc 7 pt 10  ps 3 notitle , \\" >> plot_${tag}_${ytag}_fit.gnu
                    #echo "exp($func) w l lc 7 notitle" >> plot_${tag}_${ytag}_fit.gnu
                    #echo "exp($func) w l lc 7" >> plot_${tag}_${ytag}_fit.gnu
                    echo "abs($func) w l lc 7 title \"\|$func\|\", \\" >> plot_${tag}_${ytag}_fit.gnu
                    echo "exp($funclog) w l lc -1 title \"exp\($funclog\)\"" >> plot_${tag}_${ytag}_fit.gnu
                    gnuplot plot_${tag}_${ytag}_fit.gnu
                    cd step$istep
                done

                ###if [ $Nx == "48" ]; then
                ###    xtag=$(echo "x011 x012 x013 x014 x015 x016")
                ###elif [ $Nx == "72" ]; then
                ###    xtag=$(echo "x017 x018 x019 x020 x021 x022 x023 x024")
                ###fi
                for i in $xtag; do
                    # a powrlaw fit
                    python -W ignore $exe ${tag}vsx_$i.dat > tmp.plot
                    mv tmp.pdf ${tag}vsx_${i}_ft.pdf
                    pdf2eps ${tag}vsx_${i}_ft.pdf
                    eps2png ${tag}vsx_${i}_ft.eps
                    sed -i 's/^/+/g' tmp.plot
                    #cat tmp.plot |tr -d '\n' > ${tag}vsx_fitfunc
                    func=$( cat tmp.plot |tr -d '\n' )
                    rm tmp.plot

                    # a exponential fit
                    python -W ignore $exelog ${tag}vsx_${i}.dat > tmp.plot
                    funclog=$( cat tmp.plot |tr -d '\n' )
                    rm tmp.plot

                    cd ..
                    rm plot_${tag}_${i}_fit.gnu
                    cp plot_${tag}_${i}.gnu plot_${tag}_${i}_fit.gnu
                    sed -i 's/\.eps/_fit.eps/g'  plot_${tag}_${i}_fit.gnu
                    echo "\"stepinf/${tag}vsx_${i}.dat\" u 1:2     w p lc 7 pt 11  ps 3 title \" inf\", \\" >> plot_${tag}_${i}_fit.gnu
                    echo "\"stepinf/${tag}vsx_${i}.dat\" u 1:(-\$2) w p lc 7 pt 10  ps 3 notitle , \\" >> plot_${tag}_${i}_fit.gnu
                    #echo "exp($func) w l lc 7 notitle" >> plot_${tag}_fit.gnu
                    #echo "exp($func) w l lc 7" >> plot_${tag}_fit.gnu
                    echo "abs($func) w l lc 7 title \"\|$func\|\", \\" >> plot_${tag}_${i}_fit.gnu
                    echo "exp($funclog) w l lc -1 title \"exp\($funclog\)\"" >> plot_${tag}_${i}_fit.gnu
                    gnuplot plot_${tag}_${i}_fit.gnu
                    cd step$istep
                done
            done
            ##tagarray=$( echo "Dxdbgcorr Dydbgcorr Dxydbgcorr Sdbgcorr")
            #tagarray=$( echo "Sdbgcorr")
            tagarray=
            for tag in $tagarray; do
                # a powrlaw fit
                python -W ignore $exe ${tag}_vsx.dat > tmp.plot
                mv tmp.pdf ${tag}_vsx_ft.pdf
                pdf2eps ${tag}_vsx_ft.pdf
                eps2png ${tag}_vsx_ft.eps
                sed -i 's/^/+/g' tmp.plot
                #cat tmp.plot |tr -d '\n' > ${tag}_vsx_fitfunc
                func=$( cat tmp.plot |tr -d '\n' )
                rm tmp.plot

                # a exponential fit
                python -W ignore $exelog ${tag}_vsx.dat > tmp.plot
                funclog=$( cat tmp.plot |tr -d '\n' )
                rm tmp.plot

                cd ..
                rm plot_${tag}_vsx_fit.gnu
                cp plot_${tag}_vsx.gnu plot_${tag}_vsx_fit.gnu
                sed -i 's/\.eps/_fit.eps/g'  plot_${tag}_vsx_fit.gnu
                echo "\"stepinf/${tag}_vsx.dat\" u 1:2     w p lc 7 pt 11  ps 3 title \" inf\", \\" >> plot_${tag}_vsx_fit.gnu
                echo "\"stepinf/${tag}_vsx.dat\" u 1:(-\$2) w p lc 7 pt 10  ps 3 notitle , \\" >> plot_${tag}_vsx_fit.gnu
                #echo "exp($func) w l lc 7 notitle" >> plot_${tag}_vsx_fit.gnu
                #echo "exp($func) w l lc 7" >> plot_${tag}_vsx_fit.gnu
                echo "abs($func) w l lc 7 title \"\|$func\|\", \\" >> plot_${tag}_vsx_fit.gnu
                echo "exp($funclog) w l lc -1 title \"exp\($funclog\)\"" >> plot_${tag}_vsx_fit.gnu
                gnuplot plot_${tag}_vsx_fit.gnu
                cd step$istep
                ###if [ $Nx == "48" ]; then
                ###    xtag=$(echo "x011 x012 x013 x014 x015 x016")
                ###elif [ $Nx == "72" ]; then
                ###    xtag=$(echo "x017 x018 x019 x020 x021 x022 x023 x024")
                ###    if [ $tag == "Sdbgcorr" ]; then
                ###        #xtag=$(echo "x017 x018 x019 x020 x021 x022 x023 x024 x025 x026 x027 x028 x029 x030 x031 x032 x033 x034")
                ###        xtag=$(echo "x025 x026 x027 x028 x029 x030 x031 x032 x033 x034")
                ###    fi
                ###fi
                for i in $xtag; do
                    # a powrlaw fit
                    python -W ignore $exe ${tag}_$i.dat > tmp.plot
                    mv tmp.pdf ${tag}_${i}_ft.pdf
                    pdf2eps ${tag}_${i}_ft.pdf
                    eps2png ${tag}_${i}_ft.eps
                    sed -i 's/^/+/g' tmp.plot
                    #cat tmp.plot |tr -d '\n' > ${tag}_fitfunc
                    func=$( cat tmp.plot |tr -d '\n' )
                    rm tmp.plot

                    # a exponential fit
                    python -W ignore $exelog ${tag}_$i.dat > tmp.plot
                    funclog=$( cat tmp.plot |tr -d '\n' )
                    rm tmp.plot

                    cd ..
                    rm plot_${tag}_${i}_fit.gnu
                    cp plot_${tag}_${i}.gnu plot_${tag}_${i}_fit.gnu
                    sed -i 's/\.eps/_fit.eps/g'  plot_${tag}_${i}_fit.gnu
                    echo "\"stepinf/${tag}_${i}.dat\" u 1:2     w p lc 7 pt 11  ps 3 title \" inf\", \\" >> plot_${tag}_${i}_fit.gnu
                    echo "\"stepinf/${tag}_${i}.dat\" u 1:(-\$2) w p lc 7 pt 10  ps 3 notitle , \\" >> plot_${tag}_${i}_fit.gnu
                    #echo "exp($func) w l lc 7 notitle" >> plot_${tag}_fit.gnu
                    #echo "exp($func) w l lc 7" >> plot_${tag}_fit.gnu
                    echo "abs($func) w l lc 7 title \"\|$func\|\", \\" >> plot_${tag}_${i}_fit.gnu
                    echo "exp($funclog) w l lc -1 title \"exp\($funclog\)\"" >> plot_${tag}_${i}_fit.gnu
                    gnuplot plot_${tag}_${i}_fit.gnu
                    cd step$istep
                done
            done

        fi
      #done
    done
  done
done