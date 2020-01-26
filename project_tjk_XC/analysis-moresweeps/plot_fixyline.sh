#!/bin/bash

source cal_para.sh
firststep=7
maxstep=9

WORKDIR="$PWD"
echo $WORKDIR
cd $WORKDIR
for Ny in ${Nyarray}; do
  for Nx  in ${Nxarray}; do
    for J2 in ${J2array}; do
      cd $WORKDIR
      # prepare work
      maindir=Nx${Nx}_Ny${Ny}_Jone${J1}_Jtwo${J2}_gone${gamma1}_gtwo${gamma2}_t${t1}_idop${idop}
      maintag=Nx${Nx}_Ny${Ny}_Jone${J1}_Jtwo${J2}_gone${gamma1}_gtwo${gamma2}_t${t1}_idop${idop}
      #pretag=$HOME/mycode/itensor/
      pretag=$WORKDIR
      kfac=1
      let nxa2=$Nx*$kfac+2
      let nya1=$Ny*$kfac+1
      let nyreal=$Ny*$kfac

      #tagarray=$( echo "S Sz Spom Dx Dy Dxy X")
      tagarray=$( echo "S")
      for tag in $tagarray; do
        cd $WORKDIR
        for ((istep=$firststep;istep<=$maxstep;istep++)); do
            datadir=$WORKDIR/${maindir}/step${istep}
            if [ -d $datadir ]; then
              cd $datadir

              echo " processing $datadir's ${tag} ..."

              for ((i=0; i<$nya1; i++)); do
                  awk -v iv=$i nxa2v=$nxa2'{if(NR>iv*nxa2v && NR<=(iv+1)*nxa2v) print $0}' ${tag}dbgk.dat > ${tag}fixky_line${i}.dat
                  k2tag=$( echo "$i $Ny" | awk '{print ($1-$2/2)}' )
cat>plot.tmp.gnu<<endin
set terminal postscript eps enhanced color
set output "${tag}fixky_line${i}.eps"
set bmargin 5
set lmargin 12
set xtics font ',20'
set ytics font ',20'
set xlabel font 'Times-Roman, 24' offset 0,-1
set ylabel font 'Times-Roman, 24' offset -2,0
set key Left right reverse font 'Times-Roman,22' spacing 1.5
set xrange[-0.6:0.6]
set xtics ("-1/2" -6.0/12, "" -5.0/12, "-1/3" -4.0/12, "-1/4" -3.0/12, "-1/6" -2.0/12, "" -1.0/12, "0" 0, \
             "1/2" 6.0/12, "" 5.0/12, "1/3" 4.0/12, "1/4" 3.0/12, "1/6" 2.0/12, "" 1.0/12, "0" 0)
set grid xtics
set xlabel 'k_1'
set ylabel 'S({/Times-Roman-Bold k})'
plot \\
         "${tag}fixky_line${i}.dat" u 1:3 w lp ps 1.5 lw 1.5 title "k_2=${k2tag}/${Ny}"
endin
                  gnuplot plot.tmp.gnu
                  eps2png  ${tag}fixky_line${i}.eps
              done
            fi
        done
      done
    done
  done
done
