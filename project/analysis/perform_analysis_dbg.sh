#!/bin/bash

source cal_para.sh

WORKDIR="$PWD"
echo $WORKDIR
##EXE=$HOME/mywork/itensor/project/dmrg_tri
cd $WORKDIR
for Ny in ${Nyarray}; do
  for Nx  in ${Nxarray}; do
    for J2 in ${J2array}; do
      # prepare work
      maindir=Nx${Nx}_Ny${Ny}_Jone${J1}_Jtwo${J2}_gone${gamma1}_gtwo${gamma2}
      pretag=$HOME/mycode/itensor/
      exe=analysis_SiSj_deduct_bg.py

      ##  rename several files for automatically processing
      if [ -f $pretag/project/run/${maindir}/Siz.out ]; then
          cp -f $pretag/project/run/${maindir}/Siz.out $pretag/project/run/${maindir}/Si.out
          cp -f $pretag/project/run/${maindir}/Siz.out $pretag/project/run/${maindir}/Szi.out
      fi
      if [ -f $pretag/project/run/${maindir}/SiSjzz.out ]; then
          mv $pretag/project/run/${maindir}/SiSjzz.out $pretag/project/run/${maindir}/SziSzj.out
      fi
      if [ -f $pretag/project/run/${maindir}/SiSjpm.out ]; then
          mv -f $pretag/project/run/${maindir}/SiSjpm.out $pretag/project/run/${maindir}/SpomiSpomj.out
          mv -f $pretag/project/run/${maindir}/Sip.out $pretag/project/run/${maindir}/Spomi.out
      fi

      tagarray=$( echo "S Sz Spom Dx Dy Dxy X")
      for tag in $tagarray; do
        cd $WORKDIR
        datafile=$pretag/project/run/${maindir}/${tag}i${tag}j.out
        datafile2=$pretag/project/run/${maindir}/${tag}i.out

        if [ -f $datafile ];  then
            echo " processing $maindir's ${tag} ..."

            if [ ! -d $maindir ]; then
                mkdir $maindir
            fi
            cd $maindir

cat>model_para.py<<endin
Nx = $Nx
Ny = $Ny
N = Nx*Ny
yperiodic = True
J1 = $J1
J2 = $J2
gamma1 = $gamma1
gamma2 = $gamma2
endin
            cp -f $pretag/project/analysis/$exe .
            python $exe $datafile $datafile2 ${tag} >> ${maindir}.logs

            # plot
            maxv=$( sort -nrk3 ${tag}dbgk.dat |head -1|awk '{print $3}' )
cat>plot${tag}k.gnu<<endin
set terminal postscript eps enhanced color
set output "${tag}k.eps"
set bmargin 5
set lmargin 12
set xtics font ',20'
set ytics font ',20'
set xlabel font 'Times-Roman, 24' offset 0,-1
set ylabel font 'Times-Roman, 24' offset -2,0
set key font 'Times-Roman,22'
set size ratio -1
set pm3d map
set pm3d interpolate 7,7
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
set cbrange[0:1]
splot "${tag}dbgk.dat" u 1:2:(\$3/$maxv) notitle
endin
            gnuplot plot${tag}k.gnu
        fi
    done

   done
  done
done
