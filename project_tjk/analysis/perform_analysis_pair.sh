#!/bin/bash

source cal_para.sh
firststep=1
maxstep=10
taglist=$( echo "s d+id d-id dxy dx2-y2 f p+ip p-ip px py")

WORKDIR="$PWD"
echo $WORKDIR
##EXE=$HOME/mywork/itensor/project/dmrg_tri
cd $WORKDIR
for Ny in ${Nyarray}; do
  for Nx  in ${Nxarray}; do
    for J2 in ${J2array}; do
      # prepare work
      maindir=Nx${Nx}_Ny${Ny}_Jone${J1}_Jtwo${J2}_gone${gamma1}_gtwo${gamma2}_t${t1}_idop${idop}
      pretag=$HOME/mycode/itensor/
      exe=analysis_paircorr.py
      exe2=analysis_paircorr_xy.py

      for ((istep=$firststep;istep<=$maxstep;istep++)); do
        cd $WORKDIR
        #datafile=$pretag/project_tjk/run/${maindir}/step${istep}/paircorr.out
        datafile=$pretag/project_tjk/run/${maindir}/step${istep}/paircorrfixi0.out
        datafile2=$pretag/project_tjk/run/${maindir}/step${istep}/pair_bubble.out

        if [ -f $datafile ];  then
            echo " processing $maindir's step${istep} ..."

            if [ ! -d $maindir ]; then
                mkdir $maindir
            fi
            cd $maindir
            if [ ! -d step$istep ]; then
                mkdir step$istep
            fi
            cd step$istep

cat>model_para.py<<endin
Nx = $Nx
Ny = $Ny
N = Nx*Ny
yperiodic = $yperiodic
J1 = $J1
J2 = $J2
gamma1 = $gamma1
gamma2 = $gamma2
endin
            cp -f $pretag/project_tjk/analysis/$exe .
            cp -f $pretag/project_tjk/analysis/$exe2 .
            # transform the data to one column 
            sed -E -e 's/[[:blank:]]+/\n/g' $datafile >tmp.dat
            sed -E -e 's/[[:blank:]]+/\n/g' $datafile2 >tmp2.dat
            #python $exe tmp.dat tmp2.dat pair >> ${maindir}.logs
            #python $exe tmp.dat pair >> ${maindir}.logs
            python $exe tmp2.dat pair >> ${maindir}.logs
            python $exe2 tmp.dat pair >> ${maindir}.logs

            for pairtag in ${taglist}; do
              # plot
              maxv=$( sort -nrk3 pair${pairtag}k.dat |head -1|awk '{print $3}' )
              minv=$( sort -nrk3 pair${pairtag}k.dat |tail -1|awk '{print $3}' )
              echo "maxv in $pairtag = $maxv, minv in $pairtag = $minv"
cat>plot${pairtag}k.gnu<<endin
set terminal postscript eps enhanced color
set output "pair${pairtag}k.eps"
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
#set cbrange[0:1]
#splot "pair${pairtag}k.dat" u 1:2:(\$3/$maxv) notitle
splot "pair${pairtag}k.dat" u 1:2:(\$3) notitle
endin
              gnuplot plot${pairtag}k.gnu
            done
        fi
      done
   done
  done
done
