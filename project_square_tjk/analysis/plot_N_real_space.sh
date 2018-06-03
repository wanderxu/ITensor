#!/bin/bash

source cal_para.sh
firststep=7
maxstep=8

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
      kfac=1
      let nn=$Nx*$Ny

      ##  rename several files for automatically processing

      #tagarray=$( echo "S Sz Spom Dx Dy Dxy X")
      tagarray=$( echo "Ntot Nup Ndn")
      for ((istep=$firststep;istep<=$maxstep;istep++)); do
      for tag in $tagarray; do
        cd $WORKDIR
        datafile=$pretag/project_tjk/run/${maindir}/step${istep}/${tag}.out

        if [ -f $datafile ];  then
            echo " processing $maindir's ${tag} ..."

            if [ ! -d $maindir ]; then
                mkdir $maindir
            fi
            cd $maindir
            if [ ! -d step$istep ]; then
                mkdir step$istep
            fi
            cd step$istep

            # set head
cat>plot_${tag}_real_space.gnu<<endin
set terminal postscript eps enhanced color
set output "${tag}_real_space.eps"
set bmargin 5
set lmargin 12
set xtics font ',20'
set ytics font ',20'
set xlabel font 'Times-Roman, 24' offset 0,-1
set ylabel font 'Times-Roman, 24' offset -2,0
set key Left right reverse font 'Times-Roman,22' spacing 1.5
#set xlabel 'k_1'
#set ylabel 'S({/Times-Roman-Bold k})'
set size ratio -1
unset border
unset tics
plot \\
endin

            # split
            cat $datafile |sed -E -e 's/[[:blank:]]+/\n/g'> ${tag}_column.dat
            for ((i=0; i<$nn; i++)); do
                startp=$( echo "$i $Ny $Nx" |awk '{print int($1/$2)*1.0-($1%$2)*0.5, ($1%$2)*sqrt(3.0)/2.0}' )
                endp=$( echo "$i $Ny $Nx" |awk '{print int($1/$2)*1.0-($1%$2)*0.5-0.5, ($1%$2)*sqrt(3.0)/2.0+sqrt(3.0)/2.0}' )
                #dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?$1:-$1}' ${tag}_column.dat )
                dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print 3*$1}' ${tag}_column.dat ) # make it sharper
                #echo " \"<echo '$startp \n $endp'\" with l lw $dvalue lc 1 not, \\" >> plot_${tag}_real_space.gnu
                echo " \"<echo '$startp'\" with p ps $dvalue lc 1 pt 6 not, \\" >> plot_${tag}_real_space.gnu
            done
            # plot some a2 direction lines 
            for ((ix=0; ix<$Nx; ix++)); do
                startp=$( echo "$ix $Ny $Nx" |awk '{print $1, 0}' )
                endp=$( echo "$ix $Ny $Nx" |awk '{print -($2-1)*0.5+$1, ($2-1)*sqrt(3.0)/2.0}' )
                echo " \"<echo '$startp \n $endp'\" with l lc -1 lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            done
            #### plot some horizontal lines 
            ###echo "[0:   ($Nx-1.0)] 0                 lc -1  lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            ###echo "[-0.5:($Nx-1.5)] 1.0*sqrt(3.0)/2.0 lc -1  lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            ###echo "[-1.0:($Nx-2.0)] 2.0*sqrt(3.0)/2.0 lc -1  lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            ###echo "[-1.5:($Nx-2.5)] 3.0*sqrt(3.0)/2.0 lc -1  lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            ###echo "[-2.0:($Nx-3.0)] 4.0*sqrt(3.0)/2.0 lc -1  lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            ###echo "[-2.5:($Nx-3.5)] 5.0*sqrt(3.0)/2.0 lc -1  lw 0.1 not" >> plot_${tag}_real_space.gnu
            ####echo "[-3.0:($Nx-4.0)] 6.0*sqrt(3.0)/2.0 lc -1  lw 0.1 not" >> plot_${tag}_real_space.gnu
            gnuplot plot_${tag}_real_space.gnu

        fi
      done
      done
    done
  done
done
