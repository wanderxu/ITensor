#!/bin/bash

source cal_para.sh
firststep=7
maxstep=9

WORKDIR="$PWD"
echo $WORKDIR
##EXE=$HOME/mywork/itensor/project/dmrg_tri
cd $WORKDIR
for Ny in ${Nyarray}; do
  for Nx  in ${Nxarray}; do
    for J2 in ${J2array}; do
      # prepare work
      maindir=Nx${Nx}_Ny${Ny}_Jone${J1}_Jtwo${J2}_gone${gamma1}_gtwo${gamma2}_t${t1}_${t2}_idop${idop}
      pretag=$HOME/mycode/itensor/
      kfac=1
      let nn=$Nx*$Ny

      ##  rename several files for automatically processing

      #tagarray=$( echo "S Sz Spom Dx Dy Dxy X")
      tagarray=$( echo "Ntot Nup Ndn")
      for ((istep=$firststep;istep<=$maxstep;istep++)); do
      for tag in $tagarray; do
        cd $WORKDIR
        datafile=$pretag/project_square_tjk/run/${maindir}/step${istep}/${tag}.out

        if [ -f $datafile ];  then
            echo " processing $maindir step${istep}'s ${tag} ..."

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

            # plot some a2 direction lines 
            for ((ix=0; ix<$Nx; ix++)); do
                startp=$( echo "$ix $Ny $Nx" |awk '{print $1, 0}' )
                endp=$( echo "$ix $Ny $Nx" |awk '{print $1, ($2-1)}' )
                echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            done

            #### plot some a1+a2 direction lines (1.0, 1.0)
            ###for ((ix=0; ix<($Nx-$Ny+1); ix++)); do
            ###    startp=$( echo "$ix $Ny $Nx" |awk '{print $1, 0}' )
            ###    endp=$( echo "$ix $Ny $Nx" |awk '{print  ($2-1)*1.0+$1, ($2-1)*1.0}' )
            ###    echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            ###done
            ###for ((ix=($Nx-$Ny+1); ix<$Nx; ix++)); do
            ###    startp=$( echo "$ix $Ny $Nx" |awk '{print $1, 0}' )
            ###    endp=$( echo "$ix $Ny $Nx" |awk '{print  ($2-$1+$3-$2-1)*1.0+$1, ($2-$1+$3-$2-1)*1.0}' )
            ###    echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            ###done
            ###for ((ix=-1; ix>-($Nx-$Ny-1); ix--)); do
            ###    startp=$( echo "$ix $Ny $Nx" |awk '{print $1*1.0, -$1*1.0}' )
            ###    endp=$( echo "$ix $Ny $Nx" |awk '{print  ($2-1)*1.0+$1, ($2-1)*1.0}' )
            ###    echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            ###done

            # plot some horizontal lines 
            for ((ix=0; ix<$Ny; ix++)); do
                startp=$( echo "$ix $Ny $Nx" |awk '{print -$1*0.0, $1*1.0}' )
                endp=$( echo "$ix $Ny $Nx" |awk '{print -$1*0.0+$3-1,$1*1.0}' )
                echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            done

            # split
            datafile=$pretag/project_square_tjk/run/${maindir}/step${istep}/${tag}.out
            # split
            cat $datafile |sed -E -e 's/[[:blank:]]+/\n/g'> ${tag}_column.dat
            for ((i=0; i<$nn; i++)); do
                startp=$( echo "$i $Ny $Nx" |awk '{print int($1/$2)*1.0-($1%$2)*0.0, ($1%$2)*1.0}' )
                endp=$( echo "$i $Ny $Nx" |awk '{print int($1/$2)*1.0-($1%$2)*0.0-0.0, ($1%$2)*1.0+1.0}' )
                #dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?$1:-$1}' ${tag}_column.dat )
                dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print 4*$1}' ${tag}_column.dat ) # make it sharper
                echo " \"<echo '$startp'\" with p ps $dvalue lc 1 pt 6 lw 1.5 not, \\" >> plot_${tag}_real_space.gnu
            done

            gnuplot plot_${tag}_real_space.gnu

            awk '{s+=$1}NR%4==0{print s/4;s=0}' ${tag}_column.dat > ${tag}_yave.dat
        fi
      done
      done
    done
  done
done
