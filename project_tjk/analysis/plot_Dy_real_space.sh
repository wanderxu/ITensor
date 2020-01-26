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
      tagarray=$( echo "Dy")
      for ((istep=$firststep;istep<=$maxstep;istep++)); do
      for tag in $tagarray; do
        cd $WORKDIR
        datafile=$pretag/project_tjk/run/${maindir}/step${istep}/${tag}i.out

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
            # plot some a2 direction lines 
            for ((ix=0; ix<$Nx; ix++)); do
                startp=$( echo "$ix $Ny $Nx" |awk '{print $1, 0}' )
                endp=$( echo "$ix $Ny $Nx" |awk '{print -($2-1)*0.5+$1, ($2-1)*sqrt(3.0)/2.0}' )
                echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            done

            # plot some a1+a2 direction lines (0.5, sqrt(3.0)/2.0
            for ((ix=0; ix<($Nx-$Ny+1); ix++)); do
                startp=$( echo "$ix $Ny $Nx" |awk '{print $1, 0}' )
                endp=$( echo "$ix $Ny $Nx" |awk '{print  ($2-1)*0.5+$1, ($2-1)*sqrt(3.0)/2.0}' )
                echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            done
            for ((ix=($Nx-$Ny+1); ix<$Nx; ix++)); do
                startp=$( echo "$ix $Ny $Nx" |awk '{print $1, 0}' )
                endp=$( echo "$ix $Ny $Nx" |awk '{print  ($2-$1+$3-$2-1)*0.5+$1, ($2-$1+$3-$2-1)*sqrt(3.0)/2.0}' )
                echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            done
            #for ((ix=-1; ix>-($Nx-$Ny-1); ix--)); do
            for ((ix=-1; ix>-($Ny-1); ix--)); do
                startp=$( echo "$ix $Ny $Nx" |awk '{print $1*0.5, -$1*sqrt(3.0)/2.0}' )
                endp=$( echo "$ix $Ny $Nx" |awk '{print  ($2-1)*0.5+$1, ($2-1)*sqrt(3.0)/2.0}' )
                echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            done

            # plot some horizontal lines 
            for ((ix=0; ix<$Ny; ix++)); do
                startp=$( echo "$ix $Ny $Nx" |awk '{print -$1*0.5, $1*sqrt(3.0)/2.0}' )
                endp=$( echo "$ix $Ny $Nx" |awk '{print -$1*0.5+$3-1,$1*sqrt(3.0)/2.0}' )
                echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_${tag}_real_space.gnu
            done

            # split
            datafile=$pretag/project_tjk/run/${maindir}/step${istep}/${tag}i.out
            cat $datafile |sed -E -e 's/[[:blank:]]+/\n/g'> ${tag}i_column.dat
            for ((i=0; i<$nn; i++)); do
                startp=$( echo "$i $Ny $Nx" |awk '{print int($1/$2)*1.0-($1%$2)*0.5-0.5*0.3, ($1%$2)*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.3}' )
                endp=$( echo "$i $Ny $Nx" |awk '{print int($1/$2)*1.0-($1%$2)*0.5-0.5*0.7, ($1%$2)*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.7}' )
                dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?20*$1:-20*$1}' ${tag}i_column.dat )
                #dvalue=$( awk -v iv=$i function abs(v) {return v < 0 ? -v : v} '{if(NR==(iv+1)) print 100*abs($1*$1)}' ${tag}i_column.dat )
                echo " \"<echo '$startp \n $endp'\" with l lw $dvalue lc 1 not, \\" >> plot_${tag}_real_space.gnu
            done

            gnuplot plot_${tag}_real_space.gnu

        fi
      done
      done
    done
  done
done
