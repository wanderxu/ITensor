#!/bin/bash

source cal_para.sh
firststep=9
maxstep=9

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

      for ((istep=$firststep;istep<=$maxstep;istep++)); do
        cd $WORKDIR
        datafile=$pretag/project_tjk_YC/run-moresweeps/${maindir}/step${istep}/Ntot.out

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

            # set head
cat>plot_all_real_space.gnu<<endin
set terminal postscript eps enhanced color
set output "all_real_space.eps"
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
        # odd part
        for ((ix=0; ix<$Nx; ix+=2)); do
          for ((iy=0; iy<$Ny; iy++)); do 
            startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2-0.5}' )
            endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2+0.5}' )
            echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_all_real_space.gnu
          done
        done
        # even part
        for ((ix=1; ix<$Nx; ix+=2)); do
          for ((iy=0; iy<$Ny; iy++)); do 
            startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2}' )
            endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2+1.0}' )
            echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_all_real_space.gnu
          done
        done

        # plot some a1+a2 direction lines 
        # odd part
        for ((ix=0; ix<$Nx-1; ix+=2)); do
          for ((iy=0; iy<$Ny; iy++)); do 
            startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2-0.5}' )
            endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0, $2}' )
            echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_all_real_space.gnu
          done
        done
        # even part
        for ((ix=1; ix<$Nx-1; ix+=2)); do
          for ((iy=0; iy<$Ny; iy++)); do 
            startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2}' )
            endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0, $2+0.5}' )
            echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_all_real_space.gnu
          done
        done

        # plot some a1 direction lines 
        # odd part
        for ((ix=0; ix<$Nx-1; ix+=2)); do
          for ((iy=0; iy<$Ny; iy++)); do 
            startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2-0.5}' )
            endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0, $2-1.0}' )
            echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_all_real_space.gnu
          done
        done
        # even part
        for ((ix=1; ix<$Nx-1; ix+=2)); do
          for ((iy=0; iy<$Ny; iy++)); do 
            startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2}' )
            endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0, $2-0.5}' )
            echo " \"<echo '$startp \n $endp'\" with l lc rgb 'grey' lw 0.1 not, \\" >> plot_all_real_space.gnu
          done
        done

        # density up
        tagarray=$( echo "Nup")
        for tag in $tagarray; do
            datafile=$pretag/project_tjk_YC/run-moresweeps/${maindir}/step${istep}/${tag}.out
            # split
            cat $datafile |sed -E -e 's/[[:blank:]]+/\n/g'> ${tag}_column.dat
            for ((i=0; i<$nn; i++)); do
                startp=$( echo "$i $Ny $Nx" |awk '{print int($1/$2)*sqrt(3.0)/2.0, ($1%$2)*1.0-(int($1/$2+1)%2)*0.5}' )
                #dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?$1:-$1}' ${tag}_column.dat )
                dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print 4*$1}' ${tag}_column.dat ) # make it sharper
                echo " \"<echo '$startp'\" with p ps $dvalue lc 7 pt 6 lw 1.5 not, \\" >> plot_all_real_space.gnu
            done
        done

        # density dn
        tagarray=$( echo "Ndn")
        for tag in $tagarray; do
            datafile=$pretag/project_tjk_YC/run-moresweeps/${maindir}/step${istep}/${tag}.out
            # split
            cat $datafile |sed -E -e 's/[[:blank:]]+/\n/g'> ${tag}_column.dat
            for ((i=0; i<$nn; i++)); do
                startp=$( echo "$i $Ny $Nx" |awk '{print int($1/$2)*sqrt(3.0)/2.0, ($1%$2)*1.0-(int($1/$2+1)%2)*0.5}' )
                #dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?$1:-$1}' ${tag}_column.dat )
                dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print 4*$1}' ${tag}_column.dat ) # make it sharper
                #echo " \"<echo '$startp \n $endp'\" with l lw $dvalue lc 2 not, \\" >> plot_all_real_space.gnu
                echo " \"<echo '$startp'\" with p ps $dvalue lc -1 pt 2 lw 1.5 not, \\" >> plot_all_real_space.gnu
            done
        done

        # Dy
        tagarray=$( echo "Dy")
        for tag in $tagarray; do
            datafile=$pretag/project_tjk_YC/run-moresweeps/${maindir}/step${istep}/${tag}i.out
            cat $datafile |sed -E -e 's/[[:blank:]]+/\n/g'> ${tag}i_column.dat
            let i=0
            for ((ix=0; ix<$Nx; ix++)); do
              for ((iy=0; iy<$Ny; iy++)); do 
                if ((ix%2==0)); then
                    startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2-0.5+1.0*0.3}' )
                    endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2-0.5+1.0*0.7}' )
                    dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?20*$1:-20*$1}' ${tag}i_column.dat )
                    #dvalue=$( awk -v iv=$i function abs(v) {return v < 0 ? -v : v} '{if(NR==(iv+1)) print 100*abs($1*$1)}' ${tag}i_column.dat )
                    echo " \"<echo '$startp \n $endp'\" with l lw $dvalue lc 1 not, \\" >> plot_all_real_space.gnu
                    let i+=1
                else
                    startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2+1.0*0.3}' )
                    endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0, $2+1.0*0.7}' )
                    dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?20*$1:-20*$1}' ${tag}i_column.dat )
                    #dvalue=$( awk -v iv=$i function abs(v) {return v < 0 ? -v : v} '{if(NR==(iv+1)) print 100*abs($1*$1)}' ${tag}i_column.dat )
                    echo " \"<echo '$startp \n $endp'\" with l lw $dvalue lc 1 not, \\" >> plot_all_real_space.gnu
                    let i+=1
                fi
              done
            done
        done

        # Dxy
        tagarray=$( echo "Dxy")
        for tag in $tagarray; do
            datafile=$pretag/project_tjk_YC/run-moresweeps/${maindir}/step${istep}/${tag}i.out
            cat $datafile |sed -E -e 's/[[:blank:]]+/\n/g'> ${tag}i_column.dat
            let i=0
            for ((ix=0; ix<$Nx-1; ix++)); do
              for ((iy=0; iy<$Ny; iy++)); do 
                if ((ix%2==0)); then
                    startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.3, $2-0.5+0.5*0.3}' )
                    endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.7, $2-0.5+0.5*0.7}' )
                    dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?20*$1:-20*$1}' ${tag}i_column.dat )
                    #dvalue=$( awk -v iv=$i function abs(v) {return v < 0 ? -v : v} '{if(NR==(iv+1)) print 100*abs($1*$1)}' ${tag}i_column.dat )
                    echo " \"<echo '$startp \n $endp'\" with l lw $dvalue lc 1 not, \\" >> plot_all_real_space.gnu
                    let i+=1
                else
                    startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.3, $2+0.5*0.3}' )
                    endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.7, $2+0.5*0.7}' )
                    dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?20*$1:-20*$1}' ${tag}i_column.dat )
                    #dvalue=$( awk -v iv=$i function abs(v) {return v < 0 ? -v : v} '{if(NR==(iv+1)) print 100*abs($1*$1)}' ${tag}i_column.dat )
                    echo " \"<echo '$startp \n $endp'\" with l lw $dvalue lc 1 not, \\" >> plot_all_real_space.gnu
                    let i+=1
                fi
              done
            done
        done

        # Dx
        tagarray=$( echo "Dx")
        for tag in $tagarray; do
            datafile=$pretag/project_tjk_YC/run-moresweeps/${maindir}/step${istep}/${tag}i.out
            cat $datafile |sed -E -e 's/[[:blank:]]+/\n/g'> ${tag}i_column.dat
            let i=0
            for ((ix=0; ix<$Nx-1; ix++)); do
              for ((iy=0; iy<$Ny; iy++)); do 
                if ((ix%2==0)); then
                    startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.3, $2-0.5-0.5*0.3}' )
                    endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.7, $2-0.5-0.5*0.7}' )
                    dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?20*$1:-20*$1}' ${tag}i_column.dat )
                    #dvalue=$( awk -v iv=$i function abs(v) {return v < 0 ? -v : v} '{if(NR==(iv+1)) print 100*abs($1*$1)}' ${tag}i_column.dat )
                    echo " \"<echo '$startp \n $endp'\" with l lw $dvalue lc 1 not, \\" >> plot_all_real_space.gnu
                    let i+=1
                else
                    startp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.3, $2-0.5*0.3}' )
                    endp=$( echo "$ix $iy $Ny $Nx" |awk '{print $1*sqrt(3.0)/2.0+sqrt(3.0)/2.0*0.7, $2-0.5*0.7}' )
                    dvalue=$( awk -v iv=$i '{if(NR==(iv+1)) print ($1>0)?20*$1:-20*$1}' ${tag}i_column.dat )
                    #dvalue=$( awk -v iv=$i function abs(v) {return v < 0 ? -v : v} '{if(NR==(iv+1)) print 100*abs($1*$1)}' ${tag}i_column.dat )
                    echo " \"<echo '$startp \n $endp'\" with l lw $dvalue lc 1 not, \\" >> plot_all_real_space.gnu
                    let i+=1
                fi
              done
            done
        done

        gnuplot plot_all_real_space.gnu

        fi
      done
    done
  done
done
