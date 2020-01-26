#!/bin/bash

source cal_para.sh
firststep=1
maxstep=10

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
      tagarray=$( echo "Ntot Nup Ndn Dxi Dyi Dxyi")
      for ((istep=$firststep;istep<=$maxstep;istep++)); do
      for tag in $tagarray; do
        cd $WORKDIR
        datafile=$pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/${tag}.out

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

            # split
            datafile=$pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/${tag}.out
            cat $datafile |sed -E -e 's/[[:blank:]]+/\n/g'> ${tag}_column.dat

            awk -v nyv=$Ny '{s+=$1}NR%nyv==0{print s/nyv;s=0}' ${tag}_column.dat > ${tag}_yave.dat
        fi
      done
      done
    done
  done
done
