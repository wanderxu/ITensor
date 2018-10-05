#!/bin/bash

source cal_para.sh

WORKDIR="$PWD"
echo $WORKDIR
##EXE=$HOME/mywork/itensor/project/dmrg_tri
exe=$HOME/mycode/itensor/project_tjk_XC/analysis-moresweeps/chi-square-extrapolation.py
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

      cd $WORKDIR

      echo " processing $maindir stepinf ..."

      if [ ! -d $maindir ]; then
          mkdir $maindir
      fi
      cd $maindir
      if [ ! -d stepinf ]; then
          mkdir stepinf
      fi
      cd stepinf
      tagarray=$( echo "dd11 dd12 dd13 dd22 dd23 dd33")
      for tag in $tagarray; do
          python -W ignore $exe ../step7/${tag}vsx.dat  ../step8/${tag}vsx.dat  ../step9/${tag}vsx.dat  > ${tag}vsx.dat
          if [ $Nx == "48" ]; then
              xtag=$(echo "x011 x012 x013 x014 x015 x016")
          elif [ $Nx == "72" ]; then
              xtag=$(echo "x017 x018 x019 x020 x021 x022 x023 x024")
          fi
          for i in $xtag ; do
              python -W ignore $exe ../step7/${tag}vsx_$i.dat  ../step8/${tag}vsx_$i.dat  ../step9/${tag}vsx_$i.dat  > ${tag}vsx_$i.dat
          done
      done
      tagarray=$( echo "Dxdbgcorr Dydbgcorr Dxydbgcorr")
      for tag in $tagarray; do
          python -W ignore $exe ../step7/${tag}_vsx.dat  ../step8/${tag}_vsx.dat  ../step9/${tag}_vsx.dat  > ${tag}_vsx.dat
      done
    done
  done
done
