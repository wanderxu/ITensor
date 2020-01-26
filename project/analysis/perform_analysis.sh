#!/bin/bash

source cal_para.sh

WORKDIR="$PWD"
echo $WORKDIR
##EXE=$HOME/mywork/itensor/project/dmrg_tri
cd $WORKDIR
for Ny in ${Nyarray}; do
  for Nx  in ${Nxarray}; do
    #for h in ${hxarray}; do
        cd $WORKDIR

        maindir=Nx${Nx}_Ny${Ny}_Jone${J1}_Jtwo${J2}_gone${gamma1}_gtwo${gamma2}
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

pretag=$HOME/mycode/itensor/
cp $pretag/project/analysis/analysis_SiSj.py .

datafile=$pretag/project/run/${maindir}/SiSj.out
if [ -f $datafile ];  then
    python analysis_SiSj.py $datafile ss >> ${maindir}.logs
fi

datafile=$pretag/project/run/${maindir}/DxiDxj.out
if [ -f $datafile ];  then
    python analysis_SiSj.py $datafile ddx >> ${maindir}.logs
fi

datafile=$pretag/project/run/${maindir}/DyiDyj.out
if [ -f $datafile ];  then
    python analysis_SiSj.py $datafile ddy >> ${maindir}.logs
fi

datafile=$pretag/project/run/${maindir}/DxyiDxyj.out
if [ -f $datafile ];  then
    python analysis_SiSj.py $datafile ddxy >> ${maindir}.logs
fi

datafile=$pretag/project/run/${maindir}/XiXj.out
if [ -f $datafile ];  then
    python analysis_SiSj.py $datafile xx >> ${maindir}.logs
fi

  done
done
