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
exe=analysis_SiSj_deduct_bg.py
cp $pretag/project/analysis/$exe .

datafile=$pretag/project/run/${maindir}/SiSj.out
datafile2=$pretag/project/run/${maindir}/Siz.out
if [ -f $datafile ];  then
    python $exe $datafile $datafile2 ss >> ${maindir}.logs
fi

datafile=$pretag/project/run/${maindir}/DxiDxj.out
datafile2=$pretag/project/run/${maindir}/Dxi.out
if [ -f $datafile ];  then
    python $exe $datafile $datafile2 ddx >> ${maindir}.logs
fi

datafile=$pretag/project/run/${maindir}/DyiDyj.out
datafile2=$pretag/project/run/${maindir}/Dyi.out
if [ -f $datafile ];  then
    python $exe $datafile $datafile2 ddy >> ${maindir}.logs
fi

datafile=$pretag/project/run/${maindir}/DxyiDxyj.out
datafile2=$pretag/project/run/${maindir}/Dxyi.out
if [ -f $datafile ];  then
    python $exe $datafile $datafile2 ddxy >> ${maindir}.logs
fi

datafile=$pretag/project/run/${maindir}/XiXj.out
datafile2=$pretag/project/run/${maindir}/Xi.out
if [ -f $datafile ];  then
    python $exe $datafile $datafile2 xx >> ${maindir}.logs
fi

  done
done
