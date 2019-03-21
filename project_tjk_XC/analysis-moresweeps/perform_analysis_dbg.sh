#!/bin/bash

source cal_para.sh
firststep=7
maxstep=10

WORKDIR="$PWD"
echo $WORKDIR
##EXE=$HOME/mywork/itensor/project_tjk_XC/dmrg_tri
cd $WORKDIR
for Ny in ${Nyarray}; do
  for Nx  in ${Nxarray}; do
    for J2 in ${J2array}; do
      # prepare work
      maindir=Nx${Nx}_Ny${Ny}_Jone${J1}_Jtwo${J2}_gone${gamma1}_gtwo${gamma2}_t${t1}_idop${idop}
      pretag=$HOME/mycode/itensor/
      exe=analysis_SiSj_deduct_bg.py
      exe2=analysis_SiSj_varyx0.py # no dbg case

      for ((istep=$firststep;istep<=$maxstep;istep++)); do
      ##  rename several files for automatically processing
      if [ -f $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/SiSjzz.out ]; then
          mv $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/SiSjzz.out $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/SziSzj.out
      fi
      if [ -f $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/SiSjpm.out ]; then
          mv -f $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/SiSjpm.out $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/SpomiSpomj.out
      fi
      if [ -f $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Siz.out ]; then
          ## constuct files with zeros to not deduct background for spin correlation case
          awk '{for (i = 1; i <= NF; i++) printf("%2d", 0 )}' $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Siz.out \
          > $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Si.out
          awk '{for (i = 1; i <= NF; i++) printf("%2d", 0 )}' $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Siz.out \
          > $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Szi.out
          awk '{for (i = 1; i <= NF; i++) printf("%2d", 0 )}' $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Siz.out \
          > $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Spomi.out
      fi
      if [ -f $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Ntot.out ]; then
          cp $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/Ntot.out $pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/ni.out
      fi

      #tagarray=$( echo "S Sz Spom Dx Dy Dxy X")
      #tagarray=$( echo "S Sz Spom Dx Dy Dxy")
      tagarray=$( echo "n")
      #tagarray=$( echo "S Sz Spom")
      for tag in $tagarray; do
        cd $WORKDIR
        datafile=$pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/${tag}i${tag}j.out
        datafile2=$pretag/project_tjk_XC/run-moresweeps/${maindir}/step${istep}/${tag}i.out

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

cat>model_para.py<<endin
Nx = $Nx
Ny = $Ny
N = Nx*Ny
yperiodic = "$yperiodic"
J1 = $J1
J2 = $J2
gamma1 = $gamma1
gamma2 = $gamma2
endin
            cp -f $pretag/project_tjk_XC/analysis-moresweeps/$exe .
            python $exe $datafile $datafile2 ${tag} >> ${maindir}.logs

            cp -f $pretag/project_tjk_XC/analysis-moresweeps/$exe2 .
            python $exe2 $datafile $datafile2 ${tag} >> ${maindir}.logs

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
set xtics ("-0.5" -0.5,"" -0.375,"-0.25" -0.25,"" -0.125,"0" 0,"" 0.125, "0.25" 0.25, "" 0.375, "0.5" 0.5)
set ytics ("-0.5" -0.5,"" -0.375,"-0.25" -0.25,"" -0.125,"0" 0,"" 0.125, "0.25" 0.25, "" 0.375, "0.5" 0.5)
set grid
set cbrange[0:1]
splot "${tag}dbgk.dat" u 1:2:(\$3/$maxv) notitle
endin
            gnuplot plot${tag}k.gnu

            awk '{if(NR>1) print $0}' ${tag}dbgij_xdirec.dat > ${tag}dbg0j_xdirec.dat
            awk '{if(NR>1) print $0}' ${tag}dbgij_ydirec.dat > ${tag}dbg0j_ydirec.dat
            if [ $tag == 'X' ]; then
                echo "cut Xdbgij_xdirec.dat to get Xdbgij_sametri_xdirec.dat"
                awk '{if((NR%2==1)&&(NR>1)) print $1/2, $2}' Xdbgij_xdirec.dat > Xdbg0j_sametri_xdirec.dat

                echo "cut Xdbgij_ydirec.dat to get Xdbgij_sametri_ydirec.dat"
                awk '{if((NR%2==1)&&(NR>1)) print $1/2, $2}' Xdbgij_ydirec.dat > Xdbg0j_sametri_ydirec.dat
            fi
        fi
      done
      if [ -f $pretag/project_tjk_XC/analysis-moresweeps/${maindir}/Sdbg0j_xdirec.dat ]; then
          cd $pretag/project_tjk_XC/analysis-moresweeps/${maindir}
cat>plotloglog.gnu<<endin
set terminal postscript eps enhanced color
set output "loglog.eps"
set bmargin 5
set lmargin 12
set xtics font ',20'
set ytics font ',20'
set xlabel font 'Times-Roman, 24' offset 0,-1
set ylabel font 'Times-Roman, 24' offset -2,0
set key Left right reverse font 'Times-Roman,22' spacing 1.5
set logscale x
set logscale y
set format y "10^{%L}"
set xlabel '|i-j|'
plot \\
         "Sdbg0j_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_S(i,j)", \\
        "Dxdbg0j_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_{D_x}(i,j)", \\
        "Dydbg0j_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_{D_y}(i,j)", \\
       "Dxydbg0j_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_{D_{xy}}(i,j)", \\
 "Xdbg0j_sametri_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_X(i,j)", \\
 1.0/x**4 lc -1 dt 2 lw 1.5 title "r^{-4}
endin
     gnuplot plotloglog.gnu

cat>plotlogy.gnu<<endin
set terminal postscript eps enhanced color
set output "logC.eps"
set bmargin 5
set lmargin 12
set xtics font ',20'
set ytics font ',20'
set xlabel font 'Times-Roman, 24' offset 0,-1
set ylabel font 'Times-Roman, 24' offset -2,0
set key Left right reverse font 'Times-Roman,22' spacing 1.5
set logscale y
set format y "10^{%L}"
set xlabel '|i-j|'
plot \\
         "Sdbg0j_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_S(i,j)", \\
        "Dxdbg0j_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_{D_x}(i,j)", \\
        "Dydbg0j_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_{D_y}(i,j)", \\
       "Dxydbg0j_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_{D_{xy}}(i,j)", \\
 "Xdbg0j_sametri_xdirec.dat" u 1:(abs(\$2)) w lp ps 2 lw 1.5 title "C_X(i,j)", \\
 1.0/x**4 lc -1 dt 2 lw 1.5 title "r^{-4}
endin
     gnuplot plotlogy.gnu
     fi
   done
   done
  done
done
