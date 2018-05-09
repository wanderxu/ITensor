#!/bin/bash

source cal_para.sh

WORKDIR="$PWD"
echo $WORKDIR
##EXE=$HOME/mywork/itensor/project/dmrg_tri
cd $WORKDIR
for Ny in ${Nyarray}; do
  for Nx  in ${Nxarray}; do
    for J2 in ${J2array}; do
      # prepare work
      maindir=Nx${Nx}_Ny${Ny}_Jone${J1}_Jtwo${J2}_gone${gamma1}_gtwo${gamma2}
      pretag=$HOME/mycode/itensor/
      exe=analysis_SiSj_deduct_bg.py

      ##  rename several files for automatically processing
      if [ -f $pretag/project/run/${maindir}/Siz.out ]; then
          #cp -f $pretag/project/run/${maindir}/Siz.out $pretag/project/run/${maindir}/Si.out
          #cp -f $pretag/project/run/${maindir}/Siz.out $pretag/project/run/${maindir}/Szi.out
          cp -f $pretag/project/run/${maindir}/Sip.out $pretag/project/run/${maindir}/Si.out  # as <s+>=0, use it not to deduct background
          cp -f $pretag/project/run/${maindir}/Sip.out $pretag/project/run/${maindir}/Szi.out  # as <s+>=0, use it not to deduct background
          echo ""
      fi
      if [ -f $pretag/project/run/${maindir}/SiSjzz.out ]; then
          mv $pretag/project/run/${maindir}/SiSjzz.out $pretag/project/run/${maindir}/SziSzj.out
      fi
      if [ -f $pretag/project/run/${maindir}/SiSjpm.out ]; then
          mv -f $pretag/project/run/${maindir}/SiSjpm.out $pretag/project/run/${maindir}/SpomiSpomj.out
          mv -f $pretag/project/run/${maindir}/Sip.out $pretag/project/run/${maindir}/Spomi.out
      fi
      if [ -f $pretag/project/run/${maindir}/Siz.out ]; then
          ## attention here, use Spomi.out (which is zero), equals not deduct background
          cp -f $pretag/project/run/${maindir}/Spomi.out $pretag/project/run/${maindir}/Si.out
          cp -f $pretag/project/run/${maindir}/Spomi.out $pretag/project/run/${maindir}/Szi.out
      fi

      tagarray=$( echo "S Sz Spom Dx Dy Dxy X")
      for tag in $tagarray; do
        cd $WORKDIR
        datafile=$pretag/project/run/${maindir}/${tag}i${tag}j.out
        datafile2=$pretag/project/run/${maindir}/${tag}i.out

        if [ -f $datafile ];  then
            echo " processing $maindir's ${tag} ..."

            if [ ! -d $maindir ]; then
                mkdir $maindir
            fi
            cd $maindir

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
            cp -f $pretag/project/analysis_cordk/$exe .
            python $exe $datafile $datafile2 ${tag} >> ${maindir}.logs

            # plot
            maxv=$( sort -nrk3 ${tag}dbgk.dat |head -1|awk '{print $3}' )
            cp ../Nx24_Ny8_Jone1.0_Jtwo0.2_gone1.0_gtwo1.0/bz*.dat .
cat>plot${tag}k.gnu<<endin
set terminal postscript eps enhanced color
set output "${tag}k.eps"
set title "${tag}k" font 'Times-Roman, 24'
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
#set xrange[-0.5:0.5]
#set yrange[-0.5:0.5]
set xlabel '{/Times-Roman-Italic k_x}'
set ylabel '{/Times-Roman-Italic k_y}'
set cbrange[0:1]
set contour
set samples 100
set isosamples 100
set cntrparam levels discrete 0
set linetype 1 lc -1
set linetype 2 lc -1
set linetype 3 lc -1
set linetype 4 lc 2
splot "${tag}dbgk.dat" u 1:2:((\$3/$maxv)) notitle nocontour, \\
"bz.dat" u 1:2:(\$0) w l nocontour notitle, \\
#"bz2.dat" u 1:2:(\$0) w l nocontour notitle, \\
#-2*cos(x/2)-4*cos(x/4)*cos(sqrt(3.0)*y/4) nosurf notitle lw 2, \\
#-2*cos((x-2*pi)/2)-4*cos((x-2*pi)/4)*cos(sqrt(3.0)*(y-2*pi/sqrt(3))/4) nosurf notitle lw 4 lc 3 dt 3
#-2*cos((x+2*pi)/2)-4*cos((x-2*pi)/4)*cos(sqrt(3.0)*(y+2*pi/sqrt(3))/4) nosurf notitle lw 4 lc 3 dt 3
#-2*cos((x-4*pi/3)/2)-4*cos((x-4*pi/3)/4)*cos(sqrt(3.0)*(y-4*pi/sqrt(3))/4) nosurf notitle lw 4 lc 3 dt 3, \\
#-2*cos((x+4*pi/3)/2)-4*cos((x+4*pi/3)/4)*cos(sqrt(3.0)*(y+4*pi/sqrt(3))/4) nosurf notitle lw 4 lc 3 dt 3
endin
            gnuplot plot${tag}k.gnu

            awk '{if(NR>1) print $0}' ${tag}dbgij_xdirec.dat > ${tag}dbg0j_xdirec.dat
            awk '{if(NR>1) print $0}' ${tag}dbgij_ydirec.dat > ${tag}dbg0j_ydirec.dat
            awk '{if(NR>1) print $0}' ${tag}dbgij_xydirec.dat > ${tag}dbg0j_xydirec.dat
            awk '{if(NR>1) print $0}' ${tag}dbgij_21direc.dat > ${tag}dbg0j_21direc.dat
            # positive part
            awk '{if((NR>1)&&($2>0)) print $0}' ${tag}dbgij_xdirec.dat > ${tag}dbg0j_xdirec_p.dat
            awk '{if((NR>1)&&($2>0)) print $0}' ${tag}dbgij_ydirec.dat > ${tag}dbg0j_ydirec_p.dat
            awk '{if((NR>1)&&($2>0)) print $0}' ${tag}dbgij_xydirec.dat > ${tag}dbg0j_xydirec_p.dat
            awk '{if((NR>1)&&($2>0)) print $0}' ${tag}dbgij_21direc.dat > ${tag}dbg0j_21direc_p.dat
            # negative part
            awk '{if((NR>1)&&($2<=0)) print $0}' ${tag}dbgij_xdirec.dat > ${tag}dbg0j_xdirec_n.dat
            awk '{if((NR>1)&&($2<=0)) print $0}' ${tag}dbgij_ydirec.dat > ${tag}dbg0j_ydirec_n.dat
            awk '{if((NR>1)&&($2<=0)) print $0}' ${tag}dbgij_xydirec.dat > ${tag}dbg0j_xydirec_n.dat
            awk '{if((NR>1)&&($2<=0)) print $0}' ${tag}dbgij_21direc.dat > ${tag}dbg0j_21direc_n.dat
            if [ $tag == 'X' ]; then
                echo "cut Xdbgij_xdirec.dat to get Xdbgij_sametri_xdirec.dat"
                awk '{if((NR%2==1)&&(NR>1)) print $1/2, $2}' Xdbgij_xdirec.dat > Xdbg0j_sametri_xdirec.dat
                awk '{if(NR%2==0) print $1/2, $2}' Xdbgij_xdirec.dat > Xdbgij_sametri_xdirec.dat
                awk '{if((NR%2==1)&&(NR>1)&&($2>0)) print $1/2, $2}' Xdbgij_xdirec.dat > Xdbg0j_sametri_xdirec_p.dat
                awk '{if((NR%2==1)&&(NR>1)&&($2<=0)) print $1/2, $2}' Xdbgij_xdirec.dat > Xdbg0j_sametri_xdirec_n.dat

                echo "cut Xdbgij_ydirec.dat to get Xdbgij_sametri_ydirec.dat"
                awk '{if((NR%2==1)&&(NR>1)) print $1/2, $2}' Xdbgij_ydirec.dat > Xdbg0j_sametri_ydirec.dat
                awk '{if(NR%2==0) print $1/2, $2}' Xdbgij_ydirec.dat > Xdbgij_sametri_ydirec.dat
                awk '{if((NR%2==1)&&(NR>1)&&($2>0)) print $1/2, $2}' Xdbgij_ydirec.dat > Xdbg0j_sametri_ydirec_p.dat
                awk '{if((NR%2==1)&&(NR>1)&&($2<=0)) print $1/2, $2}' Xdbgij_ydirec.dat > Xdbg0j_sametri_ydirec_n.dat

                echo "cut Xdbgij_xydirec.dat to get Xdbgij_sametri_xydirec.dat"
                awk '{if((NR%2==1)&&(NR>1)) print $1/2, $2}' Xdbgij_xydirec.dat > Xdbg0j_sametri_xydirec.dat
                awk '{if(NR%2==0) print $1/2, $2}' Xdbgij_xydirec.dat > Xdbgij_sametri_xydirec.dat
                awk '{if((NR%2==1)&&(NR>1)&&($2>0)) print $1/2, $2}' Xdbgij_xydirec.dat > Xdbg0j_sametri_xydirec_p.dat
                awk '{if((NR%2==1)&&(NR>1)&&($2<=0)) print $1/2, $2}' Xdbgij_xydirec.dat > Xdbg0j_sametri_xydirec_n.dat

                echo "cut Xdbgij_21direc.dat to get Xdbgij_sametri_21direc.dat"
                awk '{if((NR%2==1)&&(NR>1)) print $1/2, $2}' Xdbgij_21direc.dat > Xdbg0j_sametri_21direc.dat
                awk '{if(NR%2==0) print $1/2, $2}' Xdbgij_21direc.dat > Xdbgij_sametri_21direc.dat
                awk '{if((NR%2==1)&&(NR>1)&&($2>0)) print $1/2, $2}' Xdbgij_21direc.dat > Xdbg0j_sametri_21direc_p.dat
                awk '{if((NR%2==1)&&(NR>1)&&($2<=0)) print $1/2, $2}' Xdbgij_21direc.dat > Xdbg0j_sametri_21direc_n.dat
            fi
        fi
      done
      if [ -f $pretag/project/analysis_cordk/${maindir}/Sdbg0j_xdirec.dat ]; then
          cd $pretag/project/analysis_cordk/${maindir}
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
