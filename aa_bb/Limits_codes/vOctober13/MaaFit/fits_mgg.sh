#!/bin/bash

#change windows for cuts
declare -a radion=("300" "500" "700" "1000")
declare -a winl0=("250" "450" "650" "950")
declare -a winu0=("350" "550" "750" "1050")
#
declare -a winl1=("250" "450" "650" "950")
declare -a winu1=("350" "550" "750" "1050")
# doWP_tomgg.cc
# models_mtot_exp.rs

for (( i = 4 ; i < 5 ; i++ )); do # for each working point
  mkdir radlim_CSV_WP$i
  #  for (( j = 0 ; j < ${#radion[@]} ; j++ )); do # for each mass 
  for (( j = 0 ; j < 1; j++ )); do # for each mass 
  #  for (( j = 0 ; j < 1 ; j++ )); do # for each mass 
	# the signal file
	sed -i -r -e "s/WP[0-9]/WP$i/g" doWP_tomgg.cc
	sed -i -r -e "s/m[0-9]+.root/m${radion[$j]}.root/g" doWP_tomgg.cc
	sed -i -r -e "s/m[0-9]+_8TeV_nm_m[0-9]+.root/m${radion[$j]}_8TeV_nm_m${radion[$j]}.root/g" doWP_tomgg.cc
	# the legend
	sed -i -r -e "s/[0-9]+ GeV\"\);/${radion[$j]} GeV\"\);/g" doWP_tomgg.cc
        # name of files
	sed -i -r -e "s/mH[0-9]+/mH${radion[$j]}/g" doWP_tomgg.cc
	#sed -i -r -e "s/m[0-9]+\_regression/${radion[$j]}\_regression/g" doWP_tomgg.cc
	#sed -i -r -e "s/m[0-9]+\_/m${radion[$j]}\_/g" doWP_tomgg.cc
	#sed -i -r -e "s/m[0-9]+\_/m${radion[$j]}\_/g" models_mtot_exp.rs
	# the window mtot > 550 || mtot < 450
	#sed -i -r -e "s/mtot > [0-9]+ \|\| mtot < [0-9]+/mtot > ${winu[$j]} \|\| mtot < ${winl[$j]}/g" doWP_tomgg.cc
        echo WP$i MR ${radion[$j]}
        # the cuts
        #  TString cut0 = "&& mtot > 455 && mtot < 550 "; //"&& 1>0";//
        #  TString cut1 = "&& mtot > 455 && mtot < 550 "; // "&& 1>0";//
        #  TString cutj0 = "&& mjj > 90 && mjj < 170 "; //"&& 1>0";//
        #  TString cutj1 = "&& mjj > 100 && mjj < 160 "; // "&& 1>0";//
	# the signal parameters (median to fit)
	sed -i -r -e "s/TString cut0 = \"\&\& mtot > [0-9]+ \&\& mtot < [0-9]+/TString cut0 = \"\&\& mtot > ${winl0[$j]} \&\& mtot < ${winu0[$j]}/g" doWP_tomgg.cc
	sed -i -r -e "s/TString cut1 = \"\&\& mtot > [0-9]+ \&\& mtot < [0-9]+/TString cut1 = \"\&\& mtot > ${winl1[$j]} \&\& mtot < ${winu1[$j]}/g" doWP_tomgg.cc
        #
	mkdir radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	root -l -q runfits_mgg.C >> radlim_CSV_WP$i/radlim${radion[$j]}_CSV/log_radlim${radion[$j]}.txt
	mv hgg.* radlim_CSV_WP$i/radlim${radion[$j]}_CSV
	mv .txt radlim_CSV_WP$i/radlim${radion[$j]}_CSV
	#also colect the plots
	mv databkgoversig*  radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	mv sigmodel*  radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	mv remenber.txt radlim${radion[$j]}_CSV
	## create limits root files for each mass
	cd radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	combine -M Asymptotic hgg.mH${radion[$j]}_8TeV.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.txt
        mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.root
	combine -M Asymptotic hgg.mH${radion[$j]}_8TeV_onecat.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_onecat.txt
        mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_onecat.root
	combine -M Asymptotic hgg.mH${radion[$j]}_8TeV.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_nosyst.txt
        mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.root
	cd ../..
  done # mass
done # wp
