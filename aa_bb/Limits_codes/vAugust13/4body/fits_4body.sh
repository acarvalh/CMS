#!/bin/bash


#change for signal
#declare -a sigfiles=( "radion_300" "radion_500" "radion_700" "radion_1000") 
# n files of 500 events
#declare -a MR=(300 500 700) 
#change for radion mass cut
#sed -i -r -e "s/RadMass > [0-9]+\-50/RadMass \> ${MR[$j]}-50/g" ggAnalysis_full_fit_tree.C
#sed -i -r -e "s/RadMass < [0-9]+\+50/RadMass \< ${MR[$j]}+50/g" ggAnalysis_full_fit_tree.C
#change for data file
#declare -a radion=("300" "500" "700" "1000")
#declare -a winl=("270" "470" "670" "970")
#declare -a winu=("330" "530" "730" "1030")

declare -a radion=("500" "550" "600" "650" "700" "900" "1000" "1100" "1200")
declare -a winl=("450" "500" "550" "600" "650" "850" "900" "1000" "1100")
declare -a winu=("550" "600" "650" "700" "750" "950" "1100" "1200" "1300")

# R2GGBBFitter_mtot_side.cc
# models_mtot_exp.rs

for (( i = 4 ; i < 5 ; i++ )); do # for each working point
  mkdir radlim_CSV_WP$i

#  for (( j = 0 ; j < ${#radion[@]} ; j++ )); do # for each mass 
  for (( j = 0 ; j < 1 ; j++ )); do # for each mass 
#  for (( j = 0 ; j < 1 ; j++ )); do # for each mass 
	# create the datacard and the workspaces
	#check the name on R2GGBBFitter.cc!!   legmc->SetHeader("300 GeV | CIC + X jets selection");
        # TString ssignal   = "./MiniTrees/jetwin/WP1/1000.root";
        # TString ddata   = "./MiniTrees/jetwin/WP1/Data_1000.root";
	# things to choose on the runcard
	# the signal file
	sed -i -r -e "s/WP[0-9]/WP$i/g" R2GGBBFitter_mtot.cc
        sed -i -r -e "s/m[0-9]+.root/m${radion[$j]}.root/g" R2GGBBFitter_mtot.cc
        sed -i -r -e "s/m[0-9]+_8TeV_nm_m[0-9]+.root/m${radion[$j]}_8TeV_nm_m${radion[$j]}.root/g" R2GGBBFitter_mtot.cc
#	sed -i -r -e "s/[0-9]+\_minimal.root/${radion[$j]}\_minimal.root/g" R2GGBBFitter_mtot.cc
#	sed -i -r -e "s/[0-9]+\_regression/${radion[$j]}\_regression/g" R2GGBBFitter_mtot.cc
#	sed -i -r -e "s/m[0-9]+\_/m${radion[$j]}\_/g" R2GGBBFitter_mtot_side.cc
#	sed -i -r -e "s/m[0-9]+\_/m${radion[$j]}\_/g" models_mtot_exp.rs
	# the legend
	sed -i -r -e "s/[0-9]+ GeV\"\);/${radion[$j]} GeV\"\);/g" R2GGBBFitter_mtot.cc
	# the window mtot > 550 || mtot < 450
	sed -i -r -e "s/mtot > [0-9]+ \|\| mtot < [0-9]+/mtot > ${winu[$j]} \|\| mtot < ${winl[$j]}/g" R2GGBBFitter_mtot.cc
echo WP$i MR ${radion[$j]}
	# the mass to fit around
	sed -i -r -e "s/runfits\([0-9]+/runfits\(${radion[$j]}/g" runfits.C
	# the fit model
	# mtot_sig_m0[500.0, 450, 550];
	# mtot_sig_m0_cat0[500.0, 450, 550];
	# mtot_sig_m0_cat1[500.0, 450, 550];
	# mtot_sig_m0_cat2[500.0, 450, 550];
	sed -i -r -e "s/mtot_sig_m0\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_exp.rs
	sed -i -r -e "s/mtot_sig_m0_cat0\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0_cat0\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_exp.rs
	sed -i -r -e "s/mtot_sig_m0_cat1\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0_cat1\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_exp.rs
	sed -i -r -e "s/mtot_sig_m0_cat2\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0_cat2\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_exp.rs
	mkdir radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	root -l -q runfits.C >> radlim_CSV_WP$i/radlim${radion[$j]}_CSV/log_radlim${radion[$j]}.txt
	mv workspaces/hgg.* radlim_CSV_WP$i/radlim${radion[$j]}_CSV
	mv datacards/* radlim_CSV_WP$i/radlim${radion[$j]}_CSV
	#also colect the plots
	mv databkgoversig*  radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	mv sigmodel*  radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	mv remenber.txt radlim${radion[$j]}_CSV
	## create limits root files for each mass
#	cd radlim${radion[$j]}_CSV
	cd radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	#
#	combine hgg.mH${radion[$j]}.0_8TeVlnu.txt -M Asymptotic -S 0 >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_lnu.txt
#	mv higgsCombineTest.Asymptotic.mH${radion[$j]}.root higgsCombineTest.Asymptotic.mR${radion[$j]}_lnu.root
#	echo did with lnu 
	#
	#text2workspace.py hgg.mH${radion[$j]}.0_8TeVonecat.txt -o hgg.mH${radion[$j]}.0_8TeVonecat.root -L ../../GaussExp_cxx.so
        combine -M Asymptotic hgg.mH${radion[$j]}.0_8TeVonecat.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}onecat.txt -L ../../GaussExp_cxx.so
	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mR${radion[$j]}_onecat.root
	echo did with rep 2btag only
	#
	#text2workspace.py hgg.mH${radion[$j]}.0_8TeVrep.txt -o hgg.mH${radion[$j]}.0_8TeVrep.root -L ../../GaussExp_cxx.so
	#combine -M Asymptotic hgg.mH${radion[$j]}.0_8TeVrep.txt -S 0  -L ../../GaussExp_cxx.so >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_nosyst.txt
	combine -M Asymptotic hgg.mH${radion[$j]}.0_8TeVrep.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_nosyst.txt
	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mR${radion[$j]}_nosyst.root
	echo with rep no syst	
	#combine -M Asymptotic hgg.mH${radion[$j]}.0_8TeVrep.txt -L ../../GaussExp_cxx.so >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_reparametrization.txt
	combine -M Asymptotic hgg.mH${radion[$j]}.0_8TeVrep.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_reparametrization.txt
	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.root
	echo did with rep 
	# lnu do not work in 4 body
#	text2workspace.py hgg.mH${radion[$j]}.0_8TeVlnu.txt -o hgg.mH${radion[$j]}.0_8TeVlnu.root -L ../../GaussExp_cxx.so
#	combine -M Asymptotic hgg.mH${radion[$j]}.0_8TeVlnu.root -L ../../GaussExp_cxx.so >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_lnu.txt
#	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}lnu.root
#	echo did with lnu 
	cd ../..
  done # mass
#  sed -i -r -e "s/WP[0-9]/WP$i/g" Brazilianflag.cc
#  root -l -q Brazilianflag.cc
done # wp
#sed -i -r -e "s/mR300.Xcut/Mcut_cutbased/g" Brazilianflag.cc #adapt this name to the cut
#root -l Brazilianflag.cc
