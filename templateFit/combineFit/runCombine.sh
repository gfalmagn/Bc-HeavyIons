COL=${1:-"pp"} #"PbPb" or "pp"
nominalOnly=false
BASENAME=""
BASENAME+="_"${COL}"_2bins"
HISTSYST="" #meta-fit systematic, linked to the changes in InputForCombine.root
FITSYST="" #variations of the fit method, but start with the same histograms
BLINDY="" #whether to pretend having blinded or unblinded yields (default is unblinded yields, with weight*4 for blinded dataset)
FULLNAME() {
    echo ${BASENAME}${HISTSYST}${FITSYST}
}

ln -s ../../../../../InputForCombine_${COL}.root InputForCombine_${COL}.root
ln -s ../../../../../InputForCombine_${COL}_BDTbinsUp.root InputForCombine_${COL}_BDTbinsUp.root
ln -s ../../../../../InputForCombine_${COL}_BDTbinsDown.root InputForCombine_${COL}_BDTbinsDown.root
ln -s ../../../../../InputForCombine_${COL}_MbinsVar1.root InputForCombine_${COL}_MbinsVar1.root
ln -s ../../../../../InputForCombine_${COL}_MbinsVar2.root InputForCombine_${COL}_MbinsVar2.root
ln -s ../../../../../InputForCombine_${COL}_BDTuncorrFromM.root InputForCombine_${COL}_BDTuncorrFromM.root
ln -s ../../../../../InputForCombine_${COL}_regulLowStatShapes.root InputForCombine_${COL}_regulLowStatShapes.root
ln -s ../../../../../InputForCombine_${COL}_scaleSystBDTintegrated_regulLowStatShapes.root InputForCombine_${COL}_scaleSystBDTintegrated_regulLowStatShapes.root

makeWorkspace() {
    text2workspace.py datacard${BASENAME}.txt -o datacard${BASENAME}${HISTSYST}${BLINDY}.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*Kin1/BcSig:r1[1,0,6]' --PO 'map=.*Kin2/BcSig:r2[1,0,6]' --channel-masks  --keyword-value METAFITSYST=${HISTSYST} --keyword-value BLIND=${BLINDY} -v 0 > /dev/null
}
runCombineBase() {
    echo ValidateDatacards.py datacard${BASENAME}${HISTSYST}.root \
    \&\& combine -d datacard${BASENAME}${HISTSYST}.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n $(FULLNAME) --saveWithUncertainties -v 0
}
printSyst() {
    echo -e "\n*********\n   "${HISTSYST}${FITSYST}" "${BLINDY}"\n********\n"
}

makeWorkspace
echo -e "\n*********\n   nominal\n********\n"
#echo "$(runCombineBase) --minos all"
eval "$(runCombineBase) --minos all"

if [ "$nominalOnly" = false ] ; then
    FITSYST="_noBDT1"
    printSyst
    #echo "$(runCombineBase) --setParameters mask_BDT1Kin1=1,mask_BDT1Kin2=1"
    eval "$(runCombineBase) --setParameters mask_BDT1Kin1=1,mask_BDT1Kin2=1"

    if [ $COL = "PbPb" ] ; then
	FITSYST="" #for GOF test with blinded yields
	BLINDY="_blindYields"
	printSyst
	makeWorkspace
	BLINDY=""
    fi

    HISTSYST="_BDTuncorrFromM"
    makeWorkspace
    FITSYST=""
    printSyst
    eval "$(runCombineBase)"
    FITSYST="_noBDT1"
    printSyst
    eval $(runCombineBase) --setParameters mask_BDT1Kin1=1,mask_BDT1Kin2=1
    
    HISTSYST="_regulLowStatShapes"
    makeWorkspace
    FITSYST="_autoMCstatsNoBDT3"
    printSyst
    eval $(runCombineBase) --freezeParameters 'rgx{prop_binBDT3Kin.*_bin.*}'
    FITSYST="_autoMCstatsNoBDT23"
    printSyst
    eval $(runCombineBase) --freezeParameters 'rgx{prop_binBDT3Kin.*_bin.*}','rgx{prop_binBDT2Kin.*_bin.*}'
    
    HISTSYST="_scaleSystBDTintegrated_regulLowStatShapes"
    makeWorkspace
    FITSYST="_autoMCstatsNoBDT3"
    printSyst
    eval $(runCombineBase) --freezeParameters 'rgx{prop_binBDT3Kin.*_bin.*}'
    FITSYST="_autoMCstatsNoBDT23"
    printSyst
    eval $(runCombineBase) --freezeParameters 'rgx{prop_binBDT3Kin.*_bin.*}','rgx{prop_binBDT2Kin.*_bin.*}'

    FITSYST=""
    HISTSYST="_MbinsVar1"
    makeWorkspace
    printSyst
    eval $(runCombineBase)
    HISTSYST="_MbinsVar2"
    makeWorkspace
    printSyst
    eval $(runCombineBase)
    HISTSYST="_BDTbinsUp"
    makeWorkspace
    printSyst
    eval $(runCombineBase)
    HISTSYST="_BDTbinsDown"
    makeWorkspace
    printSyst
    eval $(runCombineBase)
fi


#comment this for reference
: <<'END'

combine -M GoodnessOfFit --algorithm saturated datacard_pp_2bins.root
combine -M GoodnessOfFit --algorithm saturated datacard_pp_2bins.root --toysFrequentist -t 500 --seed 2235 -n _pp

combine -M GoodnessOfFit --algorithm saturated datacard_PbPb_2bins_blindYields.root
combine -M GoodnessOfFit --algorithm saturated datacard_PbPb_2bins_blindYields.root --toysFrequentist -t 500 --seed 2235 -n _PbPb

combine -M Significance --signif datacard_pp.root #observed
combine -M Significance datacard_PbPb.root -t -1 --expectSignal=0.57 --toysFreq #expected, post-fit (otherwise the prefit shapes/norms are too wrong)

combine -d datacard_PbPb_2bins.root -M FitDiagnostics --saveNormalizations -n _PbPb_2bins_wSigToys --saveWithUncertainties --toysFreq -t 300 --setParameters r1=1.1,r2=0.46 --seed 999
combine -d datacard_pp_2bins.root -M FitDiagnostics --saveNormalizations -n _pp_2bins_wSigToys --saveWithUncertainties --toysFreq -t 300 --setParameters r1=0.85,r2=1.15 --seed 999

combineTool.py -M Impacts -d datacard_pp_2bins.root --doInitialFit -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_pp_2bins.root --doFits -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_pp_2bins.root -m 125 -o impacts_pp_2bins.json
plotImpacts.py -i impacts_pp_2bins.json -o impacts_pp_2bins

combineTool.py -M Impacts -d datacard_PbPb_2bins.root --doInitialFit -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_PbPb_2bins.root --doFits -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_PbPb_2bins.root -m 125 -o impacts_PbPb_2bins.json
plotImpacts.py -i impacts_PbPb_2bins.json -o impacts_PbPb_2bins

combine -M MultiDimFit datacard_PbPb_2bins.root --algo grid --points 625 --setParameterRanges r1=0,2.8:r2=0,1.1 -n _PbPb
combine -M MultiDimFit datacard_pp_2bins.root --algo grid --points 625 --setParameterRanges r1=0.45,1.25:r2=0.8,1.6 -n _pp

combine -d  datacard_pp_2bins.root -M FitDiagnostics -t -1 --expectSignal 0 --setParameters r1=0,r2=0 -n _pp_zeroSigCheck #r1 r2 should be 0
combine -d  datacard_PbPb_2bins.root -M FitDiagnostics -t -1 --expectSignal 0 --setParameters r1=0,r2=0 -n _PbPb_zeroSigCheck
combine -d  datacard_pp_2bins.root -M FitDiagnostics -t -1 --setParameters r1=1,r2=1 -n _pp_Sig1Check #r1 r2 should be 1
combine -d  datacard_PbPb_2bins.root -M FitDiagnostics -t -1 --setParameters r1=1,r2=1 -n _PbPb_Sig1Check

END

