COL=${1:-"pp"} #"PbPb" or "pp"
nominalOnly=false
BASENAME=""
if [ $COL = "PbPb" ]
then 
    BASENAME+="_"${COL}"_2bins"
fi
if [ $COL = "pp" ]
then 
    BASENAME+="_"${COL}"_2bins"
fi
HISTSYST="" #meta-fit systematic, linked to the changes in InputForCombine.root
FITSYST="" #variations of the fit method, but start with the same histograms
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
    text2workspace.py datacard${BASENAME}.txt -o datacard${BASENAME}${HISTSYST}.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*Kin1/BcSig:r1[1,0,6]' --PO 'map=.*Kin2/BcSig:r2[1,0,6]' --channel-masks  --keyword-value METAFITSYST=${HISTSYST} -v 0 > /dev/null
}
runCombineBase() {
    echo combine -d datacard${BASENAME}${HISTSYST}.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n $(FULLNAME) --saveWithUncertainties -v 0
}
printSyst() {
    echo -e "\n*********\n   "${HISTSYST}${FITSYST}"\n********\n"
}

makeWorkspace
echo -e "\n*********\n   nominal\n********\n"
eval $(runCombineBase) --minos all

if [ "$nominalOnly" = false ] ; then
    FITSYST="_noBDT1"
    printSyst
    eval $(runCombineBase) --setParameters mask_BDT1Kin1=1,mask_BDT1Kin2=1
    
    HISTSYST="_BDTuncorrFromM"
    makeWorkspace
    FITSYST=""
    printSyst
    eval $(runCombineBase)
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
combine -M GoodnessOfFit --algorithm saturated datacard_pp_2bins.root --seed 1111
combine -M GoodnessOfFit --algorithm saturated datacard_pp_2bins.root --toysFrequentist -t 500 --seed 2235 -n _pp

combine -M Significance --signif datacard_pp.root #observed
combine -M Significance datacard_PbPb.root -t -1 --expectSignal=0.57 --toysFreq #expected, post-fit (otherwise the prefit shapes/norms are too wrong)

combine -d datacard_PbPb_2bins.root -M FitDiagnostics --saveNormalizations -n _PbPb_2bins_wSigToys --saveWithUncertainties --toysFreq -t 300 --setParameters r1=1.15,r2=0.6 --seed 999
combine -d datacard_pp_2bins.root -M FitDiagnostics --saveNormalizations -n _pp_2bins_wSigToys --saveWithUncertainties --toysFreq -t 300 --setParameters r1=0.85,r2=1.15 --seed 999

combineTool.py -M Impacts -d datacard_BDTuncorrFromM_PbPb.root --doInitialFit -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_BDTuncorrFromM_PbPb.root --doFits -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_BDTuncorrFromM_PbPb.root -m 125 -o impacts_BDTuncorrFromM_PbPb.json
plotImpacts.py -i impacts_BDTuncorrFromM_PbPb.json -o impacts_BDTuncorrFromM_PbPb

combine -M MultiDimFit datacard_PbPb_2bins.root --algo grid --points 625 --setParameterRanges r1=0,2.7:r2=0,1.2 -n _PbPb
combine -M MultiDimFit datacard_pp_2bins.root --algo grid --points 625 --setParameterRanges r1=0.5,1.3:r2=0.7,1.5 -n _pp

END

