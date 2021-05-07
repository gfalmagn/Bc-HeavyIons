COL=${1:-"pp"} #"PbPb" or "pp"
SECONDSTEP=${2:-""} #"" or "_2ndStep"
nominalOnly=${3:-false}
BINS=${4:-"_2bins"} #or "_integrated" or "_centBins"
BASENAME=""
BASENAME+="_"${COL}${BINS}
HISTSYST="" #meta-fit systematic, linked to the changes in InputForCombine.root
FITSYST="" #variations of the fit method, but start with the same histograms
BLINDY="" #whether to pretend having blinded or unblinded yields (default is unblinded yields, with weight*4 for blinded dataset)
FULLNAME() {
    echo ${BASENAME}${SECONDSTEP}${HISTSYST}${FITSYST}
}

######## ONLY the first time
#ln -s ../../../../../InputForCombine_${COL}${SECONDSTEP}.root InputForCombine_${COL}${SECONDSTEP}.root
#ln -s ../../../../../InputForCombine_${COL}${SECONDSTEP}_BDTbinsUp.root InputForCombine_${COL}${SECONDSTEP}_BDTbinsUp.root
#ln -s ../../../../../InputForCombine_${COL}${SECONDSTEP}_BDTbinsDown.root InputForCombine_${COL}${SECONDSTEP}_BDTbinsDown.root
#ln -s ../../../../../InputForCombine_${COL}${SECONDSTEP}_MbinsVar1.root InputForCombine_${COL}${SECONDSTEP}_MbinsVar1.root
#ln -s ../../../../../InputForCombine_${COL}${SECONDSTEP}_MbinsVar2.root InputForCombine_${COL}${SECONDSTEP}_MbinsVar2.root
#ln -s ../../../../../InputForCombine_${COL}${SECONDSTEP}_BDTuncorrFromM.root InputForCombine_${COL}${SECONDSTEP}_BDTuncorrFromM.root
#ln -s ../../../../../InputForCombine_${COL}${SECONDSTEP}_regulLowStatShapes.root InputForCombine_${COL}${SECONDSTEP}_regulLowStatShapes.root
#ln -s ../../../../../InputForCombine_${COL}${SECONDSTEP}_scaleSystBDTintegrated_regulLowStatShapes.root InputForCombine_${COL}${SECONDSTEP}_scaleSystBDTintegrated_regulLowStatShapes.root

#pT or centrality bins, or integrated
BIN1=""
BIN2=""
POIs="--PO map=.*"${BIN1}"/BcSig:r1[1,0,6]"
if [ $BINS = "_2bins" ] || [ $BINS = "_centBins" ] ; then
    if [ $BINS = "_2bins" ] ; then
	BIN1="Kin1"
	BIN2="Kin2"
    elif [ $BINS = "_centBins" ] ; then
	BIN1="Cent1"
	BIN2="Cent2"
    fi
    POIs="--PO map=.*"${BIN1}"/BcSig:r1[1,0,6] --PO map=.*"${BIN2}"/BcSig:r2[1,0,6]"
fi

makeWorkspace() {
    text2workspace.py datacard${BASENAME}.txt -o datacard${BASENAME}${SECONDSTEP}${HISTSYST}${BLINDY}.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel $POIs --channel-masks --keyword-value METAFITSYST=${SECONDSTEP}${HISTSYST} --keyword-value BLIND=${BLINDY} -v 0 > /dev/null
}
runCombineBase() {
    echo ValidateDatacards.py datacard${BASENAME}${SECONDSTEP}${HISTSYST}.root \
     \&\& combine -d datacard${BASENAME}${SECONDSTEP}${HISTSYST}.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n $(FULLNAME) --saveWithUncertainties --freezeParameters ${FREEZE} -v 0
}
printSyst() {
    echo -e "\n*********\n   "${SECONDSTEP}${HISTSYST}${FITSYST}" "${BLINDY}"\n********\n"
}

if [ $COL = "pp" ] ; then
    FREEZE_OR="misIDinit,highMassInit"
else
    FREEZE_OR="misIDinit"
fi
FREEZE=${FREEZE_OR}

makeWorkspace
echo -e "\n*********\n   nominal\n********\n"
eval "$(runCombineBase) --minos all"

if [ "$nominalOnly" = false ] ; then
    FITSYST="_noBDT1"
    printSyst
    eval "$(runCombineBase) --setParameters mask_BDT1"${BIN1}"=1,mask_BDT1"${BIN2}"=1"

    #    if [ $COL = "PbPb" ] ; then
    #        FITSYST="" #for GOF test with blinded yields
    #	BLINDY="_blindYields"
    #	printSyst
    #	makeWorkspace
    #	BLINDY=""
    #    fi

    HISTSYST="_BDTuncorrFromM"
    makeWorkspace
    FITSYST=""
    printSyst
    eval "$(runCombineBase)"
    FITSYST="_noBDT1"
    printSyst
    eval "$(runCombineBase) --setParameters mask_BDT1"${BIN1}"=1,mask_BDT1"${BIN2}"=1"

    HISTSYST="_regulLowStatShapes"
    makeWorkspace
    FITSYST="_autoMCstatsNoBDT3"
    FREEZE=${FREEZE_OR}",rgx{prop_binBDT3.*_bin.*}"
    printSyst
    eval $(runCombineBase)
    FITSYST="_autoMCstatsNoBDT23"
    FREEZE=${FREEZE_OR}",rgx{prop_binBDT3.*_bin.*},rgx{prop_binBDT2.*_bin.*}"
    printSyst
    eval $(runCombineBase)

    HISTSYST="_scaleSystBDTintegrated_regulLowStatShapes"
    makeWorkspace
    FITSYST="_autoMCstatsNoBDT3"
    FREEZE=${FREEZE_OR}",rgx{prop_binBDT3.*_bin.*}"
    if [ "$SECONDSTEP" = true -a $COL = "PbPb" -a $BINS="_2bins"] ; then
	FREEZE=${FREEZE}",prop_binBDT2Kin2_bin4_FakeJpsi"
    fi
    printSyst
    eval $(runCombineBase)
    FITSYST="_autoMCstatsNoBDT23"
    FREEZE=${FREEZE_OR}",rgx{prop_binBDT3.*_bin.*},rgx{prop_binBDT2.*_bin.*}"
    printSyst
    eval $(runCombineBase)
    FREEZE=${FREEZE_OR}

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

combine -M GoodnessOfFit --algorithm saturated datacard_pp_2bins_2ndStep.root
combine -M GoodnessOfFit --algorithm saturated datacard_pp_2bins_2ndStep.root --toysFrequentist -t 500 --seed 2235 -n _pp

combine -M GoodnessOfFit --algorithm saturated datacard_PbPb_2bins_2ndStep.root
combine -M GoodnessOfFit --algorithm saturated datacard_PbPb_2bins_2ndStep.root --toysFrequentist -t 500 --seed 2235 -n _PbPb

combine -M Significance --signif datacard_pp.root #observed
combine -M Significance datacard_PbPb.root -t -1 --expectSignal=0.57 --toysFreq #expected, post-fit (otherwise the prefit shapes/norms are too wrong)

combine -d datacard_PbPb_2bins_2ndStep.root -M FitDiagnostics --saveNormalizations -n _PbPb_2bins_wSigToys --saveWithUncertainties --toysFreq -t 300 --setParameters r1=1.,r2=1. --seed 999
combine -d datacard_pp_2bins_2ndStep.root -M FitDiagnostics --saveNormalizations -n _pp_2bins_wSigToys --saveWithUncertainties --toysFreq -t 300 --setParameters r1=1.,r2=1. --seed 999

combineTool.py -M Impacts -d datacard_pp_2bins_2ndStep.root --doInitialFit -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_pp_2bins_2ndStep.root --doFits -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_pp_2bins_2ndStep.root -m 125 -o impacts_pp_2bins.json
plotImpacts.py -i impacts_pp_2bins.json -o impacts_pp_2bins

combineTool.py -M Impacts -d datacard_PbPb_2bins_2ndStep.root --doInitialFit -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_PbPb_2bins_2ndStep.root --doFits -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_PbPb_2bins_2ndStep.root -m 125 -o impacts_PbPb_2bins.json
plotImpacts.py -i impacts_PbPb_2bins.json -o impacts_PbPb_2bins

combine -M MultiDimFit datacard_PbPb_2bins_2ndStep.root --algo grid --points 900 --setParameterRanges r1=0,3.:r2=0.15,2.3 -n _PbPb
combine -M MultiDimFit datacard_pp_2bins_2ndStep.root --algo grid --points 900 --setParameterRanges r1=0.75,1.65:r2=0.6,1.5 -n _pp

combine -d  datacard_pp_2bins_2ndStep.root -M FitDiagnostics -t -1 --expectSignal 0 --setParameters r1=0,r2=0 -n _pp_zeroSigCheck #r1 r2 should be 0
combine -d  datacard_PbPb_2bins_2ndStep.root -M FitDiagnostics -t -1 --expectSignal 0 --setParameters r1=0,r2=0 -n _PbPb_zeroSigCheck
combine -d  datacard_pp_2bins_2ndStep.root -M FitDiagnostics -t -1 --setParameters r1=1,r2=1 -n _pp_Sig1Check #r1 r2 should be 1
combine -d  datacard_PbPb_2bins_2ndStep.root -M FitDiagnostics -t -1 --setParameters r1=1,r2=1 -n _PbPb_Sig1Check

combineTool.py -M FastScan -w FINAL_WORKSPACE_NAME.root:w

#Here beware to run number of points n such that n*(mean proba) > 500   (or 1500, the number of wanted 2-steps toys)
combine -M MultiDimFit datacard_PbPb_2bins.root --algo random --points 3000 --setParameterRanges r1=0,2.3:r2=0.2,0.95 -n _PbPb_POIfromNLL
combine -M MultiDimFit datacard_PbPb_2bins_2ndStep.root --algo random --points 15000 --setParameterRanges r1=0,2.5:r2=0.4,1.9 -n _PbPb_POIfromNLL_2ndStep
combine -M MultiDimFit datacard_pp_2bins.root --algo random --points 3000 --setParameterRanges r1=0.75,1.3:r2=0.85,1.4 -n _pp_POIfromNLL
combine -M MultiDimFit datacard_pp_2bins_2ndStep.root --algo random --points 15000 --setParameterRanges r1=0.75,1.3:r2=0.85,1.4 -n _pp_POIfromNLL_2ndStep

END

