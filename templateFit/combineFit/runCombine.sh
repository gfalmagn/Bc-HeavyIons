COL="PbPb"
basename=""
if [ $COL = "PbPb" ]
then 
    basename+="_NonPromptJpsi_PromptJpsi_"${COL}"_2bins"
fi
if [ $COL = "pp" ]
then 
    basename+="_NonPromptJpsi_flipJpsi_"${COL}"_2bins"
fi

ln -s ../../../../../InputForCombine_${COL}.root InputForCombine_${COL}.root
ln -s ../../../../../InputForCombine_${COL}_BDTuncorrFromM.root InputForCombine_${COL}_BDTuncorrFromM.root
ln -s ../../../../../InputForCombine_${COL}_regulLowStatShapes.root InputForCombine_${COL}_regulLowStatShapes.root
ln -s ../../../../../InputForCombine_${COL}_scaleSystBDTintegrated_regulLowStatShapes.root InputForCombine_${COL}_scaleSystBDTintegrated_regulLowStatShapes.root

#text2workspace.py datacard_${COL}_2bins.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO 'map=.*Kin1/BcSig:r1[1,0,6]' --PO 'map=.*Kin2/BcSig:r2[1,0,6]' --channel-masks  --keyword-value METAFITSYST="" -v 0
#${COL}\ nominal
echo -e "\n*********\n   nominal\n********\n"
text2workspace.py datacard_${COL}_2bins.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*Kin1/BcSig:r1[1,0,6]' --PO 'map=.*Kin2/BcSig:r2[1,0,6]' --channel-masks  --keyword-value METAFITSYST="" -v 0
combine -d datacard_${COL}_2bins.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n ${basename} --saveWithUncertainties -v 0
combine -d datacard_${COL}_2bins.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n ${basename}_noBDT1 --saveWithUncertainties --setParameters mask_BDT1Kin1=1,mask_BDT1Kin2=1 -v 0

#${COL}\ _BDTuncorrFromM
echo -e "\n*********\n   _BDTuncorrFromM\n********\n"
text2workspace.py datacard_${COL}_2bins.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*Kin1/BcSig:r1[1,0,6]' --PO 'map=.*Kin2/BcSig:r2[1,0,6]' --channel-masks  --keyword-value METAFITSYST="_BDTuncorrFromM" -v 0
combine -d datacard_${COL}_2bins.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n ${basename}_BDTuncorrFromM --saveWithUncertainties
combine -d datacard_${COL}_2bins.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n ${basename}_BDTuncorrFromM_noBDT1 --saveWithUncertainties --setParameters mask_BDT1Kin1=1,mask_BDT1Kin2=1

#${COL}\ _regulLowStatShapes
echo -e "\n*********\n   _regulLowStatShapes\n********\n"
text2workspace.py datacard_${COL}_2bins.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*Kin1/BcSig:r1[1,0,6]' --PO 'map=.*Kin2/BcSig:r2[1,0,6]' --channel-masks  --keyword-value METAFITSYST="_regulLowStatShapes" -v 0
combine -d datacard_${COL}_2bins.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n ${basename}_regulLowStatShapes_autoMCstatsNoBDT3 --saveWithUncertainties --freezeParameters 'rgx{prop_binBDT3Kin.*_bin.*}' -v 0
combine -d datacard_${COL}_2bins.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n ${basename}_regulLowStatShapes_autoMCstatsNoBDT23 --saveWithUncertainties --freezeParameters 'rgx{prop_binBDT3Kin.*_bin.*}','rgx{prop_binBDT2Kin.*_bin.*}' -v 0

#${COL}\ _scaleSystBDTintegrated_regulLowStatShapes
echo -e "\n*********\n   _scaleSystBDTintegrated_regulLowStatShapes\n********\n"
text2workspace.py datacard_${COL}_2bins.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*Kin1/BcSig:r1[1,0,6]' --PO 'map=.*Kin2/BcSig:r2[1,0,6]' --channel-masks  --keyword-value METAFITSYST="_scaleSystBDTintegrated_regulLowStatShapes" -v 0
combine -d datacard_${COL}_2bins.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n ${basename}_scaleSystBDTintegrated_regulLowStatShapes_autoMCstatsNoBDT3 --saveWithUncertainties --freezeParameters 'rgx{prop_binBDT3Kin.*_bin.*}' -v 0
combine -d datacard_${COL}_2bins.root -M FitDiagnostics --saveNormalizations --saveOverallShapes --saveShapes -n ${basename}_scaleSystBDTintegrated_regulLowStatShapes_autoMCstatsNoBDT23 --saveWithUncertainties --freezeParameters 'rgx{prop_binBDT3Kin.*_bin.*}','rgx{prop_binBDT2Kin.*_bin.*}' -v 0

#comment this for reference
: <<'END'
combine -M GoodnessOfFit --algorithm saturated datacard_pp.root --seed 1234567
combine -M GoodnessOfFit --algorithm saturated datacard_pp.root --toysFrequentist -t 1000 --seed 2233 --setParameters mask_BDT0=1

combine -M Significance --signif datacard_pp.root #observed
combine -M Significance datacard_PbPb.root -t -1 --expectSignal=0.57 --toysFreq #expected, post-fit (otherwise the prefit shapes/norms are too wrong)

combine -d datacard_BDTuncorrFromM_PbPb.root -M FitDiagnostics -n _NonPromptJpsi_PromptJpsi_BDTuncorrFromM_PbPb_toysWsig --saveWithUncertainties --toysFreq -t 400 --expectSignal 0.7
combine -d datacard_BDTuncorrFromM_PbPb.root -M FitDiagnostics -n _NonPromptJpsi_PromptJpsi_BDTuncorrFromM_PbPb_toysNoSig --saveWithUncertainties --toysFreq -t 400

combineTool.py -M Impacts -d datacard_BDTuncorrFromM_PbPb.root --doInitialFit -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_BDTuncorrFromM_PbPb.root --doFits -m 125 --robustFit 1
combineTool.py -M Impacts -d datacard_BDTuncorrFromM_PbPb.root -m 125 -o impacts_BDTuncorrFromM_PbPb.json
plotImpacts.py -i impacts_BDTuncorrFromM_PbPb.json -o impacts_BDTuncorrFromM_PbPb

combine -M MultiDimFit datacard_PbPb_2bins.root --algo grid --points 484 --setParameterRanges r1=0,2.5:r2=0,1.3
combine -M MultiDimFit datacard_pp_2bins.root --algo grid --points 484 --setParameterRanges r1=0.5,1.5:r2=0.5,1.5
END

