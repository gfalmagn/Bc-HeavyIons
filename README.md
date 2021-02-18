* `rt` stands for `root -l`
* When the argument of a script is true/false, it means it must be run for both values, corresponding to pp or PbPb
* `cd ~` means returning to the base directory

## Produce pre-selected samples
- NB: No need to select on centrality (negligible effect of centrality bin variation)
- Normalise/weight (non)prompt Jpsi MC:
  * `rt -b "GetNormalization.C(true/false,true/false)"'
- Do pre-selection:
  * `cd ~/BDT/`
  * `rt "MakeInputTrees(true/false)"`
  * `cp BDT_InputTree_pp.root BDT_InputTree_pp_copy.root`
- Add a weight for the candidates with ambiguous Jpsi dimuon choice (2 OS pairs):
  * `rt "addJpsiChoiceWeight.C(true/false,false)"`
  * `cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep1.root`
- Train BDT:
  * `rt "ClassifierSigBkg.C(true/false,false)" > BDT_trainingOutput_pp/PbPb.txt`
  * `rt "addBDTvariable.C(true/false,false)"`
  * `cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep2.root`
- Add the Jpsi choice weight, with probabilities that depend on the BDT value (purer dimuon mass distro when looking at high BDT):
  * `rt "addJpsiChoiceWeight.C(true/false,true)"`
  * `cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep3.root`
- Check dimuon mass for Jpsi choice weight:
  * `cd ~/JpsiMass/`
  * `rt "drawJpsiMass.C(true/false)"`


## Preliminaries (BDT-linked)
- Determine a less powerful version of BDT variable, but that is uncorrelated from the mass. Actually, should do this AFTER first nominal fit, for extracting postfit normalisations
  * `cd ~/BDT/`
  * `rt -b "UncorrelateBDTfromM.C(false)"' (loophole: technically, needs the input of the postfit yields, even in the first step. Works because needed only for systematics)
- Determine the BDT binning (written in a header file), for each strategy and kinematic bin. The first time, when the header file does not exist yet, have to run it with option `firstTime=true`, before the second addJpsiChoiceWeight.C running 
  * `rt "DetermineBDTcuts.C(false)"` (use true in second argument for the first time you run this)


## Acceptance and efficiency, first step (no pt spectrum correction)
- Make the acceptance and efficiency maps. Only needs the pre-selected samples, and the BDT binning. Can be run with BDTuncorrelatedFromM=true (changes BDT binning, hence BDT efficiency), to obtain efficiency of being in some BDT bin
  * `cd ~/acceptance/`
  * `rt -b "BuildAcceptanceMap.C(false,false)"`
  * `cd ~/efficiency/`
  * `rt -b "runBuildEffMap.C(false,false)"`

## Template fitting 
- Make the templates histos for combine. Includes histos with acc eff corrections
  * `cd ~/templateFit`
  * `rt "HistsForCombine.C(true/false,false,false)"`
- Make the template fit, without acc eff corrections
  * `cd ~/templateFit`
  * `cmsrel CMSSW_10_3_4; cd CMSSW_10_3_4/src'
  * `cmsenv'
  * `git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit'
  * `cp ../../FitDiagnostics_modified.cc HiggsAnalysis/CombinedLimit/src/FitDiagnostics.cc; cp ../../FitDiagnostics_modified.h HiggsAnalysis/CombinedLimit/interface/FitDiagnostics.h' (to record r2 in the toys tree in FitDiagnostics)
  * `cd HiggsAnalysis/CombinedLimit/test'
  * `ln -s /home/llr/cms/falmagne/Bc/templateFit/runCombine.sh runCombine.sh'
  * `source runCombine.sh "pp"'
  * `source runCombine.sh "PbPb"'
- Plot postfit trimuon mass shapes. The two string options correspond different fit strategies (metafit systematic variations)
  *`rt -b "plotCombineOutput.C(true/false,false,false,\"\",\"\")"'

## Yields + AccEff corrections. Correction of pT distribution of MC with corrected yield.
- Correct yields with event-by-event acc eff corrections, then store+draw them
  * `cd ~/AccEffCorr/'
  * `rt -b "runMetafitSysts.C(false)"'
- Draw all meta-fit variations and extract a systematic (with first-step AccEff)
  * `rt "Draw_metafitSyst.C(false)"'
  * `rt -b "Draw_corrYields.C(true/false,false)"'
- Make toy biases for pt spectrum correction of MC
  * `cd ~/twoSteps/'
  * `rt "MakeToyBiases.C(false)"'

## 2nd step: Run BDT, accXeff and fit, after correcting the pT distribution of signal MC
  * `cd ~/BDT/`   
- Add pT weighting on signal MC, from 2nd step
  * `rt "add2ndStepWeight.C(true/false)"`
  * `cp BDT_InputTree_pp/PbPb.root BDT_InputTree_pp/PbPb_copystep4.root`
- Train BDT + mass decorrelation
  * `rt "ClassifierSigBkg.C(true/false,true)" > BDT_trainingOutput_2ndStep_pp/PbPb.txt`
  * `rt "addBDTvariable.C(true/false,true)"`
  * `cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep5.root`
  * `rt "DetermineBDTcuts.C(true)"`
  * `rt -b "UncorrelateBDTfromM.C(true)"'
  * `rt "DetermineBDTcuts.C(true)"` (to update for decorrelated BDT bins)
- Run acceptance and efficiency, from 1st-step pT-corrected MC
  * `cd ~/acceptance/`
  * `rt -b "BuildAcceptanceMap.C(false,true)"`
  * `cd ~/efficiency/`
  * `rt -b "runBuildEffMap.C(false,true)"`
- Preliminary fit without BDT weighting, for later BDT weighting (pp only)
  * `cd ~/templateFit`
  * `rt "HistsForCombine.C(true/false,true,false)"' (doing main PbPb histos here)
  * `cd CMSSW_10_3_4/src/HiggsAnalysis/CombinedLimit/test'
  * `source runCombine.sh "pp"/"PbPb" "_2ndStep" true'
- Check+weight the BDT distribution of the summed templates
  * `cd ~/BDT/weighting'
  * `rt -b "BDTweighting.C(true)"'
- Run final (second-step) fit
  * `cd ~/templateFit`
  * `rt "HistsForCombine.C(true,true,true)"' (pp histos including the BDT distro weights)
  * `source runCombine.sh "pp"/"PbPb" "_2ndStep"'
  * Run fit tests commented in runCombine.sh
- Plot postfit trimuon mass shapes, second step
  *`rt -b "plotCombineOutput.C(true/false,true,false,\"\",\"\")"'

##Record yields, metafit systematics, deploy AccEff methods
- Correct yields with event-by-event acc eff corrections + record postfit yields
  * `cd ~/AccEffCorr/'
  * `rt -b "runMetafitSysts.C(true)"'
- Draw all meta-fit variations and extract a systematic (with first-step AccEff). Record 2nd-step corrected yield.
  * `rt "Draw_metafitSyst.C(true)"'
  * `rt -b "Draw_corrYields.C(true/false,true)"'
- Make toy biases for pt spectrum correction of MC. Run also the version for plot only (less variations are drawn).
  * `cd ~/twoSteps/'
  * `rt "MakeToyBiases.C(true)"'
  * `rt "MakeToyBiases.C(true,true)"'
- 3rd-step nominal+toys acceptance and efficiency
  * `cd ~/acceptance/`
  * `rt -b "BuildAcceptanceMap.C(true,true)"`
  * `cd ~/efficiency/`
  * `rt -b "runBuildEffMap.C(true,true)"`
- Calculate 2-steps acceptance x efficiency correction, and systematic from toy corrections
  * `cd ~/twoSteps/'
  * `rt "plotToyAccEff.C(true/false,true)"'


## Results+Systematics
- Compare corrected yields, including with 3rd step AccEff:
  * `cd ~/AccEffCorr/'
  * `rt -b "Draw_corrYields.C(true/false,true,true)"'
- Compute correlation factors for RAA from the ones for pp and PbPb
  * `cd ~/RAA/correlations/`
  * `rt -b "RAAfitErrCorr.C(true/false)"'
- Draw R_PbPb:
  * `rt Draw_XSandRAA.C'

## Automatic latex tables
- N-1 efficiencies (preselection section 3.4)
  * `cd ~/Bc/efficiency/'
  * `rt "Nmin1efficiencies.C(true/false)"
  * `python Nmin1efficiencies_texTable.py'
- N-1 efficiencies (preselection section 3.4)
  * `cd ~/Bc/RAA/'
  * `python XSRAA_texTable.py'
- BDT performance and overtraining
  * `cd ~/Bc/BDT/'
  * `python BDTperf_texTable.py'

## Various plottings
  * `cd ~/Bc/templateFit/'
  * `rt -b "drawFitChecks.C(true)"'
  * `cd ~/Bc/fakeJpsi/'
  * `python sideband_fitting.py'
  * `rt -b "drawFakeJpsi.C(true/false,true)"'
  * `cd ~/Bc/flipJpsi/'
  * `rt -b "drawFlipJpsi.C(true/false,true)"'
  * `cd ~/BDT/drawVars/'
  * `rt -b "drawVariables.C(true/false,true)"'
