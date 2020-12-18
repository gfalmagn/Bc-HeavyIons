* `rt` stands for `root -l`
* Everything is shown for pp, but has to be reproduced with PbPb by changing the first option `ispp=true` (except for the scripts gathering pp and PbPb info)
* `cd ~` means returning to the base directory

## Produce pre-selected samples

- Do pre-selection:
  * `cd ~/BDT/`
  * `rt "MakeInputTrees(true,false)"`
  * `cp BDT_InputTree_pp.root BDT_InputTree_pp_copy.root`
- Add a weight for the candidates with ambiguous Jpsi dimuon choice (2 OS pairs):
  * `rt "addJpsiChoiceWeight.C(true,false)"`
  * `cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep1.root`
- Train BDT:
  * `rt "ClassifierSigBkg.C(true)"`
  * `rt "addBDTvariable.C(true)"`
  * `cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep2.root`
- Add the Jpsi choice weight, with probabilities that depend on the BDT value (purer dimuon mass distro when looking at high BDT):
  * `rt "addJpsiChoiceWeight.C(true,true)"`
  * `cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep3.root`

## Preliminaries (BDT-linked)
- Determine a less powerful version of BDT variable, but that is uncorrelated from the mass
  * `cd ~/BDT/`
  * `rt -b "UncorrelateBDTfromM.C(true)"'
- Determine the BDT binning (written in a header file), for each strategy and kinematic bin. The first time, when the header file does not exist yet, have to run it with option `firstTime=true`, before the second addJpsiChoiceWeight.C running 
  * `rt DetermineBDTcuts.C`
- Weight the BDT distribution of the [true Jpsi + muon from different vertex] background, i.e. flipJpsi or Prompt+(uncorrelated)Non-prompt MC, to data in the mass control region (or in the signal reigon)
  * `rt -b "BDTweighting.C(true)"'

## Acceptance and efficiency
- Make the acceptance and efficiency maps. Only needs the pre-selected samples, and the BDT binning. Need to run twice, with or without BDTuncorrelatedFromM (changes BDT binning, hence BDT efficiency)
  * `cd ~/acceptance/`
  * `rt BuildAcceptanceMap.C`
  * `cd ~/efficiency/`
  * `rt -b runBuildEffMap.C`

## Template fitting 
- Make the templates histos for combine. Includes histos with acc eff corrections
  * `cd ~/templateFit`
  * `rt HistsForCombine.C`
- Make the template fit, without acc eff corrections
  * `cd ~/templateFit`
  * `cmsrel CMSSW_10_3_4; cd CMSSW_10_3_4/src'
  * `cmsenv'
  * `git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit'
  * `cp ../../FitDiagnostics_modified.cc HiggsAnalysis/CombinedLimit/src/FitDiagnostics.cc; cp ../../FitDiagnostics_modified.h HiggsAnalysis/CombinedLimit/interface/FitDiagnostics.h' (to record r2 in the toys tree in FitDiagnostics)
  * `cd HiggsAnalysis/CombinedLimit/test'
  * `ln -s /home/llr/cms/falmagne/Bc/templateFit/combineFit/runCombine.sh runCombine.sh'
  * `source runCombine.sh "pp"'
  * `source runCombine.sh "PbPb"'
- Plot postfit trimuon mass shapes. The two string options correspond different fit strategies (metafit systematic variations)
  *`rt -b "plotCombineOutput.C(true,false,\"\",\"\")"'

## Yields and errors
- Correct yields with event-by-event acc eff corrections
  * `cd ~/AccEffCorr/'
  * `rt -b runMetafitSysts.C'
- Draw corrected yields for various acc eff methods, and get a systematic from it
  * `rt "Draw_corrYields.C(true)"'
  * `rt "Draw_corrYields.C(false)"'
- Draw all meta-fit variations and extract a systematic
  * `rt Draw_metafitSyst.C'
- Compute correlation factors for RAA from the ones for pp and PbPb
  * `rt -b "RAAfitErrCorr.C(true)"'
- Draw R_PbPb:
  * `cd ~/RAA/`
  * `rt Draw_XSandRAA.C'

## Automatic latex tables
- N-1 efficiencies (preselection section 3.4)
  * `cd ~/Bc/efficiency/'
  * `python Nmin1efficiencies_texTable.py'
- N-1 efficiencies (preselection section 3.4)
  * `cd ~/Bc/RAA/'
  * `python XSRAA_texTable.py'