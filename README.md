* `rt` stands for `root -l`
* Everything is shown for pp, but has to be reproduced with PbPb by changing the first option `ispp=true` (expect for the scripts gathering pp and PbPb info)
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

## Preliminaries
- Determine the BDT binning (written in a header file), for each strategy and kinematic bin. The first time, when the header file does not exist yet, have to run it with option `firstTime=true`, before the second addJpsiChoiceWeight.C running 
  * `cd ~/BDT/`
  * `rt DetermineBDTcuts.C`

## Acceptance and efficiency
- Make the acceptance and efficiency maps. Only needs the pre-selected samples, and the BDT binning. Need to run twice, with or without BDTuncorrelatedFromM (changes BDT binning, hence BDT efficiency)
  * `cd ~/acceptance/`
  * `rt BuildAcceptanceMap.C`
  * `cd ~/efficiency/`
  * `rt BuildEffMap.C`

## Make the templates histos for combine. Includes histos with acc eff corrections
  *`cd ~/templateFit`
  *`rt HistsForCombine.C`
