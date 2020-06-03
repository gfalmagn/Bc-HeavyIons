'rt' stands for 'root -l'

************
For pp, without tracker muons:

    *** Do pre-selection
rt "MakeInputTrees(true,false)"
cp BDT_InputTree_pp.root BDT_InputTree_pp_copy.root
   *** Add a weight for the candidates with ambiguous Jpsi dimuon choice (2 OS pairs)
rt "addJpsiChoiceWeight.C(true,false)"
cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep1.root
   **** Train BDT
rt "ClassifierSigBkg.C(true)"
rt "addBDTvariable.C(true)"
cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep2.root
   *** Add the Jpsi choice weight, with probabilities that depend on the BDT value (purer dimuon mass distro when looking at high BDT)
rt "addJpsiChoiceWeight.C(true,true)"
cp BDT_InputTree_pp.root BDT_InputTree_pp_copystep3.root
