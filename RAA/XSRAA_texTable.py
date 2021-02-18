import re #regular expressions
import ROOT

#Get all info from central file
f = ROOT.TFile.Open("../AccEffCorr/corrected_yields_3rdStep.root","READ")
rhoAE_RAA = f.Get("AcceffSyst_Correlation_RAA")
XSRAA = f.Get("XSRAA_final")
metafitRelErr = f.Get("MetaFit_RelErr")
rhofit_pp = f.Get("r1r2Correlation_pp")
rhofit_PbPb = f.Get("r1r2Correlation_PbPb")
rhofit_RAA = f.Get("r1r2Correlation_RAA")
rhotot = f.Get("LinearizedCorrelationMatrix_total")
f.Close()                                                                                                                                                                                                                                 

f5 = ROOT.TFile.Open("../AccEffCorr/corrected_yields_2ndStep.root","READ")
rsig_pp = f5.Get("rsig_pp")#[pt bin +1]
rsig_PbPb = f5.Get("rsig_PbPb")
nsig_pp = f5.Get("nsig_pp")#[pt bin +1]
nsig_PbPb = f5.Get("nsig_PbPb")
rrelerr_pp = f5.Get("rsig_relerr_pp") #[pt bin+1][sym err, lo err, hi err]
rrelerr_PbPb = f5.Get("rsig_relerr_PbPb")
rhometaf_pp = f5.Get("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixpp")
rhometaf_PbPb = f5.Get("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixPbPb")
rhometaf_RAA = f5.Get("RAA_MetafitSyst_LinearizedCorrelationMatrix")
f5.Close()                                                                                                                                                                                                                                 

#Get info from two-steps AccEff
f2 = ROOT.TFile.Open("../twoSteps/AccEffFrom2ndStepToys.root","READ")
AEcorr_pp = f2.Get("InvAccEffFromCorrMC_LinearisedCorrelationFactor_pp_2ndStep")
AEcorr_PbPb = f2.Get("InvAccEffFromCorrMC_LinearisedCorrelationFactor_PbPb_2ndStep")
InvAccEff_pp = f2.Get("InvAccEffFromCorrMC_withSystErr_pp_2ndStep")
InvAccEff_PbPb = f2.Get("InvAccEffFromCorrMC_withSystErr_PbPb_2ndStep")
f2.Close()

#Final results table
fcuts = open("../helpers/Cuts.h", 'r')
tcuts = fcuts.read()
fcuts.close()
ptMin = [ float(re.findall(" \d+\.\d*, " , (re.findall("_BcPtmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_BcPtmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])] 
ptMax = [ float(re.findall(" \d+\.\d*, " , (re.findall("_BcPtmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_BcPtmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])] 

fo = open("ResultsTable.tex", "w")
col = ['pp','PbPb','RAA']
precis = [2,3,1,3,2,3]

for i in range(0,3):
    if(i!=0):
        fo.write("\\hline\n")
    fo.write("\\multirow{2}{*}{\\"+col[i]+"}}& ${:.0f}<\\pt<{:.0f}\\GeV$& ${:.{prec}f}^{{+{:.{prec}f}}}_{{-{:.{prec}f}}}$ &\\multirow{{2}}{{*}}{{{:.2f}}}\\\\\n".format(ptMin[0], ptMax[0], XSRAA[i][0][0], XSRAA[i][2][0], XSRAA[i][1][0], rhotot[i][0], prec=precis[2*i] ))
    fo.write(" & ${:.0f}<\\pt<{:.0f}\\GeV$& ${:.{prec}f}^{{+{:.{prec}f}}}_{{-{:.{prec}f}}}$ &\\\\\n".format(ptMin[1], ptMax[1], XSRAA[i][0][1], XSRAA[i][2][1], XSRAA[i][1][1], prec=precis[2*i+1] ))
    
fo.close()                                                                                                                                                                                                                                 

#Correlation table
fdef = open("../helpers/Definitions.h", 'r')
tdef = fdef.read()
rhoBcTau = float(re.findall("\d+\.\d*" , (re.findall("float _corr_BcTauSyst = \d+\.\d*;",tdef)[0]))[0])

fo2 = open("CorrelationsTable.tex", "w")
fo2.write("fit& {:.2f}&{:.2f}&{:.2f}\\\\\n".format(rhofit_pp[0],rhofit_PbPb[0],rhofit_RAA[0]))
fo2.write("\\hline\n")
fo2.write("fit method & {:.2f}&{:.2f}&{:.2f}\\\\\n".format(rhometaf_pp[0],rhometaf_PbPb[0],rhometaf_RAA[0]))
fo2.write("\\hline\n")
fo2.write("acc. and eff. correction & {:.2f}&{:.2f}&{:.2f}\\\\\n".format(AEcorr_pp[0],AEcorr_PbPb[0],rhoAE_RAA[0]))
fo2.write("\\hline\n")
fo2.write("$\\PBc\\to \\PJGy \\,\\PGt \\, \\PAGnGt$ decay &{:.1f}&{:.1f}&{:.1f}\\\\\n".format(rhoBcTau,rhoBcTau,rhoBcTau))
fo2.write("\\hline\n")
fo2.write("total&{:.2f}&{:.2f}&{:.2f}".format(rhotot[0][0],rhotot[1][0],rhotot[2][0]))
fo2.close()                                                                                                                                                                                                                                 

#Metafit errors table
fo3 = open("MetaFitErrorTable.tex", "w")
fo3.write("${:.0f}<\\pt<{:.0f}\\GeV$ ($b=1$)&{:.1f}\\%&{:.1f}\\%&{:.1f}\\%\\\\\n".format(ptMin[0], ptMax[0], 100*metafitRelErr[0][0][0], 100*metafitRelErr[1][0][0], 100*metafitRelErr[2][0][0]) )
fo3.write("\\hline\n")
fo3.write("${:.0f}<\\pt<{:.0f}\\GeV$ ($b=2$)&{:.1f}\\%&{:.1f}\\%&{:.1f}\\%\\\\\n".format(ptMin[1], ptMax[1], 100*metafitRelErr[0][0][1], 100*metafitRelErr[1][0][1], 100*metafitRelErr[2][0][1]) )
fo3.write("\\hline\n")
fo3.write("correlation $\\rho_{{b=1,2}}$&{:.2f} & {:.2f} & {:.2f}\\\\\n".format(rhometaf_pp[0],rhometaf_PbPb[0],rhometaf_RAA[0]) )
fo3.close()                                                                                                                                                                                                                                 

#Fit POI table
fo4 = open("FitPOITable.tex", "w")
fo4.write("\\pp& $r_1={:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$, $N={:.0f}\\pm{:.0f}$& $r_2={:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$, $N={:.0f}\\pm{:.0f}$ & ${:.2f}$\\\\\n".format(rsig_pp[1], rrelerr_pp[1][2]*rsig_pp[1], rrelerr_pp[1][1]*rsig_pp[1], 
                                                                                                                                                       nsig_pp[1], rrelerr_pp[1][0]*nsig_pp[1],
                                                                                                                                                       rsig_pp[2], rrelerr_pp[2][2]*rsig_pp[2], rrelerr_pp[2][1]*rsig_pp[2], 
                                                                                                                                                       nsig_pp[2], rrelerr_pp[2][0]*nsig_pp[2],
                                                                                                                                                       rhofit_pp[0]) )
fo4.write("\\PbPb& $r_1={:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$, $N={:.0f}\\pm{:.0f}$& $r_2={:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$, $N={:.0f}\\pm{:.0f}$ & ${:.2f}$\\\\\n".format(rsig_PbPb[1], rrelerr_PbPb[1][2]*rsig_PbPb[1], rrelerr_PbPb[1][1]*rsig_PbPb[1], 
                                                                                                                                                       nsig_PbPb[1], rrelerr_PbPb[1][0]*nsig_PbPb[1],
                                                                                                                                                       rsig_PbPb[2], rrelerr_PbPb[2][2]*rsig_PbPb[2], rrelerr_PbPb[2][1]*rsig_PbPb[2], 
                                                                                                                                                       nsig_PbPb[2], rrelerr_PbPb[2][0]*nsig_PbPb[2],
                                                                                                                                                       rhofit_PbPb[0]) )
fo4.close()                                                                                                                                                       

#Acc x Eff table 
facc = ROOT.TFile.Open("../acceptance/acceptanceMap.root","READ")
acc_oneB_step1 = facc.Get("acceptance_oneBinned"); #vector<vector<float>>
acc_oneB_step2 = facc.Get("acceptance_oneBinned_2ndStep");
acc_oneB_step3 = facc.Get("acceptance_oneBinned_3rdStep");

feff = ROOT.TFile.Open("../efficiency/AcceptanceEfficiencyMap.root","READ")
eff_oneB_step1_pp = feff.Get("efficiency_oneBinned_pp");
eff_oneB_step1_PbPb = feff.Get("efficiency_oneBinned_PbPb");
eff_oneB_step2_pp = feff.Get("efficiency_oneBinned_pp_2ndStep");
eff_oneB_step2_PbPb = feff.Get("efficiency_oneBinned_PbPb_2ndStep");
eff_oneB_step3_pp = feff.Get("efficiency_oneBinned_pp_3rdStep");
eff_oneB_step3_PbPb = feff.Get("efficiency_oneBinned_PbPb_3rdStep");

fo5 = open("AccEffTable.tex", "w")
fo5.write("\\multirow{{2}}{{*}}{{${:.0f}<\\pt<{:.0f}\\GeV$}}& \\pp& {:.3f}& {:.3f}& {:.4f}& {:.4f}\\\\\n".format(ptMin[0], ptMax[0], acc_oneB_step3[0][1], eff_oneB_step3_pp[1][0], 1/InvAccEff_pp[0][0] , acc_oneB_step1[0][1]*eff_oneB_step1_pp[1][0] ) )
fo5.write("&\\PbPb &{:.3f} &{:.3f}& {:.4f}& {:.4f}\\\\\n".format(acc_oneB_step3[1][1], eff_oneB_step3_PbPb[1][0], 1/InvAccEff_PbPb[0][0], acc_oneB_step1[1][1]*eff_oneB_step1_PbPb[1][0] ) )
fo5.write("\\hline\n")
fo5.write("\\multirow{{2}}{{*}}{{${:.0f}<\\pt<{:.0f}\\GeV$}}& \\pp& {:.3f}& {:.3f}& {:.4f}& {:.4f}\\\\\n".format(ptMin[1], ptMax[1], acc_oneB_step3[0][2], eff_oneB_step3_pp[2][0], 1/InvAccEff_pp[1][0], acc_oneB_step1[0][2]*eff_oneB_step1_pp[2][0]) )
fo5.write("&\\PbPb &{:.3f} &{:.3f}& {:.3f}& {:.3f}\\\\\n".format(acc_oneB_step3[1][2], eff_oneB_step3_PbPb[2][0], 1/InvAccEff_PbPb[1][0], acc_oneB_step1[1][2]*eff_oneB_step1_PbPb[2][0] ) )
fo5.close()


