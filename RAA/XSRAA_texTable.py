import re #regular expressions
import ROOT

#Get all info from central file
f = ROOT.TFile.Open("../AccEffCorr/corrected_yields_3rdStep.root","READ")
rhoAE_RAA = f.Get("AcceffSyst_Correlation_RAA")
rhoAE_RAAcent = f.Get("AcceffSyst_Correlation_RAA_centralityBins")
fullYcorr_RAA = f.Get("FullErrorCorrelation_RAA")
XSRAA = f.Get("XSRAA_final")
XSRAAcent = f.Get("XSRAA_final_centralityBins")
metafitRelErr = f.Get("MetaFit_RelErr")
metafitRelErr_cent = f.Get("MetaFit_RelErr_centralityBins")
rhofit_RAA = f.Get("r1r2Correlation_RAA")
rhofit_RAAcent = f.Get("r1r2Correlation_RAA_centralityBins")
rhotot = f.Get("LinearizedCorrelationMatrix_total")
rhotot_cent = f.Get("LinearizedCorrelationMatrix_total_centralityBins")
f.Close()                                                                                                                                                                                                                                 

f5 = ROOT.TFile.Open("../AccEffCorr/corrected_yields_2ndStep.root","READ")
rsig_pp = f5.Get("rsig_pp")#[pt bin +1]
rsig_PbPb = f5.Get("rsig_PbPb")
rsig_PbPbcent = f5.Get("rsig_centralityDep_PbPb")
nsig_pp = f5.Get("nsig_pp")#[pt bin +1]
nsig_PbPb = f5.Get("nsig_PbPb")
nsig_PbPbcent = f5.Get("nsig_centralityDep_PbPb")
rrelerr_pp = f5.Get("rsig_relerr_pp") #[pt bin+1][sym err, lo err, hi err]
rrelerr_PbPb = f5.Get("rsig_relerr_PbPb")
rrelerr_PbPbcent = f5.Get("rsig_relerr_centralityDep_PbPb")
rhofit_pp = f5.Get("r1r2Correlation_pp")
rhofit_PbPb = f5.Get("r1r2Correlation_PbPb")
rhofit_PbPbcent = f5.Get("r1r2Correlation_centralityDep_PbPb")
rhometaf_pp = f5.Get("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixpp")
rhometaf_PbPb = f5.Get("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrixPbPb")
rhometaf_RAA = f5.Get("RAA_MetafitSyst_LinearizedCorrelationMatrix")
rhometaf_PbPbcent = f5.Get("CorrectedYields_MetafitSyst_LinearizedCorrelationMatrix_CentralityDep")
rhometaf_RAAcent = f5.Get("RAA_MetafitSyst_LinearizedCorrelationMatrix_CentralityDep")
f5.Close()                                                                                                                                                                                                                                 

f6 = ROOT.TFile.Open("../AccEffCorr/corrected_yields.root","READ")
rsig1stStep_pp = f6.Get("rsig_pp")#[pt bin +1]
rsig1stStep_PbPb = f6.Get("rsig_PbPb")
rsig1stStep_PbPbcent = f6.Get("rsig_centralityDep_PbPb")
f6.Close()                                                                                                                                                                                                                                 

#Get info from two-steps AccEff
f2 = ROOT.TFile.Open("../twoSteps/AccEffFrom2ndStepToys.root","READ")
AEcorr_pp = f2.Get("InvAccEffFromCorrMC_LinearisedCorrelationFactor_pp_2ndStep")
AEcorr_PbPb = f2.Get("InvAccEffFromCorrMC_LinearisedCorrelationFactor_PbPb_2ndStep")
AEcorr_PbPbcent = f2.Get("InvAccEffFromCorrMC_LinearisedCorrelationFactor_PbPb_inCentBins_2ndStep")
fullYcorr_pp = f2.Get("CorrYieldsFromCorrMC_LinearisedCorrelationFactor_pp_2ndStep")
fullYcorr_PbPb = f2.Get("CorrYieldsFromCorrMC_LinearisedCorrelationFactor_PbPb_2ndStep")
InvAccEff_pp = f2.Get("InvAccEffFromCorrMC_withSystErr_pp_2ndStep")
InvAccEff_PbPb = f2.Get("InvAccEffFromCorrMC_withSystErr_PbPb_2ndStep")
InvAccEff_PbPbcent = f2.Get("InvAccEffFromCorrMC_withSystErr_PbPb_inCentBins_2ndStep")
f2.Close()

#Grab pt min max
fcuts = open("../helpers/Cuts.h", 'r')
tcuts = fcuts.read()
fcuts.close()
ptMin = [ float(re.findall(" \d+\.\d*, " , (re.findall("_BcPtmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_BcPtmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])] 
ptMax = [ float(re.findall(" \d+\.\d*, " , (re.findall("_BcPtmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_BcPtmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])] 
centMin = [ float(re.findall(" \d+\.\d*, " , (re.findall("_Centmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_Centmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])] 
centMax = [ float(re.findall(" \d+\.\d*, " , (re.findall("_Centmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_Centmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])] 

#Grab TnP error correlation
fdef = open("../helpers/Definitions.h", 'r')
tdef = fdef.read()
fdef.close()
corrTnPerr = float(re.findall("\d+\.\d*" , (re.findall("_corrTnPerr = \d+\.\d*;",tdef)[0]))[0])

#Final results table
fo = open("ResultsTable.tex", "w")
precis = [1,3,1,1,3,0,0,0,2,2,2,2,2]
ptrange = '${:.0f}<\\pt<{:.0f}\\GeV$'
centrange = 'centrality ${:.0f}\\textnormal{{-}}{:.0f}$\\%'
quantities = ['\\multirow{{2}}{{*}}{{$BF \\times \\frac{{\\textrm{{d}}\\sigma_{{\\pp}}}}{{\\textrm{{d}}p_{{T}}\\textrm{{d}}y}}$ [pb/\\GeV]}}',
              '$BF \\times \\sigma_{{\\pp}}$ [pb]',
              '\\multirow{{2}}{{*}}{{$BF \\times \\frac{{1}}{{N_{{MB}}^{{\\mathcal{{C}}}} T_{{\\PbPb}}^{{\\mathcal{{C}}}}}} \\frac{{\\textrm{{d}}N_{{\\PbPb}}^{{corr}}}}{{\\textrm{{d}}p_{{T}} \\textrm{{d}}y}}$ [pb/\\GeV]}}',
              '\\multirow{{3}}{{*}}{{$BF \\times \\frac{{1}}{{N_{{MB}}^{{\\mathcal{{C}}}} T_{{\\PbPb}}^{{\\mathcal{{C}}}}}} N_{{\\PbPb}}^{{corr}}(\\mathcal{{C}})$ [pb/$\\mathcal{{C}}$]}}',
              '\\multirow{{2}}{{*}}{{\\RAA(\\pt)}}', 
              '\\multirow{{2}}{{*}}{{$\\RAA(\\mathcal{{C}})$}}',
              '\\RAA','']
quantidx = [0,7,1,2,7,3,7,7,4,7,5,7,6]
systemstr = ["\\multirow{{3}}{{*}}{{\\pp}}","","","\\multirow{{5}}{{*}}{{\\PbPb}}","","","","","\\multirow{{5}}{{*}}{{\\Large$\\frac{{\\PbPb}}{{\\pp}}$}}","","","",""]
binstr = [ptrange,ptrange,"integrated",ptrange,ptrange,centrange,centrange,"integrated",ptrange,ptrange,centrange,centrange,"integrated"]
corrstr = ["\\multirow{{2}}{{*}}{{{:.2f}}}",'','$-$',"\\multirow{{2}}{{*}}{{{:.2f}}}",'',"\\multirow{{2}}{{*}}{{{:.2f}}}",'','$-$',"\\multirow{{2}}{{*}}{{{:.2f}}}",'',"\\multirow{{2}}{{*}}{{{:.2f}}}",'','$-$']
vals    = [XSRAA[0][0][0] , XSRAA[0][0][1], XSRAAcent[0][0][1], XSRAA[1][0][0] , XSRAA[1][0][1], XSRAAcent[1][0][0] , XSRAAcent[1][0][1], XSRAAcent[0][0][0], XSRAA[2][0][0] , XSRAA[2][0][1], XSRAAcent[2][0][0] , XSRAAcent[2][0][1], XSRAAcent[3][0][0] ]
valsEHi = [XSRAA[0][2][0] , XSRAA[0][2][1], XSRAAcent[0][2][1], XSRAA[1][2][0] , XSRAA[1][2][1], XSRAAcent[1][2][0] , XSRAAcent[1][2][1], XSRAAcent[0][2][0], XSRAA[2][2][0] , XSRAA[2][2][1], XSRAAcent[2][2][0] , XSRAAcent[2][2][1], XSRAAcent[3][2][0] ]
valsELo = [XSRAA[0][1][0] , XSRAA[0][1][1], XSRAAcent[0][1][1], XSRAA[1][1][0] , XSRAA[1][1][1], XSRAAcent[1][1][0] , XSRAAcent[1][1][1], XSRAAcent[0][1][0], XSRAA[2][1][0] , XSRAA[2][1][1], XSRAAcent[2][1][0] , XSRAAcent[2][1][1], XSRAAcent[3][1][0] ]
corrval = [rhotot[0][0] ,0,0, rhotot[1][0] ,0, rhotot_cent[1][0] ,0,0, rhotot[2][0] ,0, rhotot_cent[2][0] ,0,0]
valsPrint = '${:.{prec}f}^{{+{:.{prec}f}}}_{{-{:.{prec}f}}}$'
valsPrintBf = '\\textbf{{${:.{prec}f}^{{+{:.{prec}f}}}_{{-{:.{prec}f}}}$}}'

for i in range(0,13):
    rangm = []
    rangM = []
    if(binstr[i]==ptrange):
        rangm = ptMin
        rangM = ptMax
    if(binstr[i]==centrange):
        rangm = centMin
        rangM = centMax
    valsPrin = valsPrint if (i<8) else valsPrintBf

    if(i!=0):
        if(systemstr[i]!=''):
            fo.write("\\hline\n")
        if(binstr[i]=="integrated" or (corrstr[i]!='' and corrstr[i]!='$-$')):
            fo.write("\\cline{"+('4' if (i==7) else '2')+"-5}\n")

    if(corrstr[i]!='' and corrstr[i]!='$-$'):
        fo.write((systemstr[i]+" & "+binstr[i]+" & "+quantities[quantidx[i]]+" & "+valsPrin+" & "+corrstr[i]+"\\\\\n").format(rangm[0], rangM[0], vals[i], valsEHi[i], valsELo[i], corrval[i], prec=precis[i] ))
    elif(binstr[i]=="integrated"):
        fo.write((systemstr[i]+" & "+binstr[i]+" & "+quantities[quantidx[i]]+" & "+valsPrin+" & "+corrstr[i]+"\\\\\n").format(vals[i], valsEHi[i], valsELo[i], prec=precis[i] ))
    else:
        fo.write((systemstr[i]+" & "+binstr[i]+" & "+quantities[quantidx[i]]+" & "+valsPrin+" & "+corrstr[i]+"\\\\\n").format(rangm[1], rangM[1], vals[i], valsEHi[i], valsELo[i], prec=precis[i] ))

fo.close()                                                                                                                                                                                                                                 

#Correlation table
fdef = open("../helpers/Definitions.h", 'r')
tdef = fdef.read()
rhoBcTau = float(re.findall("\d+\.\d*" , (re.findall("float _corr_BcTauSyst = \d+\.\d*;",tdef)[0]))[0])

fo2 = open("CorrelationsTable.tex", "w")
fo2.write("fit& {:.2f}&{:.2f}&{:.2f}&{:.2f}&{:.2f}\\\\\n".format(rhofit_pp[0],rhofit_PbPb[0],rhofit_RAA[0],rhofit_PbPbcent[0],rhofit_RAAcent[0]))
fo2.write("\\hline\n")
fo2.write("fit method & {:.2f}&{:.2f}&{:.2f}&{:.2f}&{:.2f}\\\\\n".format(rhometaf_pp[0],rhometaf_PbPb[0],rhometaf_RAA[0],rhometaf_PbPbcent[0],rhometaf_RAAcent[0]))
fo2.write("\\hline\n")
fo2.write("acc. and eff. correction & {:.2f}&{:.2f}&{:.2f}&{:.2f}&{:.2f}\\\\\n".format(AEcorr_pp[0],AEcorr_PbPb[0],rhoAE_RAA[0],AEcorr_PbPbcent[0],rhoAE_RAAcent[0]))
fo2.write("\\hline\n")
fo2.write("$\\PBc\\to \\PJGy \\,\\PGt \\, \\PAGnGt$ decay + lumi. &{:.1f}&{:.1f}&{:.1f}&{:.1f}&{:.1f}\\\\\n".format(rhoBcTau,rhoBcTau,rhoBcTau,rhoBcTau,rhoBcTau))
fo2.write("\\hline\n")
fo2.write("tag-and-probe scale factors &{:.1f}&{:.1f}&{:.1f}&{:.1f}&{:.1f}\\\\\n".format(corrTnPerr,corrTnPerr,corrTnPerr,corrTnPerr,corrTnPerr))
fo2.write("\\hline\n")
fo2.write("full (w/o lumi.+$\\PBc\\to \\PJGy \\,\\PGt$)&{:.2f}&{:.2f}&{:.2f}&$-$&$-$\\\\\n".format(fullYcorr_pp[0],fullYcorr_PbPb[0],fullYcorr_RAA[0]))
fo2.write("\\hline\n")
fo2.write("total&\\textbf{{{:.2f}}}&\\textbf{{{:.2f}}}&\\textbf{{{:.2f}}}&\\textbf{{{:.2f}}}&\\textbf{{{:.2f}}}\\\\\n".format(rhotot[0][0],rhotot[1][0],rhotot[2][0],rhotot_cent[1][0],rhotot_cent[2][0]))
fo2.close()                                                                                                                                                                                                                                 

#Metafit errors table
fo3 = open("MetaFitErrorTable.tex", "w")
fo3.write((ptrange+" &{:.1f}\\%&{:.1f}\\%&{:.1f}\\%\\\\\n").format(ptMin[0], ptMax[0], 100*metafitRelErr[0][0][0], 100*metafitRelErr[1][0][0], 100*metafitRelErr[2][0][0]) )
fo3.write((ptrange+" &{:.1f}\\%&{:.1f}\\%&{:.1f}\\%\\\\\n").format(ptMin[1], ptMax[1], 100*metafitRelErr[0][0][1], 100*metafitRelErr[1][0][1], 100*metafitRelErr[2][0][1]) )
fo3.write("correlation $\\rho_{{1,2}}$&{:.2f} & {:.2f} & {:.2f}\\\\\n".format(rhometaf_pp[0],rhometaf_PbPb[0],rhometaf_RAA[0]) )
fo3.write("\\hline\n")
fo3.write((centrange+" &$-$&{:.1f}\\%&{:.1f}\\%\\\\\n").format(centMin[0], centMax[0], 100*metafitRelErr_cent[1][0][0], 100*metafitRelErr_cent[2][0][0]) )
fo3.write((centrange+" &$-$&{:.1f}\\%&{:.1f}\\%\\\\\n").format(centMin[1], centMax[1], 100*metafitRelErr_cent[1][0][1], 100*metafitRelErr_cent[2][0][1]) )
fo3.write("correlation $\\rho_{{1,2}}$& $-$ & {:.2f} & {:.2f}\\\\\n".format(rhometaf_PbPbcent[0],rhometaf_RAAcent[0]) )
fo3.write("\\hline\n")
fo3.write("integrated &{:.1f}\\%&{:.1f}\\%&{:.1f}\\%\\\\\n".format(100*metafitRelErr_cent[0][0][1], 100*metafitRelErr_cent[0][0][0], 100*metafitRelErr_cent[3][0][0]) )
fo3.close()                                                                                                                                                                                                                                 

#Fit POI table
fo4 = open("FitPOITable.tex", "w")
fo4.write(("\\multirow{{3}}{{*}}{{\\pp}}& "+ptrange+" & ${:.3f}$ &  ${:.3f}$ & ${:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$ & $\\textbf{{{:.0f}}}^{{+{:.0f}}}_{{-{:.0f}}}$ & \\multirow{{2}}{{*}}{{${:.2f}$}}\\\\\n").format(
    ptMin[0], ptMax[0], 
    rsig1stStep_pp[1], rsig_pp[1], 
    rsig1stStep_pp[1]*rsig_pp[1] , rrelerr_pp[1][2]*rsig1stStep_pp[1]*rsig_pp[1], rrelerr_pp[1][1]*rsig1stStep_pp[1]*rsig_pp[1],
    nsig_pp[1], rrelerr_pp[1][2]*nsig_pp[1], rrelerr_pp[1][1]*nsig_pp[1],
    rhofit_pp[0]) )
fo4.write(("& "+ptrange+" & ${:.3f}$ &  ${:.3f}$ & ${:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$ & $\\textbf{{{:.0f}}}^{{+{:.0f}}}_{{-{:.0f}}}$ & \\\\\n").format(
    ptMin[1], ptMax[1], 
    rsig1stStep_pp[2], rsig_pp[2], 
    rsig1stStep_pp[2]*rsig_pp[2] , rrelerr_pp[2][2]*rsig1stStep_pp[2]*rsig_pp[2], rrelerr_pp[2][1]*rsig1stStep_pp[2]*rsig_pp[2],
    nsig_pp[2], rrelerr_pp[2][2]*nsig_pp[2], rrelerr_pp[2][1]*nsig_pp[2]) )
fo4.write("\\cline{2-7}\n")
fo4.write("& integrated & ${:.3f}$ &  ${:.3f}$ & ${:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$ & $\\textbf{{{:.0f}}}^{{+{:.0f}}}_{{-{:.0f}}}$ & $-$\\\\\n".format(
    rsig1stStep_pp[0], rsig_pp[0], 
    rsig1stStep_pp[0]*rsig_pp[0] , rrelerr_pp[0][2]*rsig1stStep_pp[0]*rsig_pp[0], rrelerr_pp[0][1]*rsig1stStep_pp[0]*rsig_pp[0],
    nsig_pp[0], rrelerr_pp[0][2]*nsig_pp[0], rrelerr_pp[0][1]*nsig_pp[0]) )
fo4.write("\\hline\\hline\n")

fo4.write(("\\multirow{{5}}{{*}}{{\\PbPb}}& "+ptrange+" & ${:.3f}$ &  ${:.3f}$ & ${:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$ & $\\textbf{{{:.0f}}}^{{+{:.0f}}}_{{-{:.0f}}}$ & \\multirow{{2}}{{*}}{{${:.2f}$}}\\\\\n").format(
    ptMin[0], ptMax[0], 
    rsig1stStep_PbPb[1], rsig_PbPb[1], 
    rsig1stStep_PbPb[1]*rsig_PbPb[1] , rrelerr_PbPb[1][2]*rsig1stStep_PbPb[1]*rsig_PbPb[1], rrelerr_PbPb[1][1]*rsig1stStep_PbPb[1]*rsig_PbPb[1],
    nsig_PbPb[1], rrelerr_PbPb[1][2]*nsig_PbPb[1], rrelerr_PbPb[1][1]*nsig_PbPb[1],
    rhofit_PbPb[0]) )
fo4.write(("& "+ptrange+" & ${:.3f}$ &  ${:.3f}$ & ${:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$ & $\\textbf{{{:.0f}}}^{{+{:.0f}}}_{{-{:.0f}}}$ & \\\\\n").format(
    ptMin[1], ptMax[1], 
    rsig1stStep_PbPb[2], rsig_PbPb[2], 
    rsig1stStep_PbPb[2]*rsig_PbPb[2] , rrelerr_PbPb[2][2]*rsig1stStep_PbPb[2]*rsig_PbPb[2], rrelerr_PbPb[2][1]*rsig1stStep_PbPb[2]*rsig_PbPb[2],
    nsig_PbPb[2], rrelerr_PbPb[2][2]*nsig_PbPb[2], rrelerr_PbPb[2][1]*nsig_PbPb[2]) )
fo4.write("\\cline{2-7}\n")

fo4.write(("& "+centrange+" & ${:.3f}$ &  ${:.3f}$ & ${:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$ & $\\textbf{{{:.0f}}}^{{+{:.0f}}}_{{-{:.0f}}}$ & \\multirow{{2}}{{*}}{{${:.2f}$}}\\\\\n").format(
    centMin[0], centMax[0], 
    rsig1stStep_PbPbcent[1], rsig_PbPbcent[1], 
    rsig1stStep_PbPbcent[1]*rsig_PbPbcent[1] , rrelerr_PbPbcent[1][2]*rsig1stStep_PbPbcent[1]*rsig_PbPbcent[1], rrelerr_PbPbcent[1][1]*rsig1stStep_PbPbcent[1]*rsig_PbPbcent[1],
    nsig_PbPbcent[1], rrelerr_PbPbcent[1][2]*nsig_PbPbcent[1], rrelerr_PbPbcent[1][1]*nsig_PbPbcent[1],
    rhofit_PbPbcent[0]) )
fo4.write(("& "+centrange+" & ${:.3f}$ &  ${:.3f}$ & ${:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$ & $\\textbf{{{:.0f}}}^{{+{:.0f}}}_{{-{:.0f}}}$ & \\\\\n").format(
    centMin[1], centMax[1], 
    rsig1stStep_PbPbcent[2], rsig_PbPbcent[2], 
    rsig1stStep_PbPbcent[2]*rsig_PbPbcent[2] , rrelerr_PbPbcent[2][2]*rsig1stStep_PbPbcent[2]*rsig_PbPbcent[2], rrelerr_PbPbcent[2][1]*rsig1stStep_PbPbcent[2]*rsig_PbPbcent[2],
    nsig_PbPbcent[2], rrelerr_PbPbcent[2][2]*nsig_PbPbcent[2], rrelerr_PbPbcent[2][1]*nsig_PbPbcent[2]) )
fo4.write("\\cline{2-7}\n")
fo4.write("& integrated & ${:.3f}$ &  ${:.3f}$ & ${:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$ & $\\textbf{{{:.0f}}}^{{+{:.0f}}}_{{-{:.0f}}}$ & $-$\\\\\n".format(
    rsig1stStep_PbPb[0], rsig_PbPb[0], 
    rsig1stStep_PbPb[0]*rsig_PbPb[0] , rrelerr_PbPb[0][2]*rsig1stStep_PbPb[0]*rsig_PbPb[0], rrelerr_PbPb[0][1]*rsig1stStep_PbPb[0]*rsig_PbPb[0],
    nsig_PbPb[0], rrelerr_PbPb[0][2]*nsig_PbPb[0], rrelerr_PbPb[0][1]*nsig_PbPb[0]) )

fo4.write("\\hline\n")
fo4.close()                                                                                                                                                       

#Acc x Eff table 
facc = ROOT.TFile.Open("../acceptance/acceptanceMap.root","READ")
acc_oneB_step1 = facc.Get("acceptance_oneBinned"); #vector<vector<float>>
acc_oneB_step2 = facc.Get("acceptance_oneBinned_2ndStep");
acc_oneB_step3 = facc.Get("acceptance_oneBinned_3rdStep");

feff = ROOT.TFile.Open("../efficiency/AcceptanceEfficiencyMap.root","READ")
eff_oneB_step1_pp = feff.Get("efficiency_oneBinned_pp");
eff_oneB_step1_PbPb = feff.Get("efficiency_oneBinned_PbPb");
eff_oneB_step1_PbPbcent = feff.Get("efficiency_oneBinned_centDep_PbPb");
eff_oneB_step2_pp = feff.Get("efficiency_oneBinned_pp_2ndStep");
eff_oneB_step2_PbPb = feff.Get("efficiency_oneBinned_PbPb_2ndStep");
eff_oneB_step2_PbPbcent = feff.Get("efficiency_oneBinned_centDep_PbPb_2ndStep");
eff_oneB_step3_pp = feff.Get("efficiency_oneBinned_pp_3rdStep");
eff_oneB_step3_PbPb = feff.Get("efficiency_oneBinned_PbPb_3rdStep");
eff_oneB_step3_PbPbcent = feff.Get("efficiency_oneBinned_centDep_PbPb_3rdStep");

fo5 = open("AccEffTable.tex", "w")
fo5.write(("\\multirow{{3}}{{*}}{{\\pp}}& "+ptrange+" & {:.3f}& {:.2f}& {:.3f}& {:.4f}& \\textbf{{{:.4f}}}\\\\\n").format(ptMin[0], ptMax[0], acc_oneB_step3[0][1], eff_oneB_step3_pp[1][0] , acc_oneB_step1[0][1]*eff_oneB_step1_pp[1][0], acc_oneB_step2[0][1]*eff_oneB_step2_pp[1][0], 1/InvAccEff_pp[1][0]  ) )
fo5.write(("& "+ptrange+" & {:.3f}& {:.2f}& {:.3f}& {:.3f}& \\textbf{{{:.3f}}}\\\\\n").format(ptMin[1], ptMax[1], acc_oneB_step3[0][2], eff_oneB_step3_pp[2][0], acc_oneB_step1[0][2]*eff_oneB_step1_pp[2][0], acc_oneB_step2[0][2]*eff_oneB_step2_pp[2][0], 1/InvAccEff_pp[2][0]) )
fo5.write("& integrated & {:.3f}& {:.2f}& {:.3f}& {:.4f}& \\textbf{{{:.4f}}}\\\\\n".format(acc_oneB_step3[0][0], eff_oneB_step3_pp[0][0], acc_oneB_step1[0][0]*eff_oneB_step1_pp[0][0], acc_oneB_step2[0][0]*eff_oneB_step2_pp[0][0], 1/InvAccEff_pp[0][0] ) )
fo5.write("\\hline\n")

fo5.write(("\\multirow{{5}}{{*}}{{\\PbPb}}& "+ptrange+" & {:.3f}& {:.3f}& {:.4f}& {:.4f}& \\textbf{{{:.4f}}}\\\\\n").format(ptMin[0], ptMax[0], acc_oneB_step3[1][1], eff_oneB_step3_PbPb[1][0], acc_oneB_step1[1][1]*eff_oneB_step1_PbPb[1][0], acc_oneB_step2[1][1]*eff_oneB_step2_PbPb[1][0], 1/InvAccEff_PbPb[1][0]  ) )
fo5.write(("& "+ptrange+" & {:.2f}& {:.2f}& {:.3f}& {:.3f}& \\textbf{{{:.3f}}}\\\\\n").format(ptMin[1], ptMax[1], acc_oneB_step3[1][2], eff_oneB_step3_PbPb[2][0], acc_oneB_step1[1][2]*eff_oneB_step1_PbPb[2][0], acc_oneB_step2[1][2]*eff_oneB_step2_PbPb[2][0], 1/InvAccEff_PbPb[2][0]) )
fo5.write(("& "+centrange+" & {:.3f}& {:.3f}& {:.4f}& {:.4f}& \\textbf{{{:.4f}}}\\\\\n").format(centMin[0], centMax[0], acc_oneB_step3[1][0], eff_oneB_step3_PbPbcent[1][0], acc_oneB_step1[1][0]*eff_oneB_step1_PbPbcent[1][0], acc_oneB_step2[1][0]*eff_oneB_step2_PbPbcent[1][0], 1/InvAccEff_PbPbcent[1][0]  ) )
fo5.write(("& "+centrange+" & {:.2f}& {:.2f}& {:.3f}& {:.3f}& \\textbf{{{:.3f}}}\\\\\n").format(centMin[1], centMax[1], acc_oneB_step3[1][0], eff_oneB_step3_PbPbcent[2][0], acc_oneB_step1[1][0]*eff_oneB_step1_PbPbcent[2][0], acc_oneB_step2[1][0]*eff_oneB_step2_PbPbcent[2][0], 1/InvAccEff_PbPbcent[2][0]) )
fo5.write("& integrated & {:.3f}& {:.3f}& {:.4f}& {:.4f}& \\textbf{{{:.4f}}}\\\\\n".format(acc_oneB_step3[1][0], eff_oneB_step3_PbPb[0][0], acc_oneB_step1[1][0]*eff_oneB_step1_PbPb[0][0], acc_oneB_step2[1][0]*eff_oneB_step2_PbPb[0][0], 1/InvAccEff_PbPb[0][0] ) )
fo5.write("\\hline\n")
fo5.close()

