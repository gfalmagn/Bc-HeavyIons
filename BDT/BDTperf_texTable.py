import re #regular expressions
import ROOT

secondStep = True #consider performance of the second-step BDT training

#Get all info from BDT training
fpp = open("BDT_trainingOutput"+("_2ndStep" if secondStep else "")+"_pp.txt","r")
tpp = fpp.read()
rocpp = [[float(re.findall(" \d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin1_1stHalf: \d+\.\d+",tpp)[0]))[0])  ,
          float(re.findall(" \d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin1_2ndHalf: \d+\.\d+",tpp)[0]))[0])  ],
         [float(re.findall(" \d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin2_1stHalf: \d+\.\d+",tpp)[0]))[0])  ,
          float(re.findall(" \d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin2_2ndHalf: \d+\.\d+",tpp)[0]))[0])  ]]

effTest_pp = [[ re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin1_1stHalf: \d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)",tpp)[0]))   ,
                re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin1_2ndHalf: \d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)",tpp)[0]))   ],
              [ re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin2_1stHalf: \d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)",tpp)[0]))   ,
                re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin2_2ndHalf: \d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)",tpp)[0]))   ]]

fPbPb = open("BDT_trainingOutput"+("_2ndStep" if secondStep else "")+"_PbPb.txt","r")
tPbPb = fPbPb.read()
rocPbPb = [[float(re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_PbPb_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin1_1stHalf: \d+\.\d+",tPbPb)[0]))[0])  ,
            float(re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_PbPb_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin1_2ndHalf: \d+\.\d+",tPbPb)[0]))[0])  ],
           [float(re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_PbPb_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin2_1stHalf: \d+\.\d+",tPbPb)[0]))[0])  ,
            float(re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_PbPb_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin2_2ndHalf: \d+\.\d+",tPbPb)[0]))[0])  ]]

effTest_PbPb = [[ re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_PbPb_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin1_1stHalf: \d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)",tPbPb)[0]))   ,
                  re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_PbPb_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin1_2ndHalf: \d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)",tPbPb)[0]))   ],
                [ re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_PbPb_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin2_1stHalf: \d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)",tPbPb)[0]))   ,
                  re.findall("\d+\.\d+" , (re.findall("dataset\s+BDTfinerShallow_PbPb_withJpsiMC"+("_2ndStep" if secondStep else "")+"_kinBin2_2ndHalf: \d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)\s+\d+\.\d+ \(\d+\.\d+\)",tPbPb)[0]))   ]]

roc = {'pp':rocpp,'PbPb':rocPbPb}
effTest = {'pp':effTest_pp,'PbPb':effTest_PbPb}

#Grab ptmin, ptmax
fcuts = open("../helpers/Cuts.h", 'r')
tcuts = fcuts.read()
fcuts.close()
ptMin = [ float(re.findall(" \d+\.\d*, " , (re.findall("_BcPtmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_BcPtmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])] 
ptMax = [ float(re.findall(" \d+\.\d*, " , (re.findall("_BcPtmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_BcPtmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])] 

#Final results table
fo = open("TrainingPerfBDTTable.tex", "w")
col = ['pp','PbPb']
precis = [2,3,1,3,2,3]

for co in col:
    for b in range(0,2):
        fo.write("\\multirow{3}{*}{\\"+co+", ${:.0f}<\\pt<{:.0f}\\GeV$}}\n".format(ptMin[b], ptMax[b]))
        fo.write("& train half A, test half B & {:.{prec}f} & {:.{prec2}f} ({:.{prec2}f}) & {:.{prec2}f} ({:.{prec2}f})\\\\\n".format(roc[co][b][0] , float(effTest[co][b][0][0]) , float(effTest[co][b][0][1]) , float(effTest[co][b][0][2]) , float(effTest[co][b][0][3]) , prec=3, prec2=2))
        fo.write("& train B, test A & {:.{prec}f} & {:.{prec2}f} ({:.{prec2}f}) & {:.{prec2}f} ({:.{prec2}f})\\\\\n".format(roc[co][b][1] , float(effTest[co][b][1][0]) , float(effTest[co][b][1][1]) , float(effTest[co][b][1][2]) , float(effTest[co][b][1][3]) , prec=3, prec2=2))
        fo.write("& average & \\textbf{{{:.{prec}f}}} & \\textbf{{{:.{prec2}f} ({:.{prec2}f})}} & \\textbf{{{:.{prec2}f} ({:.{prec2}f})}}\\\\\n".format((roc[co][b][0]+roc[co][b][1])/2 , (float(effTest[co][b][0][0])+float(effTest[co][b][1][0]))/2 , (float(effTest[co][b][0][1])+float(effTest[co][b][1][1]))/2 , (float(effTest[co][b][0][2])+float(effTest[co][b][1][2]))/2 , (float(effTest[co][b][0][3])+float(effTest[co][b][1][3]))/2 , prec=3, prec2=2))
        if(b==1 and co=='pp'):
            fo.write("\\hline\\hline\n")
        else:
            fo.write("\\hline\n")

fo.close()                                                                                                                                                                                                                                 


