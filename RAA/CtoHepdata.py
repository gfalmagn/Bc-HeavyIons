import re #regular expressions
import ROOT
import math

#pT dependence
fraa = open("RAA.C", 'r')
traa = fraa.read()
fraa.close()
XY = []
Xerr = []
YerrTot = []
YerrUnco = []

for i in range(0,1000):
    if (re.findall("->SetPoint\("+str(i)+", \d+\.\d*, \d+\.\d*\)",traa) == [] ):
        #print 'found no such string for point #',i
        break
    XY.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPoint\("+str(i)+", \d+\.\d*, \d+\.\d*\)",traa)[0] )[0][1:-2]),
                 float(re.findall(" \d+\.\d*\)" , re.findall("->SetPoint\("+str(i)+", \d+\.\d*, \d+\.\d*\)",traa)[0] )[0][1:-2]) ] )
    Xerr.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEX\("+str(i)+", \d+\.\d*, \d+\.\d*\)",traa)[0] )[0][1:-2]),
                 float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEX\("+str(i)+", \d+\.\d*, \d+\.\d*\)",traa)[0] )[0][1:-2]) ] )
    YerrTot.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEY\("+str(i)+", 0, \d+\.\d*, \d+\.\d*\)",traa)[0] )[0][1:-2]),
                      float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEY\("+str(i)+", 0, \d+\.\d*, \d+\.\d*\)",traa)[0] )[0][1:-2]) ] )
    YerrUnco.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEY\("+str(i)+", 1, \d+\.\d*, \d+\.\d*\)",traa)[0] )[0][1:-2]),
                       float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEY\("+str(i)+", 1, \d+\.\d*, \d+\.\d*\)",traa)[0] )[0][1:-2]) ] )
rho12 = float(re.findall(" \d+\.\d*" , re.findall("rho_{1-2} = \d+\.\d*",traa)[0])[0][1:])

#Npart dependence
fraa = open("RAA_centralityBins_vsNpart.C", 'r')
traa_ct = fraa.read()
fraa.close()
XY_ct = []
YerrTot_ct = []
YerrUnco_ct = []

for i in range(0, len(re.findall("->SetPoint\(0, \d+\.\d*, \d+\.\d*\)",traa_ct))):
    XY_ct.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPoint\(0, \d+\.\d*, \d+\.\d*\)",traa_ct)[i] )[0][1:-2]),
                    float(re.findall(" \d+\.\d*\)" , re.findall("->SetPoint\(0, \d+\.\d*, \d+\.\d*\)",traa_ct)[i] )[0][1:-2]) ] )
    YerrTot_ct.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEY\("+str(i)+", 0, \d+\.\d*, \d+\.\d*\)",traa_ct)[0] )[0][1:-2]),
                         float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEY\("+str(i)+", 0, \d+\.\d*, \d+\.\d*\)",traa_ct)[0] )[0][1:-2]) ] )
    YerrUnco_ct.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEY\("+str(i)+", 1, \d+\.\d*, \d+\.\d*\)",traa_ct)[0] )[0][1:-2]),
                          float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEY\("+str(i)+", 1, \d+\.\d*, \d+\.\d*\)",traa_ct)[0] )[0][1:-2]) ] )
rho12_ct = float(re.findall(" \d+\.\d*" , re.findall("rho_{1-2} = \d+\.\d*",traa_ct)[0])[0][1:])

#pT dependent XS
fxs = open("CrossSections.C", 'r')
txs = fxs.read()
fxs.close()
XY_xspp = []
Xerr_xspp = []
YerrTot_xspp = []
YerrFit_xspp = []
XY_xsPbPb = []
Xerr_xsPbPb = []
YerrTot_xsPbPb = []
YerrFit_xsPbPb = []
npts = len(re.findall("->SetPoint\(0, \d+\.\d*, \d+\.\d*\)",txs))/2

for i in range(0, npts):
    XY_xsPbPb.append(      [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPoint\(0, \d+\.\d*, \d+\.\d*\)",txs)[1+i] )[0][1:-2]),
                           float(re.findall(" \d+\.\d*\)" , re.findall("->SetPoint\(0, \d+\.\d*, \d+\.\d*\)",txs)[1+i] )[0][1:-2]) ] )
    Xerr_xsPbPb.append(    [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEX\(0, \d+\.\d*, \d+\.\d*\)",txs)[1+i] )[0][1:-2]),
                           float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEX\(0, \d+\.\d*, \d+\.\d*\)",txs)[1+i] )[0][1:-2]) ] )
    YerrTot_xsPbPb.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEY\(0, 0, \d+\.\d*, \d+\.\d*\)",txs)[1+i] )[0][1:-2]),
                           float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEY\(0, 0, \d+\.\d*, \d+\.\d*\)",txs)[1+i] )[0][1:-2]) ] )
    YerrFit_xsPbPb.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEY\(0, 1, \d+\.\d*, \d+\.\d*\)",txs)[1+i] )[0][1:-2]),
                           float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEY\(0, 1, \d+\.\d*, \d+\.\d*\)",txs)[1+i] )[0][1:-2]) ] )

    XY_xspp.append(      [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPoint\(0, \d+\.\d*, \d+\.\d*\)",txs)[npts+1+i] )[0][1:-2]),
                             float(re.findall(" \d+\.\d*\)" , re.findall("->SetPoint\(0, \d+\.\d*, \d+\.\d*\)",txs)[npts+1+i] )[0][1:-2]) ] )
    Xerr_xspp.append(    [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEX\(0, \d+\.\d*, \d+\.\d*\)",txs)[npts+1+i] )[0][1:-2]),
                             float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEX\(0, \d+\.\d*, \d+\.\d*\)",txs)[npts+1+i] )[0][1:-2]) ] )
    YerrTot_xspp.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEY\(0, 0, \d+\.\d*, \d+\.\d*\)",txs)[npts+1+i] )[0][1:-2]),
                             float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEY\(0, 0, \d+\.\d*, \d+\.\d*\)",txs)[npts+1+i] )[0][1:-2]) ] )
    YerrFit_xspp.append( [ float(re.findall(" \d+\.\d*, " , re.findall("->SetPointEY\(0, 1, \d+\.\d*, \d+\.\d*\)",txs)[npts+1+i] )[0][1:-2]),
                             float(re.findall(" \d+\.\d*\)" , re.findall("->SetPointEY\(0, 1, \d+\.\d*, \d+\.\d*\)",txs)[npts+1+i] )[0][1:-2]) ] )
rho12_xspp = float(re.findall(" \d+\.\d*" , re.findall("rho_{1-2}\^{pp} = \d+\.\d*",txs)[0])[0][1:])
rho12_xsPbPb = float(re.findall(" \d+\.\d*" , re.findall("rho_{1-2}\^{PbPb} = \d+\.\d*",txs)[0])[0][1:])

#Grab centrality ranges 
fcuts = open("../helpers/Cuts.h", 'r')
tcuts = fcuts.read()
fcuts.close()
centMin = [ float(re.findall(" \d+\.\d*, " , (re.findall("_Centmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_Centmin\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])]
centMax = [ float(re.findall(" \d+\.\d*, " , (re.findall("_Centmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-2]),
          float(re.findall(" \d+\.\d*\}" , (re.findall("_Centmax\{\d+\.\d*, \d+\.\d*, \d+\.\d*\}",tcuts)[0]))[0][:-1])]

#HepData results table 
fo = open("HepDataTable_RAA.txt", "w")

fo.write("# RAA(Bc) versus visible pT(trimuon), with CMS 2017 PbPb and 2018 PbPb data\n")
fo.write("# x xmin xmax ycenter +stat_u -stat_d +syst_u -syst_d\n")
fo.write("# correlation factor between bins 1 and 2")
fo.write("# 'stat' uncertainty is the part of the total uncertainty that is not bin-to-bin correlated (mostly statistical). 'syst' uncertainty is such that total uncertainty is (stat^2+syst^2) \n")
for i in range(0,len(XY)):
    fo.write( "{:.3f} {:.2f} {:.2f} {:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f}\n".format(XY[i][0], XY[i][0]-Xerr[i][0], XY[i][0]+Xerr[i][1]+1e-4,
                                                                                     XY[i][1], YerrTot[i][1], YerrTot[i][0], YerrUnco[i][1], YerrUnco[i][0], math.sqrt(YerrTot[i][1] ** 2 - YerrUnco[i][1] ** 2), math.sqrt(YerrTot[i][0] ** 2 - YerrUnco[i][0] ** 2) ))
fo.write( "{:.2f}\n".format(rho12) )

fo.write("\n# RAA(Bc) versus N_part, with CMS 2017 pp and 2018 PbPb data\n")
fo.write("# x xmin xmax ycenter +stat_u -stat_d +syst_u -syst_d\n")
fo.write("# correlation factor between bins 1 and 2")
fo.write("# 'stat' uncertainty is the part of the total uncertainty that is not bin-to-bin correlated (mostly statistical). 'syst' uncertainty is such that total uncertainty is (stat^2+syst^2) \n")
for i in range(0,len(XY_ct)):
    fo.write( "{:.3f} {:.3f} {:.3f} {:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f}\n".format(XY_ct[i][0], XY_ct[i][0]-1, XY_ct[i][0]+1,
                                                                                     XY_ct[i][1], YerrTot_ct[i][1], YerrTot_ct[i][0], YerrUnco_ct[i][1], YerrUnco_ct[i][0], math.sqrt(YerrTot_ct[i][1] ** 2 - YerrUnco_ct[i][1] ** 2), math.sqrt(YerrTot_ct[i][0] ** 2 - YerrUnco_ct[i][0] ** 2) ))
fo.write( "{:.2f}\n".format(rho12_ct) )

fo.write("\n# RAA(Bc) versus centrality, with CMS 2017 pp and 2018 PbPb data\n")
fo.write("# x xmin xmax ycenter +stat_u -stat_d +syst_u -syst_d\n")
fo.write("# correlation factor between bins 1 and 2")
fo.write("# 'stat' uncertainty is the part of the total uncertainty that is not bin-to-bin correlated (mostly statistical). 'syst' uncertainty is such that total uncertainty is (stat^2+syst^2) \n")
for i in range(0,len(XY_ct)):
    fo.write( "{:.3f} {:.3f} {:.3f} {:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f}\n".format((centMin[i]+centMax[i])/2, centMin[i], centMax[i],
                                                                                     XY_ct[i][1], YerrTot_ct[i][1], YerrTot_ct[i][0], YerrUnco_ct[i][1], YerrUnco_ct[i][0], math.sqrt(YerrTot_ct[i][1] ** 2 - YerrUnco_ct[i][1] ** 2), math.sqrt(YerrTot_ct[i][0] ** 2 - YerrUnco_ct[i][0] ** 2) ))
fo.write( "{:.2f}\n".format(rho12_ct) )

fo.write("\n# Bc cross section in pp versus pT, with CMS 2017 data\n")
fo.write("# x xmin xmax ycenter +stat_u -stat_d +syst_u -syst_d\n")
fo.write("# correlation factor between bins 1 and 2")
fo.write("# 'stat' uncertainty is the part of the total uncertainty from the fit procedure (mostly statistical). 'syst' uncertainty is such that total uncertainty is (stat^2+syst^2) \n")
for i in range(0,len(XY_ct)):
    fo.write( "{:.3f} {:.3f} {:.3f} {:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f}\n".format(XY_xspp[i][0], XY_xspp[i][0]-Xerr_xspp[i][0], XY_xspp[i][0]+Xerr_xspp[i][1]+8e-3,
                                                                                     XY_xspp[i][1], YerrTot_xspp[i][1], YerrTot_xspp[i][0], YerrFit_xspp[i][1], YerrFit_xspp[i][0], math.sqrt(max(0, YerrTot_xspp[i][1] ** 2 - YerrFit_xspp[i][1] ** 2)), math.sqrt(max(0, YerrTot_xspp[i][0] ** 2 - YerrFit_xspp[i][0] ** 2)) ))
fo.write( "{:.2f}\n".format(rho12_xspp) )

fo.write("\n# Bc cross section in PbPb versus pT, with CMS 2018 data\n")
fo.write("# x xmin xmax ycenter +stat_u -stat_d +syst_u -syst_d\n")
fo.write("# correlation factor between bins 1 and 2")
fo.write("# 'stat' uncertainty is the part of the total uncertainty from the fit procedure (mostly statistical). 'syst' uncertainty is such that total uncertainty is (stat^2+syst^2) \n")
for i in range(0,len(XY_ct)):
    fo.write( "{:.3f} {:.3f} {:.3f} {:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f} +{:.3f} -{:.3f}\n".format(XY_xsPbPb[i][0], XY_xsPbPb[i][0]-Xerr_xsPbPb[i][0], XY_xsPbPb[i][0]+Xerr_xsPbPb[i][1]+8e-3,
                                                                                     XY_xsPbPb[i][1], YerrTot_xsPbPb[i][1], YerrTot_xsPbPb[i][0], YerrFit_xsPbPb[i][1], YerrFit_xsPbPb[i][0], math.sqrt(max(0, YerrTot_xsPbPb[i][1] ** 2 - YerrFit_xsPbPb[i][1] ** 2)), math.sqrt(max(0, YerrTot_xsPbPb[i][0] ** 2 - YerrFit_xsPbPb[i][0] ** 2)) ))
fo.write( "{:.2f}\n".format(rho12_xsPbPb) )

fo.close()
