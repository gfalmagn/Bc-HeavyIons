from ROOT import *
import matplotlib.pyplot as plt
import math

#to do: fit hsum independently, modify to make callable in a C++ macro

#open file
file = TFile('../BDT/BDT_InputTree_pp.root')
#get trees
bkgBcMtree = file.Get('bkgBCMASS')

#define sideband histograms
Nbin = 24
low = 3.2
high = 7.3
hleft = TH1F('hleft', 'Left Sideband;trimuon mass [GeV];counts', Nbin,low, high)
hright = TH1F('hright', 'Right Sideband;trimuon mass [GeV];counts',Nbin,low,high)

#fill histograms
for e in bkgBcMtree:
    if e.QQ_M < 3.1:
        hleft.Fill(e.Bc_M,e.weight)
    if e.QQ_M > 3.1:
        hright.Fill(e.Bc_M,e.weight)

#aesthetics- colors and markers
hleft.SetMarkerColor(2)
hleft.SetLineColor(2)
hleft.SetMarkerStyle(23)
hright.SetMarkerColor(3)
hright.SetLineColor(3)
hright.SetMarkerStyle(8)

#Define fitting function
#par = Amplitude of erf, move left/right, steepness, exponential decay parameter, exponential shift
#Define the fitting function
erfexp = TF1('erfexp', '[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3])',low,high)

#fit both histograms
#Guess
erfexp.SetParameters(250,3,4,5,4)
#fix the X intercept of erfexp
#erfexp.SetParLimits(1,3.31,3.41)
erfexp.FixParameter(1,3.28)
#ignore the exponential shift by setting it to 0
erfexp.FixParameter(4,0)
#Fit
hleft.Fit('erfexp','M','',low,high)
#guess
erfexp.SetParameters(160,3.5,5,2,4)
#fix the X intercept of erfexp
#erfexp.SetParLimits(1,3.31,3.41)
erfexp.FixParameter(1,3.28)
#ignore the exponential shift by setting it to 0
erfexp.FixParameter(4,0)
#We want to push the right sideband fit up, so we set a lower limit on the displacement of the exponential, [4]
#erfexp.SetParLimits(4,7,8)
#fit
hright.Fit('erfexp','M','',low,high)

#Access fit parameters and average them, then define a new function with the averaged parameters.  This should be a reasonable signal region fit
leftfit = hleft.GetFunction('erfexp')
rightfit = hright.GetFunction('erfexp')
rightfit.SetLineColor(3)

leftN = hleft.Integral(hleft.FindBin(3.3),hleft.FindBin(7.3))
rightN = hright.Integral(hright.FindBin(3.3),hright.FindBin(7.3))

"""#Arithmetic parameter average-unweighted
arithavgp0 = (leftfit.GetParameter(0) + rightfit.GetParameter(0))/(2)
arithavgp1 = (leftfit.GetParameter(1) + rightfit.GetParameter(1))/(2)
arithavgp2 = (leftfit.GetParameter(2) + rightfit.GetParameter(2))/(2)
arithavgp3 = (leftfit.GetParameter(3) + rightfit.GetParameter(3))/(2)
arithavgp4 = (leftfit.GetParameter(4) + rightfit.GetParameter(4))/(2)"""

#Arithmetic parameter average- weighted
arithavgp0 = (leftN*leftfit.GetParameter(0) + rightN*rightfit.GetParameter(0))/(leftN + rightN)
arithavgp1 = (leftN*leftfit.GetParameter(1) + rightN*rightfit.GetParameter(1))/(leftN + rightN)
arithavgp2 = (leftN*leftfit.GetParameter(2) + rightN*rightfit.GetParameter(2))/(leftN + rightN)
arithavgp3 = (leftN*leftfit.GetParameter(3) + rightN*rightfit.GetParameter(3))/(leftN + rightN)
arithavgp4 = (leftN*leftfit.GetParameter(4) + rightN*rightfit.GetParameter(4))/(leftN + rightN)

arithSignal_fit = TF1('arithSignal_fit', 'erfexp',low,high)
arithSignal_fit.SetParameters(arithavgp0,arithavgp1,arithavgp2,arithavgp3,arithavgp4)
arithSignal_fit.SetLineColor(4)

"""#geometric average- unweighted
avgp0 = math.sqrt(leftfit.GetParameter(0)*rightfit.GetParameter(0))
avgp1 = math.sqrt(leftfit.GetParameter(1)*rightfit.GetParameter(1))
avgp2 = math.sqrt(leftfit.GetParameter(2)*rightfit.GetParameter(2))
avgp3 = math.sqrt(leftfit.GetParameter(3)*rightfit.GetParameter(3))
avgp4 = math.sqrt(leftfit.GetParameter(4)*rightfit.GetParameter(4))"""

#geometric average- weighted
w1 = leftN/(leftN+rightN)
w2 = rightN/(leftN+rightN) 
geoavgp0 = math.exp(w1*math.log(leftfit.GetParameter(0)) + w2*math.log(rightfit.GetParameter(0)))
geoavgp1 = math.exp(w1*math.log(leftfit.GetParameter(1)) + w2*math.log(rightfit.GetParameter(1)))
geoavgp2= math.exp(w1*math.log(leftfit.GetParameter(2)) + w2*math.log(rightfit.GetParameter(2)))
geoavgp3 = math.exp(w1*math.log(leftfit.GetParameter(3)) + w2*math.log(rightfit.GetParameter(3)))
#geoavgp4 = math.exp(w1*math.log(leftfit.GetParameter(4)) + w2*math.log(rightfit.GetParameter(4)))
geoavgp4 = 0

geoSignal_fit = TF1('geoSignal_fit', 'erfexp',low,high)
geoSignal_fit.SetParameters(geoavgp0,geoavgp1,geoavgp2,geoavgp3,geoavgp4)
geoSignal_fit.SetLineColor(1)

#make a legend
legend = TLegend(0.6,0.75,0.9,0.99)
legend.AddEntry(hleft,'Left Sideband','lp')
legend.AddEntry(hright, 'Right Sideband','lp')
#legend.AddEntry(geoSignal_fit, 'Averaged parameter fit-Geometric Average', 'l')
legend.AddEntry(arithSignal_fit,'Averaged Parameter Fit', 'l')#arithmetic average

#Print chi squares to check when a fit is bonked
print('left band chisquared/NDF = ', leftfit.GetChisquare()/leftfit.GetNDF())
print('Right band chisquared/NDF = ', rightfit.GetChisquare()/rightfit.GetNDF())
print('Probability of left fit=', leftfit.GetProb())
print('Probability of right fit=', rightfit.GetProb())

#normalize the three functions
#integrate
intleft = leftfit.Integral(3.3,7.3)
intright = rightfit.Integral(3.3,7.3)
intgeosig = geoSignal_fit.Integral(3.3,7.3)
intarithsig = arithSignal_fit.Integral(3.3,7.3)

#define normalized function objects
import copy
print('left sideband fit function integral', intleft)

normleftfit = copy.deepcopy(leftfit)
normleftfit.SetParameter(0,leftfit.GetParameter(0)/intleft)
print('Integral of the normalized left fit function',normleftfit.Integral(3.3,7.3))
normrightfit = copy.deepcopy(rightfit)
normrightfit.SetParameter(0,rightfit.GetParameter(0)/intright)
normgeosigfit = copy.deepcopy(geoSignal_fit)
normgeosigfit.SetParameter(0,geoSignal_fit.GetParameter(0)/intgeosig)
normarithsigfit = copy.deepcopy(arithSignal_fit)
normarithsigfit.SetParameter(0,arithSignal_fit.GetParameter(0)/intarithsig)
#aesthetics
normleftfit.SetLineColor(2)
normrightfit.SetLineColor(3)
normgeosigfit.SetLineColor(1)
normarithsigfit.SetLineColor(4)
#create a legend for the normalized plot
normleg = TLegend(0.6,0.7,0.9,0.9)
normleg.AddEntry(normleftfit, 'Left sideband','l')
normleg.AddEntry(normrightfit, 'Right sideband','l')
#normleg.AddEntry(normgeosigfit, 'Normalized Signal-Geometric Average','l')
normleg.AddEntry(normarithsigfit, 'Arithmetic Average','l')#Signal-Arithmetic

#Create biased averages to see the evolution between the green and red curves. We use arithmetic average.  Plot them on a new canvas, since c1 has the comparison between the averaging methods
c2 = TCanvas('c2','Intermediate Curves',1500,650)
c2.Divide(2,1)
biases = [0,0.3,0.5,0.8]
colors = [3,4,5,6]
funcvect = []
normfuncvect = []
biaslegend = TLegend(0.7,0.6,1,0.8)

for bias in biases:
    biasavgp0 = (bias*leftN*leftfit.GetParameter(0) + (1-bias)*rightN*rightfit.GetParameter(0))/(bias*leftN + (1-bias)*rightN)
    #biasavgp1 = (bias*leftN*leftfit.GetParameter(1) + (1-bias)*rightN*rightfit.GetParameter(1))/(bias*leftN + (1-bias)*rightN)
    biasavgp1 = 3.3
    biasavgp2 = (bias*leftN*leftfit.GetParameter(2) + (1-bias)*rightN*rightfit.GetParameter(2))/(bias*leftN + (1-bias)*rightN)
    biasavgp3 = (bias*leftN*leftfit.GetParameter(3) + (1-bias)*rightN*rightfit.GetParameter(3))/(bias*leftN + (1-bias)*rightN)
    biasavgp4 = (bias*leftN*leftfit.GetParameter(4) + (1-bias)*rightN*rightfit.GetParameter(4))/(bias*leftN +(1-bias)*rightN)
    biasSignal_fit = TF1('biasSignal_fit', 'erfexp',low,high)
    biasSignal_fit.SetParameters(biasavgp0,biasavgp1,biasavgp2,biasavgp3,biasavgp4)
    intbias = biasSignal_fit.Integral(3.3,7.3)
    normbiasfit = copy.deepcopy(biasSignal_fit)
    normbiasfit.SetParameter(0,biasSignal_fit.GetParameter(0)/intbias)
    normbiasfit.SetLineColor(colors[biases.index(bias)])
    biasSignal_fit.SetLineColor(colors[biases.index(bias)])
    funcvect.append(biasSignal_fit)
    normfuncvect.append(normbiasfit)
    biaslegend.AddEntry(biasSignal_fit,'bias= '+str(bias),'l')

unnormlabel2 = TPaveLabel(-0.7, 0.7,-0.2,0.9, 'Unnormalized Intermediate Curves')
normlabel2 = TPaveLabel(-0.7,0.7,-0.2,0.9, 'Normalized Intermediate Curves')

c2.cd(2)
leftfit.GetYaxis().SetRangeUser(0,1.1*leftfit.GetMaximum())
leftfit.Draw()
biaslegend.AddEntry(leftfit, 'Left Sideband fit eg bias = 1','l')
unnormlabel2.Draw()
normlabel2.Draw()
for f in funcvect:
    f.Draw('SAME')

c2.cd(1)
normlabel2.Draw()
normleftfit.GetYaxis().SetRangeUser(0,1.1*normleftfit.GetMaximum())
normleftfit.Draw()
for f in normfuncvect:
    f.Draw('SAME')

biaslegend.Draw()
c2.Modified()
c2.Update()
c2.SaveAs('figs/BiasedFits.pdf')

#define canvas, pads, and labels
c1 = TCanvas ('c1', 'Sideband fitting',1000,650)
c1.Divide(2,2)
c1label1 = TPaveLabel(-0.7,0.7,-0.2,0.9, 'Left Sideband')
c1label2 = TPaveLabel(-0.7,0.7,-0.2,0.9,'Right Sideband')
c1label3 = TPaveLabel(-0.7,0.7,-0.2,0.9,'Unnormalized Average Fits')
c1label4 = TPaveLabel(-0.7,0.7,-0.2,0.9,'Normalized Average Fits')

#call draw functions, update canvas, save to pdf
c1.cd(1)
gStyle.SetOptStat(1110)
hleft.Draw()
c1label1.Draw()
c1.cd(1).Update()

c1.cd(2)
gStyle.SetOptStat(1110)
hright.Draw()
c1label2.Draw()

c1.cd(3).SetTopMargin(0.01)
#Stats box: print rms, mean with error, number of entries (no error)
gStyle.SetOptStat(0)
hleft.SetTitle(';trimuon mass [GeV];counts')
hright.SetTitle(';trimuon mass [GeV];counts')
arithSignal_fit.GetHistogram().SetTitle(';trimuon mass [GeV];counts')
hleft.Draw()
hright.Draw('SAME')
legend.Draw('SAME')
#geoSignal_fit.Draw('SAME')
arithSignal_fit.Draw('SAME')
c1label3.Draw()

c1.cd(3).Modified()
c1.cd(3).Update()
c1.cd(3).SaveAs("figs/ComparedSidebandFits.pdf")

c1.cd(4)
normleftfit.GetYaxis().SetRangeUser(0,1.1*normleftfit.GetMaximum())
normleftfit.SetTitle("Normalised fits")
normleftfit.GetHistogram().GetXaxis().SetTitle("trimuon mass [GeV]")
normleftfit.GetHistogram().GetYaxis().SetTitle("a. u.")
normleftfit.Draw()
normrightfit.Draw('SAME')
#normgeosigfit.Draw('SAME')
normarithsigfit.Draw('SAME')
normleg.Draw('SAME')
c1label4.Draw()

c1.Modified()
c1.Update()
c1.SaveAs('figs/FittedSidebands.pdf')

#Compare the final fit with the simple sum of left and right fits to see how much this whole shebang has actually gained us
c3 = TCanvas('c3','Sum method vs Averaged parameter fit method',1000,650)
c3.SetTopMargin(0.01)
sum_fit = TF1('sum_fit','[0]*TMath::Erf((x-[1])/[2])*exp(-(x-[4])/[3]) + [5]*TMath::Erf((x-[6])/[7])*exp(-(x-[8])/[9])',low,high)
sum_fit.SetParameter(0,leftfit.GetParameter(0))
sum_fit.SetParameter(1,leftfit.GetParameter(1))
sum_fit.SetParameter(2,leftfit.GetParameter(2))
sum_fit.SetParameter(3,leftfit.GetParameter(3))
sum_fit.SetParameter(4,leftfit.GetParameter(4))
sum_fit.SetParameter(5,rightfit.GetParameter(0))
sum_fit.SetParameter(6,rightfit.GetParameter(1))
sum_fit.SetParameter(7,rightfit.GetParameter(2))
sum_fit.SetParameter(8,rightfit.GetParameter(4))
sum_fit.SetParameter(9,rightfit.GetParameter(3))
sum_fit.GetHistogram().SetTitle('')

#Sum the histograms and draw this as well.
hsum = hleft+hright
datanormsigfit = copy.deepcopy(normarithsigfit)
Nleft = hleft.Integral(hleft.FindBin(3.3),hleft.FindBin(7.3))
Nright = hright.Integral(hright.FindBin(3.3),hright.FindBin(7.3))
Nsum = hsum.Integral(hsum.FindBin(3.3),hsum.FindBin(7.3))

datanormsigfit.SetParameter(0,datanormsigfit.GetParameter(0) * Nsum*hleft.GetBinWidth(3) / datanormsigfit.Integral(3.3,7.3))
sum_fit.SetParameter(0,sum_fit.GetParameter(0) * Nsum*hleft.GetBinWidth(3) / sum_fit.Integral(3.3,7.3))
#sum_fit.GetYaxis().SetRangeUser(0,1.1*sum_fit.GetMaximum())

##Fit hsum independently
#erfexp.FixParameter(1,3.3)
#erfexp.FixParameter(4,0)
#hsum.Fit('erfexp','M','',low,high)
#The function that is the fit of hsum, DISTINCT from the sumfit, which is the sum of the left and right fits
#fitsum = hsum.GetFunction('erfexp')
#fitsum.SetLineColor(kRed)
#fitsum.SetParameter(0,(Nleft+Nright)*hleft.GetBinWidth(3)/fitsum.Integral())

hsum.SetTitle(';trimuon mass [GeV];counts')
hsum.SetMarkerColor(kBlack)
hsum.SetLineColor(kBlack)
hsum.SetLineWidth(2)
hsum.Draw()
sum_fit.Draw('same')
#fitsum.Draw('SAME') #don't draw the fit of the summed sidebands
datanormsigfit.Draw('SAME')

#Add Legend
Sumlegend = TLegend(0.55,0.69,0.9,0.99)
Sumlegend.AddEntry(datanormsigfit, 'Averaged Parameter fit','l')
Sumlegend.AddEntry(sum_fit, 'Sum of sidebands fits', 'l')
#Sumlegend.AddEntry(fitsum, 'Fit of the summed histogram', 'l')
Sumlegend.AddEntry(hsum, 'Sum of sidebands histograms','lpe')
Sumlegend.Draw('SAME')

"""#normalize the sum
intsum = sum_fit.Integral(3.3,7.3)
normsumfit = copy.deepcopy(sum_fit)
normsumfit.SetParameter(0,sum_fit.GetParameter(0)/intsum)
normsumfit.SetParameter(5,sum_fit.GetParameter(5)/intsum)
normsumfit.GetYaxis().SetRangeUser(0,1.1*normsumfit.GetMaximum())

#Add Legend
normsumlegend= TLegend(0.7,0.6,1,0.9)
normsumlegend.AddEntry(normarithsigfit, 'Averaged Parameter fit', 'l')
normsumlegend.AddEntry(normsumfit, 'Sum of left and right fits', 'l')
c3.cd(2)
c3label2.Draw()
normsumfit.Draw()
normarithsigfit.Draw('SAME')
normsumlegend.Draw('SAME')"""

c3.Modified()
c3.Update()
c3.SaveAs('figs/Sidebands_NewVsOld.pdf')

#Save the final fit function (arithSignal_fit) to an external root file, 'SideBandFit.root'
outfile = TFile('SideBandFit.root','RECREATE', 'Side Band Fit Function')
outfile.Write('arithSignal_fit')
