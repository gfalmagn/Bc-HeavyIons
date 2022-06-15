from hepdata_lib import Submission
submission = Submission()

abst = open("abstract.txt", "w")
abst.write("${\rm B_c^+}$ mesons are observed in lead-lead and proton-proton collisions at a center-of-mass energy per nucleon pair of $\sqrt{s_{\rm NN}}=5.02~\mathrm{TeV}$, via the ${\rm B_c^+} \rightarrow ({\rm J}/\psi \rightarrow \mu^+ \mu^-) \mu^+ \nu_\mu$ decay and using 2017 and 2018 data from the CMS detector. The ${\rm B_c^+}$ production cross sections and nuclear modification factor are measured in coarse bins of the trimuon transverse momentum and of the collision centrality. \nThe insights gained from this first observation in heavy ion collisions will further the understanding of the interplay of suppression and enhancement mechanisms in the production of heavy-flavor mesons in the hot and dense matter created in these collisions. \nThe \mbox{${\rm B_c^+}$ meson} is observed to be less suppressed than other quarkonia and most open \mbox{heavy-flavor} mesons, hinting at effects of the hot deconfined medium that \mbox{contribute} to its \mbox{production}.")
submission.read_abstract("abstract.txt")
submission.add_link("Webpage with all figures and tables", "https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/HIN-20-004/index.html")
submission.add_link("cds record", "https://cds.cern.ch/record/2772382?ln=fr")
#submission.add_record_id(1657397, "inspire")

### descriptions and keywords
from hepdata_lib import Table
tableXS = Table("Bc cross sections")
tableXS.description = "The $B_c$ meson production (pp-equivalent) cross-section times branching fraction of the $B_c\\rightarrow (J/\psi \\rightarrow \mu\mu) \mu \\nu_{\mu}$ decay in pp and PbPb collisions. The used kinematic variables correspond to those of the trimuon final state. The two $p_T$ bins correspond to different rapidity ranges. The total uncertainty is decomposed in an uncertainty from the fit and one representing all other sources. The markers of the $p_T^{\mu\mu\mu}$ bins are placed according to the Lafferty-Wyatt prescription."
tableXS.location = "Data from Figure 3"
tableXS.keywords["observables"] = ["D2SIG/DPT/DYRAP"]
tableXS.keywords["reactions"] = ["P P --> B/C+ X", "PB PB --> B/C+ X", "B/C+ --> J/PSI(1S) MU+ NUMU", "J/PSI(1S) --> MU+ MU-"]
tableXS.keywords["phrases"] = ["Double Differential Cross Section","J/Psi","Beauty","B/C+","B/C+ Production","Transverse Momentum Dependence"] 

tableRAApt = Table("Bc RAA vs pT")
tableRAApt.description = "The $B_c$ meson nuclear modification factor in PbPb collisions, in $p_T^{\mu\mu\mu}$ bins corresponding to different trimuon rapidity ranges. The total uncertainty is decomposed in a bin-to-bin-uncorrelated uncertainty and one fully correlated along the two bins. The markers of the $p_T^{\mu\mu\mu}$ bins are placed at the average of the Lafferty-Wyatt prescriptions applied to the pp and PbPb spectra."
tableRAApt.location = "Data from Figure 4 left"
tableRAApt.keywords["observables"] = ["RAA"]
tableRAApt.keywords["reactions"] = ["P P --> B/C+ X", "PB PB --> B/C+ X", "B/C+ --> J/PSI(1S) MU+ NUMU", "J/PSI(1S) --> MU+ MU-"]
tableRAApt.keywords["phrases"] = ["B/C+ Nuclear Modification Factor","J/Psi","Beauty","B/C+","B/C+ Production","Transverse Momentum Dependence","R measurement"]

tableRAAcent = Table("Bc RAA vs centrality")
tableRAAcent.description = "The $B_c$ meson nuclear modification factor in PbPb collisions, in centrality bins, integrated over the studied kinematic range. The cut on the trimuon rapidity depends on the trimuon $p_T$. The total uncertainty is decomposed in a bin-to-bin-uncorrelated uncertainty and one fully correlated along the two bins. The centrality bin markers are placed at the minimum bias average number of participants $N_{part}$."
tableRAAcent.location = "Data from Figure 4 right"
tableRAAcent.keywords["observables"] = ["RAA"]
tableRAAcent.keywords["reactions"] = ["P P --> B/C+ X", "PB PB --> B/C+ X", "B/C+ --> J/PSI(1S) MU+ NUMU", "J/PSI(1S) --> MU+ MU-"]
tableRAAcent.keywords["phrases"] = ["B/C+ Nuclear Modification Factor","J/Psi","Beauty","B/C+","B/C+ Production","Centrality Dependence","NPART","R measurement"]

### Read the files
import numpy as np
import itertools

with open('HepDataTable_RAA.txt') as infile:
    data_ppxs = np.genfromtxt(itertools.islice(infile, 24, 26, None), dtype=float)
with open('HepDataTable_RAA.txt') as infile:
    data_PbPbxs = np.genfromtxt(itertools.islice(infile, 31, 33, None), dtype=float)
with open('HepDataTable_RAA.txt') as infile:
    data_raapt = np.genfromtxt(itertools.islice(infile, 3, 5, None), dtype=float)
with open('HepDataTable_RAA.txt') as infile:
    data_raanpart = np.genfromtxt(itertools.islice(infile, 10, 12, None), dtype=float)
with open('HepDataTable_RAA.txt') as infile:
    data_raacent = np.genfromtxt(itertools.islice(infile, 17, 19, None), dtype=float)
with open('HepDataTable_RAA.txt') as infile:
    rho12_raapt = np.genfromtxt(itertools.islice(infile, 5, 6, None), dtype=float)
with open('HepDataTable_RAA.txt') as infile:
    rho12_raacent = np.genfromtxt(itertools.islice(infile, 12, 13, None), dtype=float)
with open('HepDataTable_RAA.txt') as infile:
    rho12_ppxs = np.genfromtxt(itertools.islice(infile, 26, 27, None), dtype=float)
with open('HepDataTable_RAA.txt') as infile:
    rho12_PbPbxs = np.genfromtxt(itertools.islice(infile, 33, 34, None), dtype=float)
#print data_ppxs
#print data_PbPbxs
#print data_raapt
#print data_raacent

#Create the actual tables
from hepdata_lib import Variable,Uncertainty

############## pp and PbPb XS. pp in the first two points, PbPb in the last two
#pT
pt = Variable("trimuon $p_T$", is_independent=True, is_binned=False, units="GeV")
pt.values = np.concatenate((data_ppxs[:,0], data_PbPbxs[:,0]))

# XS values
XSpp1 = Variable("$BF \\times \\frac{d\\sigma_{pp}^{B_c}}{dp_T^{\mu\mu\mu} dy^{\mu\mu\mu}}$", is_independent=False, is_binned=False, units="pb/GeV")
XSpp1.values = [ data_ppxs[0,3], "--",  "--",  "--" ]
XSpp1.add_qualifier("RE", "B/C+ --> (J/PSI(1S) --> MU+ MU-) MU+ NUMU") #"branching fraction times pT-differential cross section"
XSpp1.add_qualifier("SQRT(S)/NUCLEON", 5020, "GeV")
XSpp1.add_qualifier("CENTRALITY", "--")
XSpp1.add_qualifier("ABS(RAP)", "1.3-2.3")
XSpp1.add_qualifier("bin-to-bin correlation factor", float(rho12_ppxs))
XSpp2 = Variable("$BF \\times \\frac{d\\sigma_{pp}^{B_c}}{dp_T^{\mu\mu\mu} dy^{\mu\mu\mu}}$", is_independent=False, is_binned=False, units="pb/GeV")
XSpp2.values = [ "--", data_ppxs[1,3], "--",  "--" ]
XSpp2.add_qualifier("RE", "B/C+ --> (J/PSI(1S) --> MU+ MU-) MU+ NUMU")
XSpp2.add_qualifier("SQRT(S)/NUCLEON", 5020, "GeV")
XSpp2.add_qualifier("CENTRALITY", "--")
XSpp2.add_qualifier("ABS(RAP)", "0-2.3")
XSpp2.add_qualifier("bin-to-bin correlation factor", float(rho12_ppxs))
XSPbPb1 = Variable("$BF \\times \\frac{1}{N_{MB} T_{PbPb}} \\frac{dN_{PbPb}^{B_c}}{dp_T^{\mu\mu\mu} dy^{\mu\mu\mu}}$", is_independent=False, is_binned=False, units="pb/GeV")
XSPbPb1.values = [ "--", "--", data_PbPbxs[0,3],  "--" ]
XSPbPb1.add_qualifier("RE", "B/C+ --> (J/PSI(1S) --> MU+ MU-) MU+ NUMU")
XSPbPb1.add_qualifier("SQRT(S)/NUCLEON", 5020, "GeV")
XSPbPb1.add_qualifier("CENTRALITY", "0-90%")
XSPbPb1.add_qualifier("ABS(RAP)", "1.3-2.3")
XSPbPb1.add_qualifier("bin-to-bin correlation factor", float(rho12_PbPbxs))
XSPbPb2 = Variable("$BF \\times \\frac{1}{N_{MB} T_{PbPb}} \\frac{dN_{PbPb}^{B_c}}{dp_T^{\mu\mu\mu} dy^{\mu\mu\mu}}$", is_independent=False, is_binned=False, units="pb/GeV")
XSPbPb2.values = [ "--", "--", "--", data_PbPbxs[1,3] ]
XSPbPb2.add_qualifier("RE", "B/C+ --> (J/PSI(1S) --> MU+ MU-) MU+ NUMU")
XSPbPb2.add_qualifier("SQRT(S)/NUCLEON", 5020, "GeV")
XSPbPb2.add_qualifier("CENTRALITY", "0-90%")
XSPbPb2.add_qualifier("ABS(RAP)", "0-2.3")
XSPbPb2.add_qualifier("bin-to-bin correlation factor", float(rho12_PbPbxs))

#pp XS uncertainty
XSpp1totE = Uncertainty("total", is_symmetric=False)
XSpp1totE.values = [ (data_ppxs[0][5],data_ppxs[0][4]) , (0,0), (0,0), (0,0) ]
XSpp1.add_uncertainty(XSpp1totE)
XSpp1uncorE = Uncertainty("from fit", is_symmetric=False)
XSpp1uncorE.values = [ (data_ppxs[0][7],data_ppxs[0][6]) , (0,0), (0,0), (0,0) ]
XSpp1.add_uncertainty(XSpp1uncorE)
XSpp1corE = Uncertainty("not from fit", is_symmetric=False)
XSpp1corE.values = [ (data_ppxs[0][9],data_ppxs[0][8]) , (0,0), (0,0), (0,0) ]
XSpp1.add_uncertainty(XSpp1corE)

XSpp2totE = Uncertainty("total", is_symmetric=False)
XSpp2totE.values = [ (0,0), (data_ppxs[1][5],data_ppxs[1][4]) , (0,0), (0,0) ]
XSpp2.add_uncertainty(XSpp2totE)
XSpp2uncorE = Uncertainty("from fit", is_symmetric=False)
XSpp2uncorE.values = [ (0,0), (data_ppxs[1][7],data_ppxs[1][6]) , (0,0), (0,0) ]
XSpp2.add_uncertainty(XSpp2uncorE)
XSpp2corE = Uncertainty("not from fit", is_symmetric=False)
XSpp2corE.values = [ (0,0), (data_ppxs[1][9],data_ppxs[1][8]) , (0,0), (0,0) ]
XSpp2.add_uncertainty(XSpp2corE)

XSPbPb1totE = Uncertainty("total", is_symmetric=False)
XSPbPb1totE.values = [ (0,0), (0,0), (data_PbPbxs[0][5],data_PbPbxs[0][4]) , (0,0) ]
XSPbPb1.add_uncertainty(XSPbPb1totE)
XSPbPb1uncorE = Uncertainty("from fit", is_symmetric=False)
XSPbPb1uncorE.values = [ (0,0), (0,0), (data_PbPbxs[0][7],data_PbPbxs[0][6]) , (0,0) ]
XSPbPb1.add_uncertainty(XSPbPb1uncorE)
XSPbPb1corE = Uncertainty("not from fit", is_symmetric=False)
XSPbPb1corE.values = [ (0,0), (0,0), (data_PbPbxs[0][9],data_PbPbxs[0][8]) , (0,0) ]
XSPbPb1.add_uncertainty(XSPbPb1corE)

XSPbPb2totE = Uncertainty("total", is_symmetric=False)
XSPbPb2totE.values = [ (0,0), (0,0), (0,0), (data_PbPbxs[1][5],data_PbPbxs[1][4]) ]
XSPbPb2.add_uncertainty(XSPbPb2totE)
XSPbPb2uncorE = Uncertainty("from fit", is_symmetric=False)
XSPbPb2uncorE.values = [ (0,0), (0,0), (0,0), (data_PbPbxs[1][7],data_PbPbxs[1][6]) ]
XSPbPb2.add_uncertainty(XSPbPb2uncorE)
XSPbPb2corE = Uncertainty("not from fit", is_symmetric=False)
XSPbPb2corE.values = [ (0,0), (0,0), (0,0), (data_PbPbxs[1][9],data_PbPbxs[1][8]) ]
XSPbPb2.add_uncertainty(XSPbPb2corE)

#Add variables in tables
tableXS.add_variable(pt)
tableXS.add_variable(XSpp1)
tableXS.add_variable(XSpp2)
tableXS.add_variable(XSPbPb1)
tableXS.add_variable(XSPbPb2)

#Add plot
#tableXS.add_image("CrossSections.pdf")


############## RAA vs pT
#pT
ptraa = Variable("trimuon $p_T$", is_independent=True, is_binned=False)
ptraa.values = data_raapt[:,0]

# RAA values
RAApt1 = Variable("$R_{PbPb}(B_c)$", is_independent=False, is_binned=False)
RAApt1.values = [ data_raapt[0,3], "--" ]
RAApt1.add_qualifier("RE", "B/C+ --> (J/PSI(1S) --> MU+ MU-) MU+ NUMU")
RAApt1.add_qualifier("SQRT(S)/NUCLEON", 5020, "GeV")
RAApt1.add_qualifier("CENTRALITY", "0-90%")
RAApt1.add_qualifier("ABS(RAP)", "1.3-2.3")
RAApt1.add_qualifier("bin-to-bin correlation factor", float(rho12_raapt))
RAApt2 = Variable("$R_{PbPb}(B_c)$", is_independent=False, is_binned=False)
RAApt2.values = [ "--", data_raapt[1,3] ]
RAApt2.add_qualifier("RE", "B/C+ --> (J/PSI(1S) --> MU+ MU-) MU+ NUMU")
RAApt2.add_qualifier("SQRT(S)/NUCLEON", 5020, "GeV")
RAApt2.add_qualifier("CENTRALITY", "0-90%")
RAApt2.add_qualifier("ABS(RAP)", "0-2.3")
RAApt2.add_qualifier("bin-to-bin correlation factor", float(rho12_raapt))

# Uncertainties
RAApt1totE = Uncertainty("total", is_symmetric=False)
RAApt1totE.values = [ (data_raapt[0][5],data_raapt[0][4]) , (0,0) ]
RAApt1.add_uncertainty(RAApt1totE)
RAApt1uncorE = Uncertainty("uncorrelated", is_symmetric=False)
RAApt1uncorE.values = [ (data_raapt[0][7],data_raapt[0][6]) , (0,0) ]
RAApt1.add_uncertainty(RAApt1uncorE)
RAApt1corE = Uncertainty("correlated", is_symmetric=False)
RAApt1corE.values = [ (data_raapt[0][9],data_raapt[0][8]) , (0,0) ]
RAApt1.add_uncertainty(RAApt1corE)

RAApt2totE = Uncertainty("total", is_symmetric=False)
RAApt2totE.values = [ (0,0), (data_raapt[1][5],data_raapt[1][4]) ]
RAApt2.add_uncertainty(RAApt2totE)
RAApt2uncorE = Uncertainty("uncorrelated", is_symmetric=False)
RAApt2uncorE.values = [ (0,0), (data_raapt[1][7],data_raapt[1][6]) ]
RAApt2.add_uncertainty(RAApt2uncorE)
RAApt2corE = Uncertainty("correlated", is_symmetric=False)
RAApt2corE.values = [ (0,0), (data_raapt[1][9],data_raapt[1][8]) ]
RAApt2.add_uncertainty(RAApt2corE)

# Add variables in tables
tableRAApt.add_variable(ptraa)
tableRAApt.add_variable(RAApt1)
tableRAApt.add_variable(RAApt2)


############## RAA vs cent
#cent
centraa = Variable("CENTRALITY", is_independent=True, is_binned=True,units="%")
centraa.values = [(lo,hi) for (lo,hi) in zip(data_raacent[:,1],data_raacent[:,2]) ]
#npart
npartraa = Variable("NPART", is_independent=False, is_binned=False)
npartraa.values = data_raanpart[:,0]

# RAA values
RAAcent1 = Variable("$R_{PbPb}(B_c)$", is_independent=False, is_binned=False)
RAAcent1.values = [ data_raacent[0,3], "--" ]
RAAcent1.add_qualifier("RE", "B/C+ --> (J/PSI(1S) --> MU+ MU-) MU+ NUMU")
RAAcent1.add_qualifier("SQRT(S)/NUCLEON", 5020, "GeV")
RAAcent1.add_qualifier("PT", "6-35", "GeV")
RAAcent1.add_qualifier("ABS(RAP)", "1.3-2.3 if (PT<11 GeV) else 0-2.3")
RAAcent1.add_qualifier("bin-to-bin correlation factor", float(rho12_raacent))
RAAcent2 = Variable("$R_{PbPb}(B_c)$", is_independent=False, is_binned=False)
RAAcent2.values = [ "--", data_raacent[1,3] ]
RAAcent2.add_qualifier("RE", "B/C+ --> (J/PSI(1S) --> MU+ MU-) MU+ NUMU")
RAAcent2.add_qualifier("SQRT(S)/NUCLEON", 5020, "GeV")
RAAcent2.add_qualifier("PT", "6-35", "GeV")
RAAcent2.add_qualifier("ABS(RAP)", "1.3-2.3 if (PT<11 GeV) else 0-2.3")
RAAcent2.add_qualifier("bin-to-bin correlation factor", float(rho12_raacent))

# Uncertainties
npartraaE = Uncertainty("", is_symmetric=True)
npartraaE.values = [ abs(data_raanpart[0][1] - data_raanpart[0][0]) , abs(data_raanpart[1][1] - data_raanpart[1][0]) ]
npartraa.add_uncertainty(npartraaE)

RAAcent1totE = Uncertainty("total", is_symmetric=False)
RAAcent1totE.values = [ (data_raacent[0][5],data_raacent[0][4]) , (0,0) ]
RAAcent1.add_uncertainty(RAAcent1totE)
RAAcent1uncorE = Uncertainty("uncorrelated", is_symmetric=False)
RAAcent1uncorE.values = [ (data_raacent[0][7],data_raacent[0][6]) , (0,0) ]
RAAcent1.add_uncertainty(RAAcent1uncorE)
RAAcent1corE = Uncertainty("correlated", is_symmetric=False)
RAAcent1corE.values = [ (data_raacent[0][9],data_raacent[0][8]) , (0,0) ]
RAAcent1.add_uncertainty(RAAcent1corE)

RAAcent2totE = Uncertainty("total", is_symmetric=False)
RAAcent2totE.values = [ (0,0), (data_raacent[1][5],data_raacent[1][4]) ]
RAAcent2.add_uncertainty(RAAcent2totE)
RAAcent2uncorE = Uncertainty("uncorrelated", is_symmetric=False)
RAAcent2uncorE.values = [ (0,0), (data_raacent[1][7],data_raacent[1][6]) ]
RAAcent2.add_uncertainty(RAAcent2uncorE)
RAAcent2corE = Uncertainty("correlated", is_symmetric=False)
RAAcent2corE.values = [ (0,0), (data_raacent[1][9],data_raacent[1][8]) ]
RAAcent2.add_uncertainty(RAAcent2corE)

# Add variables in tables
tableRAAcent.add_variable(centraa)
tableRAAcent.add_variable(RAAcent1)
tableRAAcent.add_variable(RAAcent2)
tableRAAcent.add_variable(npartraa)


########### Add tables in submission
submission.add_table(tableXS)
submission.add_table(tableRAApt)
submission.add_table(tableRAAcent)

for table in submission.tables:
    table.keywords["cmenergies"] = [5020]

# Create archive
outdir = "BcHepData"
submission.create_files(outdir)


####### Fix the binned pT variable with low and high bin limits
import re,os
#search for string to replace
for filename in ['bc_cross_sections','bc_raa_vs_pt','bc_raa_vs_centrality']:
    with open(outdir+'/'+filename+'.yaml', 'r') as f :
        fdata = f.read()

    if filename!='bc_raa_vs_centrality':
        fdata2 = fdata.partition("independent_variables")[2] #keep everything after independent_variables
        ValsToChange = re.findall("- value: \d+\.\d*",fdata2)
        ValsChanged = []
        for i,l in enumerate(ValsToChange):
            datatmp = data_raapt if filename=='bc_raa_vs_pt' else (data_ppxs if i<2 else data_PbPbxs)
            ValsChanged.append( l[:2]+"{"+l[2:] # insert "{" 
                                + ", low: "+str(round( datatmp[i%2 ,1] ,1))+", high: "+str(round( datatmp[i%2 ,2] ,1))+"}" )  # add the low and high bin limits
        # replace string 
        for oldv, newv in zip(ValsToChange,ValsChanged):
            fdata = fdata.replace(oldv, newv)

    # remove errors for empty values #!!!!!! Needs to be updated if names of error type change!
    fdata = fdata.replace("errors:\n    - asymerror:\n        minus: 0.0\n        plus: 0.0\n      label: total\n    - asymerror:\n        minus: 0.0\n        plus: 0.0\n      label: from fit\n    - asymerror:\n        minus: 0.0\n        plus: 0.0\n      label: not from fit\n    ", "")
    fdata = fdata.replace("errors:\n    - asymerror:\n        minus: 0.0\n        plus: 0.0\n      label: total\n    - asymerror:\n        minus: 0.0\n        plus: 0.0\n      label: uncorrelated\n    - asymerror:\n        minus: 0.0\n        plus: 0.0\n      label: correlated\n    ","")
    # write updated file
    with open(outdir+'/'+filename+'.yaml', 'w') as f:
        f.write(fdata)

# Recreate archive
os.system("tar -czvf submission.tar.gz "+outdir)

