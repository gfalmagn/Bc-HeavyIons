import re #regular expressions

fcuts = open("../helpers/Cuts.h", 'r')
fi = [open("Nmin1efficiencies_pp.txt", 'r') , open("Nmin1efficiencies_PbPb.txt", 'r')]
fo = [open("N-1effTable_pp.tex", "w") , open("N-1effTable_PbPb.tex", "w")]
li = [f.readlines() for f in fi]

cuts = []
vart = []
tcuts = fcuts.read()
for col in [0,1]:
    cuts.append( [ re.findall("\d+\.\d+" , (re.findall("_vtxProb_cut = \d+\.\d+",tcuts)[0]))[0] ,
                   re.findall("\d+\.\d+" , (re.findall("_QQvtxProb_cut = \d+\.\d+",tcuts)[0]))[0] ,                   
                   re.findall("\d+\.\d+:" if (col==0) else "\d+\.\d+;", (re.findall("float _alpha3D_cut\(bool ispp=true\)\{return ispp\?\d+\.\d+:\d+\.\d+;\}",tcuts)[0]))[0][:-1],
                   re.findall("\d+\.\d+:" if (col==0) else "\d+\.\d+;", (re.findall("float _alpha_cut\(bool ispp=true\)\{return ispp\?\d+\.\d+:\d+\.\d+;\}",tcuts)[0]))[0][:-1],
                   re.findall("\d+\.\d+" , (re.findall("_ctauSignif3D_cut = \d+\.\d+",tcuts)[0]))[0] ,                   
                   re.findall("\d+\.\d+" , (re.findall("_ctauSignif_cut = \d+\.\d+",tcuts)[0]))[0] ,                   
                   re.findall("\d+\.\d+" , (re.findall("_QQdca_cut = \d+\.\d+",tcuts)[0]))[0] ,                   
                   re.findall("\d+:" if (col==0) else "\d+;", (re.findall("float _BcCorrM_cut\(bool ispp=true\)\{return ispp\?\d+:\d+;\}",tcuts)[0]))[0][:-1],
               ] );
    
print cuts
    
vart.append( [
    ("TrimuVProb" , 
     "\\rule{{0pt}}{{8mm}}(1) trimuon VtxProb & $>"+cuts[0][0]+"$ & {:.{prec0}f} &\\shortstack{{(from $prob>0.005$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & \\shortstack{{(from $prob>0.005$)\\\\$>${:.{prec3}f}}}\\\\\n"),
    ("DimuVProb",
     "\\rule{{0pt}}{{8mm}}(2) \\PJGy VtxProb & $>"+cuts[0][1]+"$ & {:.{prec0}f} & \\shortstack{{(from $prob>0.002$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & \\shortstack{{(from $prob>0.002$)\\\\$>${:.{prec3}f}}}\\\\\n" ),
    ("TrimuVProbDimuVProb",
     "\\multicolumn{{2}}{{c||}}{{(1) and (2)}} & {:.{prec0}f} & $>${:.{prec1}f} & {:.{prec2}f} & $>${:.{prec3}f}\\\\\n" ),
    ("alpha3D",
     "\\rule{{0pt}}{{8mm}}(3) $\\alpha_{{3D}}$ [rad] & $<"+cuts[0][2]+"$ & {:.{prec0}f} & \\shortstack{{(from $\\alpha_{{3D}}<1.57$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & \\shortstack{{(from $\\alpha_{{3D}}<1.57$)\\\\$>${:.{prec3}f}}}\\\\\n" ),
    ("alpha2D",
     "\\rule{{0pt}}{{8mm}}(4) $\\alpha_{{2D}}$ [rad] & $<"+cuts[0][3]+"$ & {:.{prec0}f} & \\shortstack{{(from $\\alpha_{{2D}}<1.37$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & \\shortstack{{(from $\\alpha_{{2D}}<1.37$)\\\\$>${:.{prec3}f}}}\\\\\n" ),
    ("alpha2Dalpha3D",
     "\\multicolumn{{2}}{{c||}}{{(3) and (4)}} & {:.{prec0}f} & $>${:.{prec1}f} & {:.{prec2}f} & $>${:.{prec3}f}\\\\\n" ),
    ("tauSignif3D",
     "\\rule{{0pt}}{{8mm}}(5) $\\tau_{{3D}}/\\sigma_{{\\tau_{{3D}}}}$ & $>"+cuts[0][4]+"$ & {:.{prec0}f} & \\shortstack{{(from $\\tau/\\sigma>0$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & \\shortstack{{(from $\\tau/\\sigma>0$)\\\\$>${:.{prec3}f}}}\\\\[2pt]\n" ),
    ("tauSignif2D",
     "\\rule{{0pt}}{{8mm}}(6) $\\tau_{{2D}}/\\sigma_{{\\tau_{{2D}}}}$ & $>"+cuts[0][5]+"$ & {:.{prec0}f} & \\shortstack{{(from $\\tau/\\sigma>0$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & \\shortstack{{(from $\\tau/\\sigma>0$)\\\\$>${:.{prec3}f}}}\\\\[2pt]\n" ),
    ("tauSignif3DtauSignif2D",
     "\\multicolumn{{2}}{{c||}}{{(5) and (6)}} & {:.{prec0}f} & $>${:.{prec1}f} & {:.{prec2}f} & $>${:.{prec3}f}\\\\\n" ),
    ("dca",
     "(7) $dca(\\PJGy)$ & $<"+cuts[0][6]+"$~mm & {:.{prec0}f} &{:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("dz",
     "(8) $max_i(d_z(\\mu_i))$ & $<6$~mm & {:.{prec0}f} & {:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ), 
    ("CorrM",
     "(9) $m_{{corr}}$ & $<"+cuts[0][7]+"\\GeV$ & {:.{prec0}f} & {:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("3global",
     "\\multicolumn{{2}}{{c||}}{{(10) 3 \\textit{{hybrid-soft}} in loose acc.}} & {:.{prec0}f} & {:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("TrimuVProbDimuVProbTrimuVProbDimuVProbalpha3Dalpha2Dalpha2Dalpha3DtauSignif3DtauSignif2DtauSignif3DtauSignif2DdcaCorrMdz",
     "\\multicolumn{{2}}{{c||}}{{cuts (1) to (9)}} & \\textbf{{{:.{prec0}f}}} & \\textbf{{$>${:.{prec1}f}}} & \\textbf{{{:.{prec2}f}}} & \\textbf{{$>${:.{prec3}f}}}\\\\\n" )
])

vart.append( [
    ("TrimuVProb" , 
     "\\rule{{0pt}}{{8mm}}(1) trimuon VtxProb & $>"+cuts[1][0]+"$ & {:.{prec0}f} &\\shortstack{{(from $prob>0.005$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & {:.{prec3}f}\\\\\n"),
    ("DimuVProb",
     "\\rule{{0pt}}{{8mm}}(2) \\PJGy VtxProb & $>"+cuts[1][1]+"$ & {:.{prec0}f} & \\shortstack{{(from $prob>0.002$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("TrimuVProbDimuVProb",
     "\\multicolumn{{2}}{{c||}}{{(1) and (2)}} & {:.{prec0}f} & $>${:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("alpha3D",
     "\\rule{{0pt}}{{8mm}}(3) $\\alpha_{{3D}}$ [rad] & $<"+cuts[1][2]+"$ & {:.{prec0}f} & \\shortstack{{(from $\\alpha_{{3D}}<1.16$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("alpha2D",
     "\\rule{{0pt}}{{8mm}}(4) $\\alpha_{{2D}}$ [rad] & $<"+cuts[1][3]+"$ & {:.{prec0}f} & \\shortstack{{(from $\\alpha_{{2D}}<1.57$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("alpha2Dalpha3D",
     "\\multicolumn{{2}}{{c||}}{{(3) and (4)}} & {:.{prec0}f} & $>${:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("tauSignif3D",
     "\\rule{{0pt}}{{8mm}}(5) $\\tau_{{3D}}/\\sigma_{{\\tau_{{3D}}}}$ & $>"+cuts[1][4]+"$ & {:.{prec0}f} & \\shortstack{{(from $\\tau/\\sigma>0$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & {:.{prec3}f}\\\\[2pt]\n" ),
    ("tauSignif2D",
     "\\rule{{0pt}}{{8mm}}(6) $\\tau_{{2D}}/\\sigma_{{\\tau_{{2D}}}}$ & $>"+cuts[1][5]+"$ & {:.{prec0}f} & \\shortstack{{(from $\\tau/\\sigma>0$)\\\\$>${:.{prec1}f}}} & {:.{prec2}f} & {:.{prec3}f}\\\\[2pt]\n" ),
    ("tauSignif3DtauSignif2D",
     "\\multicolumn{{2}}{{c||}}{{(5) and (6)}} & {:.{prec0}f} & $>${:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("dca",
     "(7) $dca(\\PJGy)$ & $<"+cuts[1][6]+"$~mm & {:.{prec0}f} &{:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("dz",
     "(8) $max_i(d_z(\\mu_i))$ & $<6$~mm & {:.{prec0}f} & {:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ), 
    ("CorrM",
     "(9) $m_{{corr}}$ & $<"+cuts[1][7]+"\\GeV$ & {:.{prec0}f} & {:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("3global",
     "\\multicolumn{{2}}{{c||}}{{(10) 3 \\textit{{hybrid-soft}} in loose acc.}} & {:.{prec0}f} & {:.{prec1}f} & {:.{prec2}f} & {:.{prec3}f}\\\\\n" ),
    ("TrimuVProbDimuVProbTrimuVProbDimuVProbalpha3Dalpha2Dalpha2Dalpha3DtauSignif3DtauSignif2DtauSignif3DtauSignif2DdcaCorrMdz",
     "\\multicolumn{{2}}{{c||}}{{cuts (1) to (9)}} & \\textbf{{{:.{prec0}f}}} & \\textbf{{$>${:.{prec1}f}}} & \\textbf{{{:.{prec2}f}}} & \\textbf{{{:.{prec3}f}}}\\\\\n" )
])

samp = [["signal","fakeJpsi","bMCcorr","flipJpsi"],#pp
        ["signal","fakeJpsi","bMC","PromptMC"]]#pbpb

for f,fout, lines, sam, var in zip(fi,fo,li,samp,vart):
    lookhere = [False]*4
    eff = []
    for v, text in var:
        res = [0]*4
        
        for line in lines:
            for i in range(4):
                if sam[i] in line:
                    lookhere[i] = True
                elif (sam[i] and ("sample" in f)):
                    lookhere[i] = False #end of interesting text block
                if (lookhere[i] and ((" "+v+" ") in line)):
                    r = re.findall("\d+\.\d+$", line)
                    if (len(r)==0 and  (" 0" in line or "e-0" in line)):
                        r = [0]
                    res[i] = 100*float(r[0]) #find number with decimals at end of this line
                    lookhere[i] = False #security
                
        eff.append((v, text, res))
        print v, res

    f.close()

    for v, text, e in eff:
        #two decimals if value <1, else 1 decimal 
        fout.write(text.format( e[0],e[1],e[2],e[3], prec0=1+(e[0]<1),prec1=1+(e[1]<1),prec2=1+(e[2]<1),prec3=1+(e[3]<1)) )
        fout.write("\\hline")
    
    fout.close()


