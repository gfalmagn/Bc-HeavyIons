#!/bin/bash

#
mkdir -p plotRAA

#
g++ cmsRAAdraw.C $(root-config --cflags --libs) -g -o cmsRAAdraw.exe || return 1 

CENTPbPbMIN=0
CENTPbPbMAX=100
OUTPUTFILERAAMB="ROOTfilesCent0100/outputRAAMB.root"
OUTPUTFILERAA="ROOTfilesCent0100/outputRAA.root"
set -x
#                 -                   -                 centmin         centmax         iD  ihad  iB  injs  ind  ipjs iBc  iBs  psi2 iUps -
./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  0   0     0   0     0    1    1    0    1    1    #  Bc/promptJpsi/Psi2S/UpsilonNS
./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   1     1   0     0    0    1    1    0    0    #  Bc/Bs/D/B/charged
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   1     1   1     1    0    0    0    0    0    #  charged/charm/beauty
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   0     0   0     0    0    0    0    0    0    #  charm
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  0   0     1   1     1    0    0    0    0    0    #  beauty
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   0     1   1     1    0    0    0    0    0    #  charm/beauty
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   1     0   0     0    0    0    0    0    0    #  charged/charm
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  0   1     1   1     1    0    0    0    0    0    #  charged/beauty
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   0     0   1     0    0    0    0    0    0    #  charm/non-promptjpsi
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   0     0   0     1    0    0    0    0    0    #  charm/non-promptD
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   0     0   0     0    1    0    0    0    0    #  charm/promptjpsi
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   1     1   0     1    0    0    0    0    0    #  charged/charm/B/npD
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   1     1   1     0    0    0    0    0    0    #  charged/charm/B/npjpsi
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   1     1   0     0    0    0    0    0    0    #  charged/charm/B
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   1     0   0     1    0    0    0    0    0    #  charged/charm/npD
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   1     0   1     0    0    0    0    0    0    #  charged/charm/npjpsi
set +x

CENTPbPbMIN=0
CENTPbPbMAX=10
OUTPUTFILERAAMB="ROOTfilesCent010/outputRAAMB.root"
OUTPUTFILERAA="ROOTfilesCent010/outputRAA.root"
set -x
#                 -                   -                 centmin         centmax         iD  ihad  iB  injs  ind  ipjs iBc  iBs  -
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   0     0   0     0    0    0    0    0    0    #  charm
#./cmsRAAdraw.exe  "$OUTPUTFILERAAMB"  "$OUTPUTFILERAA"  "$CENTPbPbMIN"  "$CENTPbPbMAX"  1   1     0   0     0    0    0    0    0    0    #  charged/charm
set +x

rm cmsRAAdraw.exe
