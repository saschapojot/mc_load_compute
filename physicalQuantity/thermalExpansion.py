import pickle
import numpy as np
from datetime import datetime
import sys
import re
import glob
import os
import json
from pathlib import Path
import matplotlib.pyplot as plt
#this script compite alpha

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()


funcName=sys.argv[1]
rowName=sys.argv[2]

TFolderRoot="../dataAll/"+funcName+"/"+rowName+"/jsonOutAll/"


#search directory
TVals=[]
TFileNames=[]
TStrings=[]
for TFile in glob.glob(TFolderRoot+"/T*"):
    # print(TFile)
    matchT=re.search(r"T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",TFile)
    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))
        TStrings.append("T"+matchT.group(1))

#sort T values
sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]
sortedTStrings=[TStrings[ind] for ind in sortedInds]






def compute_alpha(oneTFolder):
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFolder)
    TVal=float(matchT.group(1))
    LFilePath=oneTFolder+"/jsonData/jsonL/LData.json"
    with open (LFilePath,"r") as fptr:
        LData=json.load(fptr)
    LVec=LData["L"]

    UFilePath=oneTFolder+"/jsonData/jsonU/UData.json"
    with open (UFilePath,"r") as fptr:
        UData=json.load(fptr)

    UVec=UData["U"]
    LVec=np.array(LVec)

    UVec=np.array(UVec)

    LUProd=LVec*UVec

    LUMean=np.mean(LUProd)

    LMean=np.mean(LVec)

    UMean=np.mean(UVec)

    alphaVal=1/(TVal**2*LMean)*(LUMean-LMean*UMean)

    return alphaVal


tStart=datetime.now()

alphaAll=[]

for oneTFile in sortedTFiles:
    alphaTmp=compute_alpha(oneTFile)
    alphaAll.append(alphaTmp)
tEnd=datetime.now()
alphaAll=np.array(alphaAll)
print("alpha time: ",tEnd-tStart)
sortedTVals=np.array(sortedTVals)
TInds=np.where(sortedTVals<12)
TToPlt=sortedTVals[TInds]
interpolatedTVals=np.linspace(np.min(TToPlt)*0.9,np.max(TToPlt)*1.1,30)
plt.figure()

plt.plot(interpolatedTVals,[0]*len(interpolatedTVals),color="black",label="theory")
plt.scatter(TToPlt,alphaAll[TInds],color="red",label="mc")
plt.title("Thermal expansion")
plt.xlabel("$T$")
plt.ylabel("$\\alpha$")
plt.ylim((-0.025,0.025))
plt.legend(loc="best")
pathData=TFolderRoot
plt.savefig(pathData+"/alpha.png")
plt.close()