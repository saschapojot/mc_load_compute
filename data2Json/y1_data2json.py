import pickle
import numpy as np
from datetime import datetime
import sys
import re
import glob
import os
import json
from pathlib import Path


#this script extracts effective data for y1

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()

funcName=sys.argv[1]
rowName=sys.argv[2]

TFolderRoot="../dataAll/"+funcName+"/"+rowName+"/"
#name of observable
obs_name="y1"
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

def sort_data_files_by_lpStart(oneTFolder):
    """

    :param oneTFolder: Txxx
    :return: pkl data files sorted by loopStart
    """
    dataFolderName=oneTFolder+"/data_files/"+obs_name+"_AllPickle/"
    # print(dataFolderName)
    dataFilesAll=[]
    loopStartAll=[]

    for oneDataFile in glob.glob(dataFolderName+"/*.pkl"):
        dataFilesAll.append(oneDataFile)
        matchStart=re.search(r"loopStart(\d+)",oneDataFile)
        if matchStart:
            loopStartAll.append(int(matchStart.group(1)))

    startInds=np.argsort(loopStartAll)
    # loopStartSorted=[loopStartAll[i] for i in startInds]
    sortedDataFiles=[dataFilesAll[i] for i in startInds]

    return sortedDataFiles


def parseSummary(oneTFolder):
    """

    :param oneTFolder:
    :return:
    """

    startingFileInd=-1
    startingVecPosition=-1
    lag=-1
    smrFile=oneTFolder+"/summary_"+obs_name+"/summaryFile_"+obs_name+".txt"
    summaryFileExists=os.path.isfile(smrFile)
    if summaryFileExists==False:
        return startingFileInd,startingVecPosition

    # eq=False
    with open(smrFile,"r") as fptr:
        lines=fptr.readlines()
    for oneLine in lines:
        # #match equilibrium
        # matchEq=re.search(r"equilibrium",oneLine)
        # if matchEq:
        #     eq=True

        #match startingFileInd
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))
        #match startingVecPosition
        matchStartingVecPosition=re.search(r"startingVecPosition=(\d+)",oneLine)
        if matchStartingVecPosition:
            startingVecPosition=int(matchStartingVecPosition.group(1))
        #match lag
        matchLag=re.search(r"lag=(\d+)",oneLine)
        if matchLag:
            lag=int(matchLag.group(1))

    return startingFileInd, startingVecPosition,lag


def data2jsonForOneT(oneTFolder,oneTStr):
    """

    :param oneTFolder: Txxx
    :return:
    """
    TRoot=oneTFolder
    sortedDataFilesToRead=sort_data_files_by_lpStart(TRoot)
    parsedStartingFileInd,parsedStartingVecPosition,lag=parseSummary(oneTFolder)
    if parsedStartingFileInd>=0:
        startingFileInd=parsedStartingFileInd
    if parsedStartingFileInd<0:
        print(oneTStr+": Equilibrium not known")
        return
    print("ind="+str(startingFileInd))
    startingFileName=sortedDataFilesToRead[startingFileInd]
    #read starting data file
    with open(startingFileName,"rb") as fptr:
        vec=np.array(pickle.load(fptr))

        startingVecPosition=parsedStartingVecPosition

        vecTruncated=vec[startingVecPosition:]

        print("startingVecPosition="+str(startingVecPosition))
    #read the rest of the data files
    for inFile in sortedDataFilesToRead[(startingFileInd+1):]:
        with open(inFile,"rb") as fptr:
            inVec=np.array(pickle.load(fptr))
            vecTruncated=np.r_[vecTruncated,inVec]

    vecValsSelected=vecTruncated[::lag]

    outJsonDataRoot=TFolderRoot+"/jsonOutAll/"
    outJsonFolder=outJsonDataRoot+"/"+oneTStr+"/jsonData/json"+obs_name+"/"
    Path(outJsonFolder).mkdir(parents=True, exist_ok=True)
    outJsonFile=outJsonFolder+"/"+obs_name+"Data.json"

    dataOut={obs_name:list(vecValsSelected)}
    with open(outJsonFile,"w+") as fptr:
        json.dump(dataOut,fptr,indent=4)
    print(obs_name+": "+oneTStr+", dataNum="+str(len(vecValsSelected)))

for k in range(0,len(sortedTFiles)):
    tStart=datetime.now()
    oneTFolder=sortedTFiles[k]
    oneTStr=sortedTStrings[k]
    data2jsonForOneT(oneTFolder,oneTStr)
    tEnd=datetime.now()
    print("processed T="+str(sortedTVals[k])+": ",tEnd-tStart)
