import pickle
import numpy as np
from datetime import datetime
import sys
import re
import glob
import os
import json
from pathlib import Path

#this script extracts effective data for L, y0,z0,y1
#they should be loaded together

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()

funcName=sys.argv[1]
rowName=sys.argv[2]

TFolderRoot="../dataAll/"+funcName+"/"+rowName+"/"

#name of observable
obs_name0="L"
obs_name1="y0"
obs_name2="z0"
obs_name3="y1"
obs_name4="U"
obsNamesAll=[obs_name0,obs_name1,obs_name2,obs_name3,obs_name4]

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


def parseSummary(oneTFolder,obs_name):
    """

    :param oneTFolder:
    :param obs_name:
    :return:
    """
    startingFileInd=-1
    startingVecPosition=-1
    lag=-1
    smrFile=oneTFolder+"/summary_"+obs_name+"/summaryFile_"+obs_name+".txt"
    summaryFileExists=os.path.isfile(smrFile)
    if summaryFileExists==False:
        return startingFileInd,startingVecPosition,-1

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

def sort_data_files_by_lpStart(oneTFolder,obs_name):
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


def data2jsonForOneT(oneTFolder,oneTStr,obs_name,startingFileInd,startingVecPosition,lag):
    """

    :param oneTFolder:
    :param oneTStr:
    :param startingFileInd:
    :param startingVecPosition:
    :return:
    """
    TRoot=oneTFolder
    sortedDataFilesToRead=sort_data_files_by_lpStart(TRoot,obs_name)
    print("ind="+str(startingFileInd))
    print("startingVecPosition="+str(startingVecPosition))
    startingFileName=sortedDataFilesToRead[startingFileInd]
    #read starting data file
    with open(startingFileName,"rb") as fptr:
        vec=np.array(pickle.load(fptr))
        vecTruncated=vec[startingVecPosition:]
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

    startingFileIndAll=[]
    startingVecPositionAll=[]
    lagAll=[]
    for oneObs in obsNamesAll:
        flIndTmp,vecIndTmp,lagTmp=parseSummary(oneTFolder,oneObs)
        if flIndTmp<0:
            print("summary file does not exist for "+oneTStr+" "+oneObs)
            continue
        startingFileIndAll.append(flIndTmp)
        startingVecPositionAll.append(vecIndTmp)
        lagAll.append(lagTmp)
    maxStartingFileInd=np.max(startingFileIndAll)
    maxStartingVecPosition=np.max(startingVecPositionAll)
    lag=np.max(lagAll)

    for oneObs in obsNamesAll:
        data2jsonForOneT(oneTFolder,oneTStr,oneObs,maxStartingFileInd,maxStartingVecPosition,lag)


    tEnd=datetime.now()
    print("processed T="+str(sortedTVals[k])+": ",tEnd-tStart)