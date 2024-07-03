import pickle
import numpy as np
from datetime import datetime
import statsmodels.api as sm
import sys
import re
import warnings
from scipy.stats import ks_2samp
import glob
from pathlib import Path
import os
#This script checks if  U values reach equilibrium and writes summary file of U

if (len(sys.argv)!=4):
    print("wrong number of arguments")
    exit()

funcName=sys.argv[1]
rowName=sys.argv[2]
TStr=sys.argv[3]

TFolder="../dataAll/"+funcName+"/"+rowName+"/"

#name of observable
obs_name="U"
dataNumRequired=2000
# #search directory
# TVals=[]
# TFileNames=[]
# for TFile in glob.glob(TFolder+"/T*"):
#     # print(TFile)
#     matchT=re.search(r"T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",TFile)
#     if matchT:
#         TFileNames.append(TFile)
#         TVals.append(float(matchT.group(1)))
#
#
# #sort T values
# sortedInds=np.argsort(TVals)
# sortedTVals=[TVals[ind] for ind in sortedInds]
# sortedTFiles=[TFileNames[ind] for ind in sortedInds]


# print(sortedTFiles)
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

    return startingFileInd, startingVecPosition









def checkDataFilesForOneT(oneTFolder):
    """

    :param oneTFolder: Txxx
    :return:
    """
    TRoot=oneTFolder
    sortedDataFilesToRead=sort_data_files_by_lpStart(TRoot)
    if len(sortedDataFilesToRead)==0:
        print("no data.")
        exit(10)
    parsedStartingFileInd,parsedStartingVecPosition=parseSummary(oneTFolder)
    if parsedStartingFileInd>0:
        startingFileInd=parsedStartingFileInd
    else:
        startingFileInd=int(len(sortedDataFilesToRead)*1/2)#we guess that the equilibrium starts within this data file
    print("ind="+str(startingFileInd))
    startingFileName=sortedDataFilesToRead[startingFileInd]
    startingVecPosition=0
    #read starting data file
    with open(startingFileName,"rb") as fptr:
        vec=np.array(pickle.load(fptr))
        lengthTmp=len(vec)
        if parsedStartingVecPosition>0:
            startingVecPosition=parsedStartingVecPosition
        else:
            startingVecPosition=int(lengthTmp/2)#we guess that the equilibrium starts at this position
        vecTruncated=vec[startingVecPosition:]

        print("startingVecPosition="+str(startingVecPosition))

    #read the rest of the data files
    for inFile in sortedDataFilesToRead[(startingFileInd+1):]:
        with open(inFile,"rb") as fptr:
            inVec=np.array(pickle.load(fptr))
            vecTruncated=np.r_[vecTruncated,inVec]

    NLags=int(np.ceil(len(vecTruncated)*3/4))#maximum lag value
    # print("NLags="+str(NLags))
    eps=1e-3
    lagVal=0
    same=False
    summaryFolder=TRoot+"/summary_"+obs_name+"/"
    Path(summaryFolder).mkdir(parents=True, exist_ok=True)

    summaryFile=summaryFolder+"/summaryFile_"+obs_name+".txt"

    #compute auto-correlations
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVec=sm.tsa.acf(vecTruncated,nlags=NLags)
    except Warning as w:
        same=True

    #all values are the same, exit with code 3
    if same==True:
        with open(summaryFile,"w+") as fptr:
            msg="error: same\n"
            fptr.writelines(msg)
            return

    if np.min(np.abs(acfOfVec))<eps:
        lagVal=np.where(np.abs(acfOfVec)<=eps)[0][0]
        vecValsSelected=vecTruncated[::lagVal]
        lengthTmp=len(vecValsSelected)
        if lengthTmp%2==1:
            lengthTmp-=1
        vecValsToCompute=vecValsSelected[-lengthTmp:]
        lenPart=int(len(vecValsToCompute)/2)
        selectedFromPart0=vecValsToCompute[:lenPart]
        selectedFromPart1=vecValsToCompute[lenPart:]
        result = ks_2samp(selectedFromPart0, selectedFromPart1)
        numDataPoints=len(selectedFromPart0)+len(selectedFromPart1)
        if result.pvalue>0.1 and numDataPoints>=200:
            if numDataPoints>=dataNumRequired:
                newDataPointNum=0
            else:
                newDataPointNum=dataNumRequired-numDataPoints
            msg="equilibrium\n" \
                +"lag="+str(lagVal)+"\n" \
                +"K-S statistic: "+str(result.statistic)+"\n" \
                +"P-value: "+str(result.pvalue)+"\n" \
                +"numDataPoints="+str(numDataPoints)+"\n" \
                +"startingFileInd="+str(startingFileInd)+"\n" \
                +"startingVecPosition="+str(startingVecPosition)+"\n" \
                +"newDataPointNum="+str(newDataPointNum)+"\n"
            with open(summaryFile,"w+") as fptr:
                fptr.writelines(msg)
            return


        continueMsg="continue\n"
        if result.pvalue<=0.1:
            #not the same distribution
            continueMsg+="p value: "+str(result.pvalue)+"\n"
        if numDataPoints<200:
            #not enough data number
            continueMsg+="numDataPoints="+str(numDataPoints)+" too low\n"
        with open(summaryFile,"w+") as fptr:
            fptr.writelines(continueMsg)
        return


    else:
        msg="high correlation: "+str(np.min(np.abs(acfOfVec)))
        with open(summaryFile,"w+") as fptr:
            fptr.writelines(msg)
        return


# print(sortedTFiles)
# print("checking "+str(obs_name))
# for k in range(0,len(sortedTFiles)):
#     tStart=datetime.now()
#     oneTFolder=sortedTFiles[k]
#     # print(oneTFolder)
#     checkDataFilesForOneT(oneTFolder)
#     tEnd=datetime.now()
#     print("processed T="+str(sortedTVals[k])+": ",tEnd-tStart)

tStart=datetime.now()
oneTFolder=TFolder+"T"+TStr+"/"
checkDataFilesForOneT(oneTFolder)
tEnd=datetime.now()
print(obs_name+", processed T="+str(TStr)+": ",tEnd-tStart)