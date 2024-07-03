import os
import subprocess
import sys
import re
import numpy as np
from pathlib import Path
from decimal import Decimal
oneTCheckDir="./oneTCheckObservables/"

obsNames=["L","U","y0","y1","z0"]

def format_using_decimal(value):
    # Convert the float to a Decimal
    decimal_value = Decimal(value)
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit(4)
TStr=format_using_decimal(sys.argv[1])
print(TStr)
# TStr=sys.argv[1]

funcName="quadratic"
rowName="row0"

TFolder="./dataAll/"+funcName+"/"+rowName+"/T"+TStr+"/"

def readSummary(obs_name):
    """

    :param obs_name:
    :return:
    """
    smrFile=TFolder+"/summary_"+obs_name+"/summaryFile_"+obs_name+".txt"

    summaryFileExists=os.path.isfile(smrFile)
    eq=False
    startingFileInd=-1
    startingVecPosition=-1
    newDataPointNum=-1
    lagVal=-1
    if summaryFileExists==False:
        return eq, startingFileInd,startingVecPosition,newDataPointNum,lagVal
    with open(smrFile,"r") as fptr:
        lines=fptr.readlines()

    for oneLine in lines:
        matchEq=re.search(r"equilibrium",oneLine)
        if matchEq:
            eq=True
        matchLag=re.search(r"lag=(\d+)",oneLine)
        if matchLag:
            lagVal=int(matchLag.group(1))
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))

        matchStartingVecPosition=re.search(r"startingVecPosition=(\d+)",oneLine)
        if matchStartingVecPosition:
            startingVecPosition=int(matchStartingVecPosition.group(1))

        matchNewDataPointNum=re.search(r"newDataPointNum=(\d+)",oneLine)
        if matchNewDataPointNum:
            newDataPointNum=int(matchNewDataPointNum.group(1))

    return eq, startingFileInd,startingVecPosition,newDataPointNum,lagVal



active=True
execName="./run_mc_load_and_compute"
params=[TStr,"U"]
counter=0
while active:
    try:
        resultcpp=subprocess.run([execName]+params,capture_output=True, text=True, check=True)

        print("Output from C++ executable:")
        print(resultcpp.stdout)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        print(f"Return code: {e.returncode}")
        print(f"Output: {e.output}")

    os.chdir(oneTCheckDir)
    for obs_name in obsNames:
        scriptName="check"+obs_name+"OneT.py"
        subprocess.run(["python3", scriptName] + [funcName,rowName,TStr], capture_output=True, text=True, check=True)


    os.chdir("../")

    eqAll=[]
    startingFileIndAll=[]
    startingVecPositionAll=[]
    newDataPointNumAll=[]
    lagAll=[]
    for obs_name in obsNames:
        eqTmp, startingFileIndTmp,startingVecPositionTmp,newDataPointNumTmp,lagTmp=readSummary(obs_name)
        eqAll.append(eqTmp)
        startingFileIndAll.append(startingFileIndTmp)
        startingVecPositionAll.append(startingVecPositionTmp)
        newDataPointNumAll.append(newDataPointNumTmp)
        lagAll.append(lagTmp)

    AllEq=True
    for eq in eqAll:
        AllEq=AllEq and eq
    if AllEq==False:
        index_of_first_false = next((index for index, value in enumerate(eqAll) if not value), -1)
        params=[TStr,obsNames[index_of_first_false]]
        print("continue")
        continue
    else:
        if np.max(newDataPointNumAll)==0:
            active=False
        else:
            lagComb=np.max(lagAll)
            startingFileIndComb=np.max(startingFileIndAll)
            startingVecPositionComb=np.max(startingVecPositionAll)
            newDataPointNumComb=np.max(newDataPointNumAll)

            new_obs="comb"
            outCombDir=TFolder+"/summary_"+new_obs+"/"
            Path(outCombDir).mkdir(exist_ok=True,parents=True)

            outCombFile=outCombDir+"/summaryFile_"+new_obs+".txt"

            lines2Comb=["equilibrium\n",
                        "lag="+str(lagComb)+"\n",
                        "startingFileInd="+str(startingFileIndComb)+"\n",
                        "startingVecPosition="+str(startingVecPositionComb)+"\n",
                        "newDataPointNum="+str(newDataPointNumComb)+"\n"]


            with open(outCombFile,"w+") as fptr:
                fptr.writelines(lines2Comb)
            params=[TStr,new_obs]

    counter+=1
    if counter>=5:
        break
