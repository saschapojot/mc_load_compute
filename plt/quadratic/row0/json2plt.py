import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd


#This script loads json data and plot

# if (len(sys.argv)!=3):
#     print("wrong number of arguments")
#     exit()

funcName="quadratic"
rowName="row0"
jsonFolderRoot="../../../dataAll/"+funcName+"/"+rowName+"/jsonOutAll/"

inDataFile="../../../inputData/"+funcName+"/quadraticCoeffs.txt"

line=""

with open(inDataFile,"r") as fptr:
    lines=fptr.readlines()
for oneLine in lines:
    matchRow=re.search(rowName,oneLine)
    if matchRow:
        line=oneLine
        break

#match consts
#match a1
match_a1=re.search(r"a1=([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
if match_a1:
    a1=float(match_a1.group(1))

#match a2
match_a2=re.search(r"a2=([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
if match_a2:
    a2=float(match_a2.group(1))


#match c1
match_c1=re.search(r"c1=([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
if match_c1:
    c1=float(match_c1.group(1))

#match c2
match_c2=re.search(r"c2=([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
if match_c2:
    c2=float(match_c2.group(1))

TVals=[]
TFileNames=[]

for TFile in glob.glob(jsonFolderRoot+"/T*"):

    matchT=re.search(r"T(\d+(\.\d+)?)",TFile)
    # if float(matchT.group(1))<1:
    #     continue

    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))


sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]



def pltU(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: U plots, U mean, U var
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    UFilePath=oneTFile+"/jsonData/jsonU/UData.json"
    with open(UFilePath, 'r') as fptr:
        data = json.load(fptr)
    UVec=np.array(data["U"])
    print("T="+str(TVal)+", data num="+str(len(UVec)))
    meanU=np.mean(UVec)

    varU=np.var(UVec,ddof=1)
    sigmaU=np.sqrt(varU)

    nbins=500
    fig=plt.figure()
    axU=fig.add_subplot()
    (n0,_,_)=axU.hist(UVec,bins=nbins)
    meanU=np.round(meanU,4)
    sigmaU=np.round(sigmaU,4)

    axU.set_title("T="+str(TVal))
    axU.set_xlabel("$U$")
    axU.set_ylabel("#")
    xPosUText=(np.max(UVec)-np.min(UVec))*1/2+np.min(UVec)
    yPosUText=np.max(n0)*2/3
    axU.text(xPosUText,yPosUText,"mean="+str(meanU)+"\nsd="+str(sigmaU))
    plt.axvline(x=meanU,color="red",label="mean")
    axU.text(meanU*1.1,0.5*np.max(n0),str(meanU)+"$\pm$"+str(sigmaU),color="red")
    axU.hlines(y=0,xmin=meanU-sigmaU,xmax=meanU+sigmaU,color="green",linewidth=15)

    plt.legend(loc="best")

    EHistOut="T"+str(TVal)+"UHist.png"
    plt.savefig(oneTFile+"/"+EHistOut)

    plt.close()

    ### test normal distribution for mean U
    #block mean
    USelectedAll=UVec
    def meanPerBlock(length):
        blockNum=int(np.floor(len(USelectedAll)/length))
        UMeanBlock=[]
        for blkNum in range(0,blockNum):
            blkU=USelectedAll[blkNum*length:(blkNum+1)*length]
            UMeanBlock.append(np.mean(blkU))
        return UMeanBlock
    fig=plt.figure(figsize=(20,20))
    fig.tight_layout(pad=5.0)
    lengthVals=[2,5,7,10]
    for i in range(0,len(lengthVals)):
        l=lengthVals[i]
        UMeanBlk=meanPerBlock(l)
        ax=fig.add_subplot(2,2,i+1)
        (n,_,_)=ax.hist(UMeanBlk,bins=100,color="aqua")
        xPosTextBlk=(np.max(UMeanBlk)-np.min(UMeanBlk))*1/7+np.min(UMeanBlk)
        yPosTextBlk=np.max(n)*3/4
        meanTmp=np.mean(UMeanBlk)
        meanTmp=np.round(meanTmp,3)
        sdTmp=np.sqrt(np.var(UMeanBlk))
        sdTmp=np.round(sdTmp,3)
        ax.set_title("L="+str(l))
        ax.text(xPosTextBlk,yPosTextBlk,"mean="+str(meanTmp)+", sd="+str(sdTmp))
    fig.suptitle("T="+str(TVal))
    plt.savefig(oneTFile+"/T"+str(TVal)+"UBlk.png")
    # plt.savefig(EBlkMeanDir+"/T"+str(TTmp)+"EBlk.png")
    plt.close()
    return [meanU,varU]


def plt_dist(oneTFile):
    """

    :param oneTFile: corresponds to one temperature

    :return: y1Mean,y1Var,LMean,LVar,d0A0BMean,d1A1BMean
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))

    LFilePath=oneTFile+"/jsonData/jsonL/LData.json"
    with open(LFilePath,"r") as fptr:
        LVec=np.array(json.load(fptr)["L"])


    y0FilePath=oneTFile+"/jsonData/jsony0/y0Data.json"
    with open(y0FilePath,"r") as fptr:
        y0Vec=np.array(json.load(fptr)["y0"])


    z0FilePath=oneTFile+"/jsonData/jsonz0/z0Data.json"
    with open(z0FilePath,"r") as fptr:
        z0Vec=np.array(json.load(fptr)["z0"])


    y1FilePath=oneTFile+"/jsonData/jsony1/y1Data.json"
    with open(y1FilePath,"r") as fptr:
        y1Vec=np.array(json.load(fptr)["y1"])
    # print("LVec="+str(LVec))
    LMean=np.mean(LVec)
    LVar=np.var(LVec,ddof=1)

    # print("len(LVec)="+str(len(LVec)))
    # print("len(y0Vec)="+str(len(y0Vec)))
    # print("len(z0Vec)="+str(len(z0Vec)))
    # print("len(y1Vec)="+str(len(y1Vec)))
    lengths=[len(LVec),len(y0Vec),len(z0Vec),len(y1Vec)]
    minLen=np.min(lengths)
    LVecSelected=LVec[-minLen:]
    y0VecSelected=y0Vec[-minLen:]
    z0VecSelected=z0Vec[-minLen:]
    y1VecSelected=y1Vec[-minLen:]

    d0A0BVec=[]#y0
    d0B1AVec=[]#z0
    d1A1BVec=[]#y1
    d1B0AVec=[]#z1

    for i in range(0,len(y1VecSelected)):
        LTmp=LVecSelected[i]
        y0Tmp=y0VecSelected[i]
        z0Tmp=z0VecSelected[i]
        y1Tmp=y1VecSelected[i]
        d0A0BVec.append(y0Tmp)
        d0B1AVec.append(z0Tmp)
        d1A1BVec.append(y1Tmp)
        d1B0AVec.append(LTmp-y0Tmp-z0Tmp-y1Tmp)


    z1VecSelected=np.array(d1B0AVec)
    # #summary of distance between neighboring points
    d0A0BMean=np.mean(d0A0BVec)
    d0B1AMean=np.mean(d0B1AVec)
    d1A1BMean=np.mean(d1A1BVec)
    d1B0AMean=np.mean(d1B0AVec)

    d0A0BVar=np.var(d0A0BVec,ddof=1)
    d0B1AVar=np.var(d0B1AVec,ddof=1)
    d1A1BVar=np.var(d1A1BVec,ddof=1)
    d1B0AVar=np.var(d1B0AVec,ddof=1)

    d0A0BSd=np.sqrt(d0A0BVar)
    d0B1ASd=np.sqrt(d0B1AVar)
    d1A1BSd=np.sqrt(d1A1BVar)
    d1B0ASd=np.sqrt(d1B0AVar)

    y1Mean=np.mean(y1VecSelected)
    y1Var=np.var(y1VecSelected,ddof=1)

    smrDistFileName=oneTFile+"/distSummary.txt"
    #
    dist=[d0A0BMean,d0B1AMean,d1A1BMean,d1B0AMean]
    varsAll=[d0A0BVar,d0B1AVar,d1A1BVar,d1B0AVar]
    sdAll=[d0A0BSd,d0B1ASd,d1A1BSd,d1B0ASd]
    fptrTxt=open(smrDistFileName,"w")
    fptrTxt.write("T="+str(TVal)+"\n")
    fptrTxt.write("distances: "+str(dist)+"\n")
    fptrTxt.write("var: "+str(varsAll)+"\n")
    fptrTxt.write("sd: " +str(sdAll)+"\n")
    fptrTxt.close()

    rowNames=["0A0B","0B1A","1A1B","1B0A"]
    dfOut=pd.DataFrame({"dist":dist,"var":varsAll},index=rowNames)

    dfOut.to_csv(oneTFile+"/distVar.csv")

    #cov functions
    y0z0Mean=np.mean(y0VecSelected*z0VecSelected)
    y0y1Mean=np.mean(y0VecSelected*y1VecSelected)
    z0z1Mean=np.mean(z0VecSelected*z1VecSelected)
    y0LMean=np.mean(y0VecSelected*LVecSelected)
    y1LMean=np.mean(y1VecSelected*LVecSelected)

    prodData={"y0z0Mean":y0z0Mean,"y0y1Mean":y0y1Mean,"z0z1Mean":z0z1Mean,"y0LMean":y0LMean
        ,"y1LMean":y1LMean}
    prodFile=oneTFile+"/prod.json"
    with open(prodFile, 'w') as json_file:
        json.dump(prodData, json_file, indent=4)

    return [y1Mean,y1Var,LMean,LVar,d0A0BMean,d0B1AMean,d1A1BMean,d1B0AMean,d0A0BVar,d0B1AVar,d1A1BVar,d1B0AVar]


UMeanValsAll=[]
UVarValsAll=[]
y1MeanValsAll=[]
y1VarValsAll=[]
LMeanValsAll=[]
LVarValsAll=[]
d0A0BMeanValsAll=[]
d0B1AMeanValsAll=[]
d1A1BMeanValsAll=[]
d1B0AMeanValsAll=[]

d0A0BVarValsAll=[]
d0B1AVarValsAll=[]
d1A1BVarValsAll=[]
d1B0AVarValsAll=[]

tStatsStart=datetime.now()
# TToPlt=[]
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    # if sortedTVals[k]>60:
    #     continue
    # TToPlt.append(sortedTVals[k])
    UMeanTmp,UVarTmp=pltU(oneTFile)
    UMeanValsAll.append(UMeanTmp)
    UVarValsAll.append(UVarTmp)
    y1Mean,y1Var,LMean,LVar,d0A0BMean,d0B1AMean,d1A1BMean,d1B0AMean,d0A0BVar,d0B1AVar,d1A1BVar,d1B0AVar=plt_dist(oneTFile)
    y1MeanValsAll.append(y1Mean)
    y1VarValsAll.append(y1Var)
    LMeanValsAll.append(LMean)
    LVarValsAll.append(LVar)
    d0A0BMeanValsAll.append(d0A0BMean)
    d0B1AMeanValsAll.append(d0B1AMean)
    d1A1BMeanValsAll.append(d1A1BMean)
    d1B0AMeanValsAll.append(d1B0AMean)

    d0A0BVarValsAll.append(d0A0BVar)
    d0B1AVarValsAll.append(d0B1AVar)
    d1A1BVarValsAll.append(d1A1BVar)
    d1B0AVarValsAll.append(d1B0AVar)


tStatsEnd=datetime.now()
print("stats total time: ",tStatsEnd-tStatsStart)

UMeanValsAll=np.array(UMeanValsAll)
UVarValsAll=np.array(UVarValsAll)
y1MeanValsAll=np.array(y1MeanValsAll)
y1VarValsAll=np.array(y1VarValsAll)
LMeanValsAll=np.array(LMeanValsAll)
LVarValsAll=np.array(LVarValsAll)
d0A0BMeanValsAll=np.array(d0A0BMeanValsAll)
d0B1AMeanValsAll=np.array(d0B1AMeanValsAll)
d1A1BMeanValsAll=np.array(d1A1BMeanValsAll)
d1B0AMeanValsAll=np.array(d1B0AMeanValsAll)




sortedTVals=np.array(sortedTVals)
print(sortedTVals)
TInds=np.where(sortedTVals<12)
TToPlt=sortedTVals[TInds]
interpolatedTVals=np.linspace(np.min(TToPlt)*0.9,np.max(TToPlt)*1.1,30)


def EV(T):
    '''

    :param T: temperature
    :return: asymptotic value of E(V)
    '''
    return 2*T
plt.figure()
EVVals=[EV(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,EVVals,color="green",label="theory")
plt.scatter(TToPlt,UMeanValsAll[TInds],color="darkred",label="mc")
plt.title("E(V)")
plt.xlabel("$T$")
plt.legend(loc="best")
plt.savefig(jsonFolderRoot+"/EV.png")
plt.close()


def varV(T):
    """

    :param T: temperature
    :return:  asymptotic value of var(V)
    """
    return 2*T**2
plt.figure()
plt.scatter(TToPlt,UVarValsAll[TInds],color="violet",label="mc")
varVVals=[varV(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,varVVals,color="navy",label="theory")
plt.title("var(V)")
plt.xlabel("$T$")
plt.legend(loc="best")
plt.savefig(jsonFolderRoot+"/varV.png")
plt.close()

plt.figure()
plt.scatter(TToPlt,d0A0BMeanValsAll[TInds],color="deeppink",label="d0A0B")
plt.scatter(TToPlt,d1A1BMeanValsAll[TInds],color="aqua",label="d1A1B",s=1)
plt.title("Intracell distance")
plt.xlabel("$T$")
plt.ylabel("distance")
plt.ylim((0.5,1.5))
plt.legend(loc="best")
plt.xscale("log")
plt.savefig(jsonFolderRoot+"/intracellDistance.png")
plt.close()



#compare d0A0B and d1A1B
#d0A0B
d0A0BVarValsAll=np.array(d0A0BVarValsAll)
d0A0BSdValsAll=np.sqrt(d0A0BVarValsAll)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
ax1.errorbar(TToPlt, d0A0BMeanValsAll[TInds], yerr=d0A0BSdValsAll[TInds], fmt='o',color="black", ecolor='r', capsize=5, label='d0A0B')
ax1.set_xlabel('$T$')
ax1.set_ylabel('E(d0A0B)')
ax1.set_title('Intracell distance between 0A and 0B')
ax1.legend()



d1A1BVarValsAll=np.array(d1A1BVarValsAll)
d1A1BSdValsAll=np.sqrt(d1A1BVarValsAll)
ax2.errorbar(TToPlt, d1A1BMeanValsAll[TInds], yerr=d1A1BSdValsAll[TInds], fmt='o',color="black", ecolor='b', capsize=5, label='d1A1B')
ax2.set_xlabel('$T$')
ax2.set_ylabel('E(d1A1B)')
ax2.set_title('Intracell distance between 1A and 1B')
ax2.legend()
plt.tight_layout()
plt.savefig(jsonFolderRoot+"/d0A0Bd1A1B.pdf")
plt.close()

#compare d0B1A and d1B0A
#d0B1A
d0B1AVarValsAll=np.array(d0B1AVarValsAll)
d0B1ASdValsAll=np.sqrt(d0B1AVarValsAll)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
ax1.errorbar(TToPlt, d0B1AMeanValsAll[TInds], yerr=d0B1ASdValsAll[TInds], fmt='o',color="black", ecolor='r', capsize=5, label='d0B1A')
ax1.set_xlabel('$T$')
ax1.set_ylabel('E(0B1A)')
ax1.set_title('Intracell distance between 0B and 1A')
ax1.legend()
#d1B0A
d1B0AVarValsAll=np.array(d1B0AVarValsAll)
d1B0ASdValsAll=np.sqrt(d1B0AVarValsAll)
ax2.errorbar(TToPlt, d1B0AMeanValsAll[TInds], yerr=d1B0ASdValsAll[TInds], fmt='o',color="black", ecolor='b', capsize=5, label='d1B0A')
ax2.set_xlabel('$T$')
ax2.set_ylabel('E(1B0A)')
ax2.set_title('Intracell distance between 1B and 0A')
ax2.legend()
plt.tight_layout()
plt.savefig(jsonFolderRoot+"/d0B1Ad1B0A.pdf")
plt.close()


def EL(T):
    """

    :param T: temperature
    :return: asymptotic value of E(L)
    """
    return 2*a1+2*a2

def varL(T):
    """

    :param T: temperature
    :return: asymptotic value of var(L)
    """

    return  c1**(-1)*c2**(-1)*(c1+c2)*T


plt.figure()
plt.scatter(TToPlt,LMeanValsAll[TInds],color="red",label="mc")
ELVals=[EL(T) for T in interpolatedTVals]
ELVals=np.array(ELVals)
plt.plot(interpolatedTVals,ELVals,color="blue",label="theory")
plt.title("E(L)")
plt.ylabel("E(L)")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(jsonFolderRoot+"/EL.png")
plt.close()

plt.figure()
plt.scatter(TToPlt,LVarValsAll[TInds],color="magenta",label="mc")

varLVals=[varL(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,varLVals,color="green",label="theory")
plt.title("var(L)")
plt.ylabel("var(L)")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(jsonFolderRoot+"/varL.png")
plt.close()

#L, bar plot

LSdValsAll=np.sqrt(LVarValsAll)
fig, ax = plt.subplots()
plt.plot(interpolatedTVals,ELVals,color="blue",label="theory")

ax.errorbar(TToPlt,LMeanValsAll[TInds],yerr=LSdValsAll[TInds],fmt='o',color="red", ecolor='aqua', capsize=5, label='mc')
ax.set_xlabel('$T$')
ax.set_ylabel('E(L)')
ax.set_title('E(L) with error bar')
ax.legend()
plt.legend(loc="best")
plt.savefig(jsonFolderRoot+"/ELErrBar.png")

plt.close()
def Ey1(T):
    """

    :param T: temperature
    :return: asymptotic value of E(y1)
    """

    return a1



plt.figure()
plt.scatter(TToPlt,y1MeanValsAll[TInds],color="red",label="mc")
Ey1Vals=[Ey1(T) for T in interpolatedTVals]

plt.plot(interpolatedTVals,Ey1Vals,color="blue",label="theory")

plt.title("E(y1)")
plt.ylabel("E(y1)")
plt.xlabel("$T$")
plt.legend(loc="best")

# plt.ylim((0,1.2))

plt.savefig(jsonFolderRoot+"/Ey1.png")
plt.close()


def vary1(T):
    """

    :param T: temperature
    :return: asymptotic value of var(y1)
    """


    val=1/4*c1**(-1)*(c1+c2)**(-1)*(2*c1+c2)*T \
        +1/4*c1**(-1)*c2*(c1+c2)**(-1)*T
    return val



plt.figure()
plt.scatter(TToPlt,y1VarValsAll[TInds],color="magenta",label="mc")
vary1Vals=[vary1(T) for T in interpolatedTVals]
plt.plot(interpolatedTVals,vary1Vals,color="green",label="theory")

plt.title("var(y1)")
plt.ylabel("var(y1)")
plt.xlabel("$T$")
plt.legend(loc="best")

# plt.ylim((0,1.5))
plt.savefig(jsonFolderRoot+"/vary1.png")
plt.close()


#E(y1) with error bar
y1SdValsAll=np.sqrt(y1VarValsAll)
fig, ax = plt.subplots()
plt.plot(interpolatedTVals,Ey1Vals,color="blue",label="theory")
ax.errorbar(TToPlt,y1MeanValsAll[TInds],yerr=y1SdValsAll[TInds],fmt='o',color="red", ecolor='aqua', capsize=5, label='cyan')
ax.set_xlabel('$T$')
ax.set_ylabel('$E(y_{1})$')
ax.set_title('$E(y_{1})$ with error bar')
ax.legend()
plt.legend(loc="best")
plt.savefig(jsonFolderRoot+"/Ey1ErrBar.png")
#combine all distVar.csv files

combined_d0A0BMeanValsAll=[]
combined_d0B1AMeanValsAll=[]
combined_d1A1BMeanValsAll=[]
combined_d1B0AMeanValsAll=[]

combined_d0A0BVarValsAll=[]
combined_d0B1AVarValsAll=[]
combined_d1A1BVarValsAll=[]
combined_d1B0AVarValsAll=[]

combinedTValsAll=[]

for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    distFile=oneTFile+"/distVar.csv"
    dfDist=pd.read_csv(distFile,header=0)
    combinedTValsAll.append(sortedTVals[k])


    distTmp=dfDist.loc[:,"dist"]
    combined_d0A0BMeanValsAll.append(distTmp[0])
    combined_d0B1AMeanValsAll.append(distTmp[1])
    combined_d1A1BMeanValsAll.append(distTmp[2])
    combined_d1B0AMeanValsAll.append(distTmp[3])

    varTmp=dfDist.loc[:,"var"]
    combined_d0A0BVarValsAll.append(varTmp[0])
    combined_d0B1AVarValsAll.append(varTmp[1])
    combined_d1A1BVarValsAll.append(varTmp[2])
    combined_d1B0AVarValsAll.append(varTmp[3])


combined_d0A0BMeanValsAll=np.array(combined_d0A0BMeanValsAll)
combined_d0B1AMeanValsAll=np.array(combined_d0B1AMeanValsAll)
combined_d1A1BMeanValsAll=np.array(combined_d1A1BMeanValsAll)
combined_d1B0AMeanValsAll=np.array(combined_d1B0AMeanValsAll)


combined_d0A0BVarValsAll=np.array(combined_d0A0BVarValsAll)
combined_d0B1AVarValsAll=np.array(combined_d0B1AVarValsAll)
combined_d1A1BVarValsAll=np.array(combined_d1A1BVarValsAll)
combined_d1B0AVarValsAll=np.array(combined_d1B0AVarValsAll)

combined_d0A0BSdValsAll=np.sqrt(combined_d0A0BVarValsAll)
combined_d0B1ASdValsAll=np.sqrt(combined_d0B1AVarValsAll)
combined_d1A1BSdValsAll=np.sqrt(combined_d1A1BVarValsAll)
combined_d1B0ASdValsAll=np.sqrt(combined_d1B0AVarValsAll)

df0A0BCombinedOut=pd.DataFrame({"T":combinedTValsAll,"0A0Bdist":combined_d0A0BMeanValsAll,"sd":combined_d0A0BSdValsAll,"var":combined_d0A0BVarValsAll})
df0A0BCombinedOut.to_csv(jsonFolderRoot+"/0A0B.csv",index=False)


df0B1ACombinedOut=pd.DataFrame({"T":combinedTValsAll,"0B1Adist":combined_d0B1AMeanValsAll,"sd":combined_d0B1ASdValsAll,"var":combined_d0B1AVarValsAll})
df0B1ACombinedOut.to_csv(jsonFolderRoot+"/0B1A.csv",index=False)

df1A1BCombinedOut=pd.DataFrame({"T":combinedTValsAll,"1A1Bdist":combined_d1A1BMeanValsAll,"sd":combined_d1A1BSdValsAll,"var":combined_d1A1BVarValsAll})
df1A1BCombinedOut.to_csv(jsonFolderRoot+"/1A1B.csv",index=False)


df1B0ACombinedOut=pd.DataFrame({"T":combinedTValsAll,"1B0Adist":combined_d1B0AMeanValsAll,"sd":combined_d1B0ASdValsAll,"var":combined_d1B0AVarValsAll})
df1B0ACombinedOut.to_csv(jsonFolderRoot+"/1B0A.csv",index=False)



#correlation functions

y0z0MeanValsAll=[]
y0y1MeanValsAll=[]
z0z1MeanValsAll=[]
y0LMeanValsAll=[]
y1LMeanValsAll=[]
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    prodFile=oneTFile+"/prod.json"
    with open(prodFile,"r") as fptr:
        prodData=json.load(fptr)

    y0z0Tmp=prodData["y0z0Mean"]
    y0y1Tmp=prodData["y0y1Mean"]
    z0z1Tmp=prodData["z0z1Mean"]
    y0LTmp=prodData["y0LMean"]
    y1LTmp=prodData["y1LMean"]

    y0z0MeanValsAll.append(y0z0Tmp)
    y0y1MeanValsAll.append(y0y1Tmp)
    z0z1MeanValsAll.append(z0z1Tmp)
    y0LMeanValsAll.append(y0LTmp)
    y1LMeanValsAll.append(y1LTmp)

y0z0MeanValsAll=np.array(y0z0MeanValsAll)
y0y1MeanValsAll=np.array(y0y1MeanValsAll)
z0z1MeanValsAll=np.array(z0z1MeanValsAll)
y0LMeanValsAll=np.array(y0LMeanValsAll)
y1LMeanValsAll=np.array(y1LMeanValsAll)


y0MeanValsAll=combined_d0A0BMeanValsAll
z0MeanValsAll=combined_d0B1AMeanValsAll
y1MeanValsAll=combined_d1A1BMeanValsAll
z1MeanValsAll=combined_d1B0AMeanValsAll


#cov(y0,z0)
cov_y0z0=y0z0MeanValsAll-y0MeanValsAll*z0MeanValsAll

plt.figure()
plt.plot(interpolatedTVals,[0]*len(interpolatedTVals),color="green",label="theory")

plt.scatter(TToPlt,cov_y0z0[TInds],color="red",label="mc")
plt.title("cov($y_{0}$,$z_{0}$)")
plt.ylabel("cov")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(jsonFolderRoot+"/covy0z0.png")
plt.close()


#cov(y0,y1)
cov_y0y1=y0y1MeanValsAll-y0MeanValsAll*y1MeanValsAll
plt.figure()
plt.plot(interpolatedTVals,[0]*len(interpolatedTVals),color="green",label="theory")

plt.scatter(TToPlt,cov_y0y1[TInds],color="red",label="mc")
plt.title("cov($y_{0}$,$y_{1}$)")
plt.ylabel("cov")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(jsonFolderRoot+"/covy0y1.png")
plt.close()



#cov(z0,z1)
cov_z0z1=z0z1MeanValsAll-z0MeanValsAll*z1MeanValsAll
plt.figure()

plt.scatter(TToPlt,cov_z0z1[TInds],color="red",label="mc")
plt.plot(interpolatedTVals,[0]*len(interpolatedTVals),color="blue",label="theory")
plt.title("cov($z_{0}$,$z_{1}$)")
plt.ylabel("cov")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(jsonFolderRoot+"/covz0z1.png")
plt.close()


#cov(y0,L)
LMeanValsAll=np.array(LMeanValsAll)
cov_y0L=y0LMeanValsAll-y0MeanValsAll*LMeanValsAll
plt.figure()
plt.plot(interpolatedTVals,vary1Vals,color="green",label="theory")
plt.scatter(TToPlt,cov_y0L[TInds],color="red",label="mc")
plt.title("cov(y0,L)")
plt.ylabel("cov")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(jsonFolderRoot+"/covy0L.png")
plt.close()


#cov(y1,L)
cov_y1L=y1LMeanValsAll-y1MeanValsAll*LMeanValsAll
plt.figure()
plt.plot(interpolatedTVals,vary1Vals,color="blue",label="theory")

plt.scatter(TToPlt,cov_y1L[TInds],color="magenta",label="mc")
plt.title("cov(y1,L)")
plt.ylabel("cov")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(jsonFolderRoot+"/covy1L.png")
plt.close()