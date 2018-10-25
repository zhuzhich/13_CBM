#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2018: Zhicheng Zhu. All right reserved.

#Info:
#back-tracing + enumeration, no solver
#
#Last update: 10/14/2018
#!/usr/bin/python

from __future__ import print_function
import sys
import itertools
import time
from scipy.stats import gamma

#####################################
#class info
#####################################
class component_info():
	def transProb(self, stateFrom, stateTo, inspItvl):
		if stateFrom > stateTo:
			return 0;
		stepSize = self.failTsh/(self.nStates - 1); 	#step size for normal states 
		degFrom = stateFrom * stepSize;			#degradation lower bound of the state
		degToU = (stateTo + 1) * stepSize;		#degradation upper bound of the state
		degToL = stateTo * stepSize;			#degradation lower bound of the state

		if stateTo >= self.nStates - 1:
			deltaDeg = self.failTsh - degFrom;
			prob = 1 - gamma.cdf(deltaDeg, self.gammaAlpha*inspItvl, scale=self.gammaBeta);
		else:
			deltaDeg1 = degToU - degFrom;
			prob1 = gamma.cdf(deltaDeg1, self.gammaAlpha*inspItvl, scale=self.gammaBeta);
			deltaDeg2 = degToL - degFrom;
			prob2 = gamma.cdf(deltaDeg2, self.gammaAlpha*inspItvl, scale=self.gammaBeta);
			prob = prob1 - prob2;
		return prob;
		
		
	'''	
	def state2lv():
		crtState = self.crtState;
		bound = [];
		bound.append(0);#put it here for now..
		bound.append(1);
		return bound;
	'''
	def __init__(self, idx, gam_a, gam_b, states, S, \
					initState,cCM, cPM):
		self.idx = idx;
		self.gammaAlpha = gam_a;
		self.gammaBeta = gam_b;
		self.nStates = states;		# 0 ... nStates - 1. nStates - 1 is failure states.
		self.failTsh = S;	#failure threshold
		self.initState = initState;
		#self.crtState = initState;
		#self.crtDgLvRange = self.state2lv();
		self.cCM = cCM;
		self.cPM = cPM;
		
#system information
#parameters
class system_info():
	def add_com(self, comInfo):
		self.comInfoAll.append(comInfo);
	def __init__(self, N, T, inspInterval, cS, cInsp):
		self.nComponents = N;
		self.nStages = T;
		self.inspItvl = inspInterval;
		self.cS = cS;
		self.cInsp = cInsp;
		self.comInfoAll = [];
		
	

#######################################
#1. initialization, START FROM HERE!!!.
#######################################
#init system parameter
nComponents = 2;
nStages = 3;
inspInterval = 10;
cS = 20;				#setup cost
cInsp = 1;
sysInfo = system_info(nComponents, nStages, inspInterval, cS, cInsp);

#init component parameter
#gamma distribution is assumed.
nStates = 3;					#number of states for components, 0 - (m-1);
gam_a =   	[1]*nComponents;
gam_b =   	[5]*nComponents;
S = 		[60]*nComponents;	#failure threshold
initState = [1, 2];
cCM = [20]*nComponents;
cPM = [5]*nComponents;

for i in range(nComponents):
	comInfo = component_info(i, gam_a[i], gam_b[i], nStates,\
							S[i], initState[i], cCM[i], cPM[i]);
	sysInfo.add_com(comInfo);


########################################
#2. build two-stage DEF model and run 
########################################
start_time = time.clock();

# 2.1 solve the last stage problem:
omega = [];
for i in itertools.product(list(range(nStates)), repeat = sysInfo.nComponents):
	omega.append(list(i));

lastStageObj = [];		
lastStageX = [];
lastStageY = [];
lastStageZ = [];		
for w in range(len(omega)):
	scenState = omega[w];
	tmpObj = 0;
	z = 0;
	lastStageXw = [];
	lastStageYw = [];
	for i in range(sysInfo.nComponents):
		y = 0;
		if scenState[i] == nStates - 1:
			y = 1;
			z = 1;
		x = y;
		lastStageXw.append(x);
		lastStageYw.append(y);
		tmpObj = tmpObj + sysInfo.comInfoAll[i].cPM * x + \
				(sysInfo.comInfoAll[i].cCM - sysInfo.comInfoAll[i].cPM) * y;
	tmpObj = tmpObj + sysInfo.cS * z;
	lastStageZ.append(z);
	lastStageX.append(lastStageXw);
	lastStageY.append(lastStageYw);
	lastStageObj.append(tmpObj);

objValues = [];
solutionX = [];
solutionY = [];
solutionZ = [];
objValues.append(lastStageObj);
solutionX.append(lastStageX);
solutionY.append(lastStageY);
solutionZ.append(lastStageZ);


#2.2 start from nStages - 1 (nStages-2 in python), calculate the obj value of two stage problem.
for stageIdx in range(nStages-2 , -1, -1):	#nStages-2 , nStages-3, ..., 1, 0
	#for each stage, we have m**n unrepeated nodes.
	objValuesStage = [];
	solutionXStage = [];
	solutionYStage = [];
	solutionZStage = [];
	
	for w1 in range(len(omega)):		#equivalent to first stage. 
		breakW1 = False;
		if stageIdx > 0:
			comStatesFrom = omega[w1];
		else:
			comStatesFrom = [];
			breakW1 = True;
			for i in range(sysInfo.nComponents):
				comStatesFrom.append(sysInfo.comInfoAll[i].initState);
		#
		
		objValuesScen = float("inf");
		solutionXScen = [];
		solutionYScen = [];
		solutionZScen = 0;		
		
		solutionSpace = [];
		childNodesObj = objValues[-1];

		#get solution space:
		for i in itertools.product([0, 1], repeat = sysInfo.nComponents):
			solutionSpace.append(list(i));
		
		#find out the optimal solution
		for i in range(len(solutionSpace)):
			soluitonSpaceI = solutionSpace[i];
			objTmp = 0;
			for w2 in range(len(omega)):
				comStatesTo = omega[w2];
				objTmp1 = 1;
				infeasibile = False;
				if childNodesObj[w2] == 0:
					continue;
				for j in range(sysInfo.nComponents):
				#check feasibility first
					comIFrom = comStatesFrom[j];
					comITo = comStatesTo[j];
					solIX = soluitonSpaceI[j]; 
					
					if (solIX == 0) and (comIFrom == sysInfo.comInfoAll[j].nStates - 1):
						infeasibile = True;
						break;
					prob1 = sysInfo.comInfoAll[j].transProb(comIFrom, comITo, sysInfo.inspItvl);
					prob2 = sysInfo.comInfoAll[j].transProb(0, comITo, sysInfo.inspItvl);
					objTmp1 = objTmp1 * (prob1 * (1 - solIX) + prob2 * solIX);
					if objTmp1 == 0:
						break;
				if infeasibile == True:
					break;
				objTmp = objTmp + objTmp1*childNodesObj[w2];			
			if infeasibile == True:
				objTmp = float("inf");	

			#add first stage
			solutionZScenTmp = 0;
			solutionYScenTmp = [];
			for ii in range(sysInfo.nComponents):
				objTmp = objTmp + sysInfo.comInfoAll[ii].cPM * soluitonSpaceI[ii];
				if soluitonSpaceI[ii] == 1:
					solutionZScenTmp = 1;
				if comStatesFrom[ii] == sysInfo.comInfoAll[ii].nStates - 1:
					solutionYScenTmp.append(1);
					objTmp = objTmp + sysInfo.comInfoAll[ii].cCM - sysInfo.comInfoAll[ii].cPM;
				else:
					solutionYScenTmp.append(0);
			objTmp = objTmp + solutionZScenTmp * sysInfo.cS;
				
			if objTmp < objValuesScen:
				objValuesScen = objTmp;
				solutionXScen = soluitonSpaceI;
				solutionZScen = solutionZScenTmp;
				solutionYScen = solutionYScenTmp;				
				
		objValuesStage.append(objValuesScen);
		solutionXStage.append(solutionXScen);
		solutionYStage.append(solutionYScen);
		solutionZStage.append(solutionZScen);
		if breakW1 == True:
			break;
	########for stages
	objValues.append(objValuesStage);
	solutionX.append(solutionXStage);
	solutionY.append(solutionYStage);
	solutionZ.append(solutionZStage);

	
########################################
#3. Result handling
########################################
end_time = time.clock();

time_elapsed = end_time - start_time;


f = open("log2.txt", "w");
old = sys.stdout;
sys.stdout = f;

print ("\n===============================main_2stage_enumeration, (m, n, t)=(%d,%d,%d)============" 
		%(nStates, sysInfo.nComponents, nStages));

print ("calculation time is %f"  %time_elapsed);

for stageIdx in range(nStages):
	curStage = nStages - stageIdx - 1;
	for w1 in range(len(objValues[stageIdx])):
		print ("=======(stage, scen) = (%d, %d)========" %(curStage,w1));
		print ("objValues:");
		print (objValues[stageIdx][w1]);
		print ("solutionX:");
		print (solutionX[stageIdx][w1]);
		print ("solutionY:");
		print (solutionY[stageIdx][w1]);
		print ("solutionZ:");
		print (solutionZ[stageIdx][w1]);
		print ("===================\n");

## 4. end of file 
sys.stdout = old;
f.close();				
		