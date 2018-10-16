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
nComponents = 5;
nStages = 2;
inspInterval = 10;
cS = 20;				#setup cost
cInsp = 1;
sysInfo = system_info(nComponents, nStages, inspInterval, cS, cInsp);

#init component parameter
#gamma distribution is assumed.
nStates = 4;					#number of states for components, 0 - (m-1);
gam_a =   	[1, 1, 1, 1, 1];
gam_b =   	[5, 5, 5, 5, 5, 5];
S = 		[60, 60, 60, 60, 60];	#failure threshold
initState = [0, 0, 0, 0, 0];
cCM = [20, 20, 20, 20, 20];
cPM = [5, 5, 5, 5, 5];

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
	lastStageZw = [];	
	for i in range(sysInfo.nComponents):
		y = 0;
		lastStageYw.append(y);
		if scenState[i] == nStates - 1:
			y = 1;
			z = 1;
		x = y;
		lastStageXw.append(x);
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
		childNodesObj = objValues[nStages-2];

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
					prob = sysInfo.comInfoAll[j].transProb(comIFrom, comITo, sysInfo.inspItvl);
					objTmp1 = objTmp1 * prob;
				if infeasibile == True:
					break;
				objTmp = objTmp + objTmp1*childNodesObj[w2];			
			if infeasibile == True:
				objTmp = float("inf");	
			if objTmp < objValuesScen:
				objValuesScen = objTmp;
				solutionXScen = soluitonSpaceI;
		
		#find the optimal obj value
		for i in range(sysInfo.nComponents):
			objValuesScen = objValuesScen + sysInfo.comInfoAll[i].cPM * solutionXScen[i];
			if solutionXScen[i] == 1:
				solutionZScen = 1;
			if comStatesFrom[i] == sysInfo.comInfoAll[i].nStates - 1:
				solutionYScen.append(1);
				objValuesScen = objValuesScen + sysInfo.comInfoAll[i].cCM - sysInfo.comInfoAll[i].cPM;
			else:
				solutionYScen.append(0);
		objValuesScen = objValuesScen + solutionZScen * sysInfo.cS;
		
		objValuesStage = objValuesScen;
		solutionXStage = solutionXScen;
		solutionYStage = solutionYScen;
		solutionZStage = solutionZScen;
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

		
		
		