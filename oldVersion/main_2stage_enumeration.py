#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2019: Zhicheng Zhu. All right reserved.

#Info:
# enumeration of two-stage, no solver
#
#Last update: 02/17/2019
#!/usr/bin/python

from __future__ import print_function
import sys
import itertools
import time
import class_info 


#######################################
#1. initialization, START FROM HERE!!!.
#######################################
#init system parameter
def main(sysInfo):

	start_time = time.clock();

	# 2.1 solve the last stage problem:
	omega = [];
	for i in itertools.product(list(range(sysInfo.comInfoAll[0].nStates)), \
					repeat = sysInfo.nComponents):
		omega.append(list(i));

	objValue = float("inf");
	solutionX = [];
	solutionY = [];
	solutionZ = [];
	secondStageObj = [];
	
	#get second-stage objective values.
	for w in range(len(omega)):
		scenState = omega[w];
		tmpObj = 0;
		z = 0;
		for i in range(sysInfo.nComponents):
			y = 0;
			if scenState[i] >= sysInfo.nStates - 1:
				y = 1;
				z = 1;
			tmpObj += sysInfo.comInfoAll[i].cCM * y;
		tmpObj = tmpObj + sysInfo.cS * z;
		secondStageObj.append(tmpObj);	
	
	#get solution space:
	for i in itertools.product([0, 1], repeat = sysInfo.nComponents):
		solutionSpace.append(list(i));
	
	#find out the optimal solution
	for i in range(len(solutionSpace)):
		soluitonSpaceI = solutionSpace[i];
		objTmp = 0;
		for w2 in range(len(omega)):
			objTmp1 = 1;
			infeasibile = False;
			if secondStageObj[w2] == 0:
				continue;
			for j in range(sysInfo.nComponents):
			#check feasibility first
				solIX = soluitonSpaceI[j]; 
				if (solIX == 0) and \
					(sysInfo.comInfoAll[j].initState >= \
					sysInfo.comInfoAll[j].nStates - 1):
					infeasibile = True;
					break;
				prob1 = sysInfo.comInfoAll[j].currentToFail;
				prob2 = sysInfo.comInfoAll[j].newToFail;
				objTmp1 = objTmp1 * (prob1 * (1 - solIX) + prob2 * solIX);
				if objTmp1 == 0:
					break;
			if infeasibile == True:
				break;
			objTmp = objTmp + objTmp1*secondStageObj[w2];			
		if infeasibile == True:
			objTmp = float("inf");	
		else:
			#add first stage
			solutionZTmp = 0;
			solutionYTmp = [];
			for ii in range(sysInfo.nComponents):
				objTmp = objTmp + sysInfo.comInfoAll[ii].cPM * soluitonSpaceI[ii];
				if soluitonSpaceI[ii] == 1:
					solutionZTmp = 1;
				if sysInfo.comInfoAll[ii].initState >= \
					sysInfo.comInfoAll[ii].nStates - 1:
					solutionYTmp.append(1);
					objTmp = objTmp + sysInfo.comInfoAll[ii].cCM - \
							sysInfo.comInfoAll[ii].cPM;
				else:
					solutionYTmp.append(0);
			objTmp = objTmp + solutionZTmp * sysInfo.cS;
			
		if objTmp < objValue:
			objValue = objTmp;
			solutionX = soluitonSpaceI;
			solutionZ = solutionZTmp;
			solutionY = solutionYTmp;				
				
	
		
	########################################
	#3. Result handling
	########################################
	end_time = time.clock();
	time_elapsed = end_time - start_time;
	N0 = [];
	N1 = [];
	for i in range(sysInfo.nComponents):
		if solutionX[i] == 0:
			N0.append(i);
		else:
			N1.append(i);
	sysInfo.time = time_elapsed;
	sysInfo.N0 = N0;
	sysInfo.N1 = N1;
	
	


