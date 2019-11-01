#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2019: Zhicheng Zhu. All right reserved.

#Info:
#main file to solve deterministic equivalent form of CBM model, two-stage only....
#
#Last update: 02/04/2019
#!/usr/bin/python

from __future__ import print_function
import sys
import cplex
import itertools
import time
import class_info

		
	
def main(sysInfo):
	start_time = time.clock();

	setS = [];	#start from j = 2
	for j in range(2, 2):	#sysInfo.nComponents >=2;
		setSj = [];
		for i in itertools.combinations(list(range(sysInfo.nComponents)), j):
			setSj.append(list(i));
		setS.append(setSj);
	


	coeA = [];	#1 * nComponents
	coeB = [];	#1 * nComponents
	coeCurrentToFail = [];
	coeNewToFail = [];
	for i in range(sysInfo.nComponents):
		tmp = sysInfo.comInfoAll[i].currentToFail;
		coeCurrentToFail.append(tmp);
		coeB.append(1-tmp);
		tmp1 = sysInfo.comInfoAll[i].newToFail;
		coeNewToFail.append(tmp1)
		coeA.append(tmp - tmp1);

	#get coefficient of X
	coeX = [];		
	for i in range(sysInfo.nComponents):
		tmpi = (coeNewToFail[i] - coeCurrentToFail[i])\
			   *sysInfo.comInfoAll[i].cCM;
		tmpr = 1;
		for r in range(sysInfo.nComponents):
			if r != i:
				tmpr = tmpr * coeB[r];
			else:
				tmpr = tmpr * coeA[r];
		tmpi = tmpi - tmpr*sysInfo.cS + sysInfo.comInfoAll[i].cPM;
		coeX.append(tmpi);

	#get coefficient of u
	coeU = [];
	if len(setS) > 0:
		for j in range(len(setS)):
			setSj = setS[j];
			tmpj = [];
			for k in range(len(setSj)):
				setSjk = setSj[k];
				tmpk = 1;
				for i in range(sysInfo.nComponents):
					if i in setSjk:
						tmpk = tmpk * coeA[i];
					else:
						tmpk = tmpk * coeB[i];
				tmpj.append(-tmpk*sysInfo.cS);
			coeU.append(tmpj);

	#constant term
	consTerm = 0;
	tmp1 = 1;
	tmp2 = 0;
	for i in range(sysInfo.nComponents):
		tmp1 = tmp1*coeB[i];
		tmp2 = coeCurrentToFail[i]*sysInfo.comInfoAll[i].cCM + tmp2;
	consTerm = tmp2 + sysInfo.cS - tmp1*sysInfo.cS;

	#init cplex
	cpx = cplex.Cplex();	
	cpx.objective.set_sense(cpx.objective.sense.minimize);
	cpx.set_log_stream(None);
	cpx.set_error_stream(None);
	cpx.set_warning_stream(None);
	cpx.set_results_stream(None);
	#decision variables:
	#x
	varNameX = [];
	varX = [];
	for i in range(sysInfo.nComponents):
		varNameX.append("x"+str(i));
		varX.append(cpx.variables.get_num());
		cpx.variables.add(obj = [coeX[i]],
						lb = [0.0], ub=[1.0], types=["B"],
						names=[varNameX[i]]);
	#y
	varNameY = [];
	varY = [];
	for i in range(sysInfo.nComponents):
		varNameY.append("y"+str(i));
		varY.append(cpx.variables.get_num());
		cpx.variables.add(obj = [sysInfo.comInfoAll[i].cCM - sysInfo.comInfoAll[i].cPM],
						lb = [0.0], ub=[1.0], types=["B"],
						names=[varNameY[i]]);	

	#z
	varNameZ = "z";
	varZ = [cpx.variables.get_num()];
	cpx.variables.add(obj = [sysInfo.cS],
			lb = [0.0], ub=[1.0], types=["B"],
			names=[varNameZ]);	

	#u
	varNameU = [];
	varU = [];
	for j in range(len(coeU)):
		for k in range(len(coeU[j])):
			varNameU.append("u" + str(j) + str(k));
			varU.append(cpx.variables.get_num());
			cpx.variables.add(obj = [coeU[j][k]],
						lb = [0.0], ub=[1.0], types=["B"],
						names=[varNameU[-1]]);	

	#constraints:
	#1
	for i in range(sysInfo.nComponents):
		cpx.linear_constraints.add(
			lin_expr=[cplex.SparsePair([varNameX[i], varNameZ], [1, -1])],
			senses=["L"], 
			range_values=[0.0],
			rhs=[0]);			

	#2
	for i in range(sysInfo.nComponents):
		comIFrom = sysInfo.comInfoAll[i].initState;
		cpx.linear_constraints.add(
			lin_expr=[cplex.SparsePair([varNameY[i]],[-comIFrom])],
			senses=["L"], 
			range_values=[0.0],
			rhs=[sysInfo.comInfoAll[i].nStates-2-comIFrom]);
		cpx.linear_constraints.add(
			lin_expr=[cplex.SparsePair([varNameY[i], varNameX[i]],[1, -1])],
			senses=["L"], 
			range_values=[0.0],
			rhs=[0.0]);
	
	#new constraints: for u
	idxTmp = -1;
	for j in range(len(setS)):
		for k in range(len(setS[j])):
			idxTmp += 1;
			nameVec = []
			coeVec = [];
			for i in setS[j][k]:
				nameVec.append(varNameX[i]);
				coeVec.append(1);
				cpx.linear_constraints.add(
					lin_expr=[cplex.SparsePair([varNameU[idxTmp], varNameX[i]],[1, -1])],
					senses=["L"], 
					range_values=[0.0],
					rhs=[0.0]);
			nameVec.append(varNameU[idxTmp]);
			coeVec.append(-1);
			cpx.linear_constraints.add(
					lin_expr=[cplex.SparsePair(nameVec, coeVec)],
					senses=["L"], 
					range_values=[0.0],
					rhs=[len(setS[j][k]) - 1]);	# j+2-1
	end_time = time.clock();
	time1 = end_time - start_time;
	start_time = time.clock();					
	cpx.solve();
	solution = cpx.solution;
	solutionAll = solution.get_values();
	#get X
	minTmp = varX[0];
	maxTmp = varX[-1] + 1;
	solutionXStage = solutionAll[minTmp:maxTmp];


	#get Y
	minTmp = varY[0];
	maxTmp = varY[-1] + 1;
	solutionYStage = solutionAll[minTmp:maxTmp];

	#get Z
	minTmp = varZ[0];
	maxTmp = varZ[-1] + 1;
	solutionZStage = solutionAll[minTmp:maxTmp];
	
	#get u
#	minTmp = varU[0];
#	maxTmp = varU[-1] + 1;
	#solutionUStage = solutionAll[minTmp:maxTmp];

	objValue = solution.get_objective_value() + consTerm;
	end_time = time.clock();
	sysInfo.time = end_time - start_time + time1;
	sysInfo.objValue = objValue;

	#print (solutionXStage);
	#print (solutionYStage);
	#print (solutionZStage);
	#print (solutionUStage);
	#print (objValue);
	#print (time1,sysInfo.time);
	
	N1 = [];
	N0 = [];


	for i in range(sysInfo.nComponents):
		if solutionXStage[i] == 0:
			N0.append(i);
		else:
			N1.append(i);
			
	sysInfo.N0 = N0;
	sysInfo.N1 = N1;




















