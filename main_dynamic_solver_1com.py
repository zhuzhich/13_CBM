#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2018: Zhicheng Zhu. All right reserved.

#Info:
#main file to solve deterministic equivalent form of CBM model by back-tracking
#
#Last update: 10/08/2018
#!/usr/bin/python

from __future__ import print_function
import sys
import cplex
import itertools
import time
from scipy.stats import gamma


def main(sysInfo):

	start_time = time.clock();

	# 2.1 solve the last stage problem:
	omega = [];
	for i in itertools.product(list(range(sysInfo.comInfoAll[0].nStates)), repeat = sysInfo.nComponents):
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
			if scenState[i] == sysInfo.comInfoAll[0].nStates - 1:
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
	for stageIdx in range(sysInfo.nStages-2 , -1, -1):	#nStages-2 , nStages-3, ..., 1, 0
		#for each stage, we have m**n unrepeated nodes.
		objValuesStage = [];
		solutionXStage = [];
		solutionYStage = [];
		solutionZStage = [];
		
		for w1 in range(len(omega)):		#equivalent to first stage. 
			cpx = cplex.Cplex();	
			cpx.objective.set_sense(cpx.objective.sense.minimize);
			cpx.set_log_stream(None);
			cpx.set_error_stream(None);
			cpx.set_warning_stream(None);
			cpx.set_results_stream(None);
			breakW1 = False;

			if stageIdx > 0:
				comStatesFrom = omega[w1];
			else:
				breakW1 = True;
				comStatesFrom = [];
				for i in range(sysInfo.nComponents):
					comStatesFrom.append(sysInfo.comInfoAll[i].initState);

			#get subsets which cardinality = j
			setS = [];	#start from j = 2
			for j in range(2, sysInfo.nComponents + 1):	#sysInfo.nComponents >=2;
				setSj = [];
				for i in itertools.combinations(list(range(sysInfo.nComponents)), j):
					setSj.append(list(i));
				setS.append(setSj);
			
			#get coefficient of a and b
			coeA = [];	#n*scen
			coeB = [];	#n*scen
			for i in range(sysInfo.nComponents):
				biw = [];
				aiw = [];
				for w2 in range(len(omega)):
					comStatesTo = omega[w2];
					comIFrom = comStatesFrom[i];
					comITo = comStatesTo[i];
					tmp = sysInfo.comInfoAll[i].transProb(comIFrom, comITo, sysInfo.inspItvl);
					biw.append(tmp);
					aiw.append(sysInfo.comInfoAll[i].transProb(0, comITo, sysInfo.inspItvl) - tmp);
				coeA.append(aiw);
				coeB.append(biw);
			#if stageIdx == 1 and w1 == 0:
			#	print ("!!!!coeB = " + str(coeB));
				
			#get coefficient of X
			coeX = [];		
			childNodesObj = objValues[-1];
			for i in range(sysInfo.nComponents):
				tmpi = 0;
				for w2 in range(len(omega)):
					tmpr = coeA[i][w2];
					for r in range(sysInfo.nComponents):
						if r != i:
							tmpr = tmpr * coeB[r][w2];
					tmpi = tmpi + childNodesObj[w2]*tmpr;
				tmpi = tmpi + sysInfo.comInfoAll[i].cPM;
				coeX.append(tmpi);
				
			#get coefficient of u
			coeU = [];
			childNodesObj = objValues[-1];
			if len(setS) > 0:
				for j in range(len(setS)):
					setSj = setS[j];
					tmpj = [];
					for k in range(len(setSj)):
						setSjk = setSj[k];
						tmpk = 0;
						for w2 in range(len(omega)):
							tmpw = childNodesObj[w2];
							if tmpw == 0:
								continue;
							for i in range(sysInfo.nComponents):
								if i in setSjk:
									tmpw = tmpw * coeA[i][w2];
								else:
									tmpw = tmpw * coeB[i][w2];
							tmpk = tmpk + tmpw;
						tmpj.append(tmpk);
					coeU.append(tmpj);
			
			#constant term
			consTerm = 0;
			for w2 in range(len(omega)):
				tmp1 = 1;
				for i in range(sysInfo.nComponents):
					tmp1 = tmp1*coeB[i][w2];
				#if stageIdx == 1 and w1 == 0:
					#print ("!!!!tmp1 = " + str(tmp1));
				consTerm = consTerm + tmp1 * childNodesObj[w2];
			
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
				#print (varNameY[i]);
				#print (sysInfo.comInfoAll[i].cCM - sysInfo.comInfoAll[i].cPM);
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
			#original constraints:
			#1
			for i in range(sysInfo.nComponents):
				cpx.linear_constraints.add(
					lin_expr=[cplex.SparsePair([varNameX[i], varNameZ], [1, -1])],
					senses=["L"], 
					range_values=[0.0],
					rhs=[0]);			


				
			#2
			for i in range(sysInfo.nComponents):
				comIFrom = comStatesFrom[i];
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
												
			cpx.solve();
			solution = cpx.solution;
			
			#get objective values
			objValuesStage.append(solution.get_objective_value() + consTerm);
			
			solutionAll = solution.get_values();

			#get X
			minTmp = varX[0];
			maxTmp = varX[-1] + 1;
			solutionXStage.append(solutionAll[minTmp:maxTmp]);
			
			#get Y
			minTmp = varY[0];
			maxTmp = varY[-1] + 1;
			solutionYStage.append(solutionAll[minTmp:maxTmp]);
			
			#get Z
			minTmp = varZ[0];
			maxTmp = varZ[-1] + 1;
			solutionZStage.append(solutionAll[minTmp:maxTmp]);
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
	sysInfo.time = time_elapsed;

	#f = open("log1.txt", "w");
	#old = sys.stdout;
	#sys.stdout = f;

	#print ("\n===============================main_2stage_solver, (m, n, t)=(%d,%d,%d)============" 
	#		%(sysInfo.comInfoAll[0].nStates, sysInfo.nComponents, sysInfo.nStages));

	#print ("calculation time is %f"  %time_elapsed);
	
	print ("init state");
	print (sysInfo.comInfoAll[0].initState);
	for stageIdx in range(sysInfo.nStages):
		curStage = sysInfo.nStages - stageIdx - 1;
		for w1 in range(len(objValues[curStage])):
			if solutionX[curStage][w1][0] == 1:
				print ("=======(stage, scen) = (%d, %d)========" %(stageIdx,w1));
				break;
			#print ("component state");
			#print (w1%sysInfo.comInfoAll[0].nStates);
			
			#print ("objValues:");
			#print (objValues[-1][0]);
			#sysInfo.objValue = objValues[-1][0];
			#print ("solutionX:");
			#print (solutionX[curStage][w1]);
			#print ("solutionY:");
			#print (solutionY[curStage][w1]);
			#print ("solutionZ:");
			#print (solutionZ[curStage][w1]);
	N1 = [];
	N0 = [];

	solutionXStage = solutionX[-1][0];
	for i in range(sysInfo.nComponents):
		if solutionXStage[i] == 0:
			N0.append(i);
			sysInfo.comInfoAll[i].x = 0;
		else:
			N1.append(i);
			sysInfo.comInfoAll[i].x = 1;
	
	
	sysInfo.N0 = set(N0);
	sysInfo.N1 = set(N1);
	
	
	#print ("===================\n");

	#sys.stdout = old;
	#f.close();		
		
		