#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2019: Zhicheng Zhu. All right reserved.

#Info:
#script to control other main file for comparison
#
#Last update: 02/05/2018
#!/usr/bin/python


import main_2stage_only
import main_dynamic_solver
import main_myAlg
import class_info
import copy
import random
# some fixed parameters
nComV = [12];
for nComponents in nComV:
	#nComponents = 5;		#fix component numbers	
	nStates = 10;					#number of states for components, 0 - (m-1);

	cS = 50;				#fix setup cost
	inspInterval = 5;
	cInsp = 1;
	nStages = 2;			### do not change this!!
	#init component parameter
	#gamma distribution is assumed.
	gam_a =   	[1]*nComponents;
	gam_b =   	[2]*nComponents;
	S = 		[60]*nComponents;	#failure threshold
	cPM = [5]*nComponents;

	idxFrom = 79;
	idxTo = 80;
	idxInt = idxTo - idxFrom
	
	
	corError = 0;
	corErrorIdx = [];
	objError = 0;
	objErrorIdx = [];
	
	iterInfoLen = {};
	iterInfoMax = {};
	timeAll_1 = 0;
	timeAll_2 = 0;
	maxTime_1 = 0;
	maxTime_2 = 0;
	for idx in range(idxFrom, idxTo):			#control replicates, random seed
		#init system parameter
		sysInfo = class_info.system_info(nComponents, nStages, inspInterval, cS, cInsp);
		initState = [];
		cCM = [];
		for i in range(nComponents):
			baseSeedState = 100000;
			baseSeedCM = baseSeedState * 10;
			#for initial state
			ranSeed = (i+1)*baseSeedState + idx;			#random seeds!!!!!!!!!!!!!remember to change it
			random.seed(ranSeed);
			tmp = random.randint(0, nStates-1); #0<= state <= nStates - 1
			initState.append(tmp);
			#for CM
			ranSeed = (i+1)*baseSeedCM + idx;				#random seeds!!!!!!!!!!!!!remember to change it
			random.seed(ranSeed);
			tmp = random.randint(20, 100); 
			cCM.append(tmp);
		
		for i in range(nComponents):
			comInfo = class_info.component_info(i, gam_a[i], gam_b[i], nStates,\
									S[i], initState[i], inspInterval, cCM[i], cPM[i], cS);
			sysInfo.add_com(comInfo);

		
		#sysInfo1 = copy.deepcopy(sysInfo);
		#main_2stage_only.main(sysInfo1);

		sysInfo2 = copy.deepcopy(sysInfo);
		main_myAlg.main(sysInfo2, nComponents);	
		
		timeAll_2 += sysInfo2.time;
		maxTime_2 = max(maxTime_2, sysInfo2.time);
		if len(sysInfo2.iterInfo) in iterInfoLen:
			iterInfoLen[len(sysInfo2.iterInfo)] += 1;
		else:
			iterInfoLen[len(sysInfo2.iterInfo)] = 1;

		if max(sysInfo2.iterInfo) in iterInfoMax:
			iterInfoMax[max(sysInfo2.iterInfo)] += 1;
		else:
			iterInfoMax[max(sysInfo2.iterInfo)] = 1;
	print ("++++++++++++++");
	print ("overall stats");
	print ("++++++++++++++");
	print ("From To");
	print (	idxFrom, idxTo);
	print ("Number of replicates:");
	print (idxInt);
	
	print ("average time:");
	print (timeAll_2/idxInt)
	print ("max time:");
	print (maxTime_2);
	#averge converge j
	print ("AVERAGE number of j");
	for key in iterInfoLen:
		print (str(key)+":"+str(iterInfoLen[key]));
	print ("MAX number of j");
	for key in iterInfoMax:
		print (str(key)+":"+str(iterInfoMax[key]));


			
		'''
		#main_dynamic_solver.main(sysInfo1);
		## result collection
		#isN0Subset = (sysInfo2.N0 <= sysInfo1.N0);
		#isN1Subset = (sysInfo2.N1 <= sysInfo1.N1);
		isN0Equal = (set(sysInfo2.N0) == set(sysInfo1.N0));
		isN1Equal = (set(sysInfo2.N1) == set(sysInfo1.N1));
		isObjEqual = (round(sysInfo2.objValue,2) == round(sysInfo1.objValue,2));
		timeAll_1 += sysInfo1.time;
		timeAll_2 += sysInfo2.time;
		maxTime_1 = max(maxTime_1, sysInfo1.time);
		maxTime_2 = max(maxTime_2, sysInfo2.time);

		#print (idx)
		#print (sysInfo2.iterInfo)
		if len(sysInfo2.iterInfo) in iterInfoLen:
			iterInfoLen[len(sysInfo2.iterInfo)] += 1;
		else:
			iterInfoLen[len(sysInfo2.iterInfo)] = 1;

		if max(sysInfo2.iterInfo) in iterInfoMax:
			iterInfoMax[max(sysInfo2.iterInfo)] += 1;
		else:
			iterInfoMax[max(sysInfo2.iterInfo)] = 1;

		if isN0Equal == False or isN1Equal == False:
			corError += 1;
			corErrorIdx.append(idx);
		if isObjEqual == False:
			objError += 1;
			objErrorIdx.append(idx);#neither 0 or 1

			
		if isN0Equal == False or isN1Equal == False or\
			isObjEqual == False or len(sysInfo2.Nu) > 0:
		print (idx);
		print ("-----")
		
		print ("result comparison");
		print ("--------------------");
		print ("benchmark1 objValues:");
		print (sysInfo1.objValue);
		print ("benchmark objValues:");
		print (sysInfo2.objValue);
		print ("time info(benchmark,algorithm, dif):");
		print (round(sysInfo1.time, 1),\
			round(sysInfo2.time, 1),\
			round(sysInfo1.time - sysInfo2.time,1));
		print (" is N0 equal, is N1 equal");
		print ( set(sysInfo2.N0) == set(sysInfo1.N0),\
				set(sysInfo2.N1) == set(sysInfo1.N1));
		print ("N0=(Benchmark, algorithm)");
		print (sysInfo1.N0);
		print (sysInfo2.N0);
		print ("N1=(Benchmark, algorithm)");
		print (sysInfo1.N1);
		print (sysInfo2.N1);
		
		print ("Nu");
		print (sysInfo2.Nu);
		print ("===============================================");
	#overall output
	if (1):
		print ("++++++++++++++");
		print ("overall stats");
		print ("++++++++++++++");
		print ("From To");
		print (	idxFrom, idxTo);
		print ("Number of replicates:");
		print (idxInt);
		print ("corollary violation, rate and index:");
		print (corError, corError/idxInt);
		print (corErrorIdx);
		print ("objective violation, rate and index:");
		print (objError, objError/idxInt);
		print (objErrorIdx);
		print ("average time:");
		print (timeAll_1/idxInt, timeAll_2/idxInt)
		print ("max time:");
		print (maxTime_1, maxTime_2);
		#averge converge j
		print ("AVERAGE number of j");
		for key in iterInfoLen:
			print (str(key)+":"+str(iterInfoLen[key]));
		print ("MAX number of j");
		for key in iterInfoMax:
			print (str(key)+":"+str(iterInfoMax[key]));

	'''
	
	
	
	