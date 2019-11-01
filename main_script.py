#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2019: Zhicheng Zhu. All right reserved.

#Info:
#script to control other main file for comparison
#
#Last update: 03/04/2019
#!/usr/bin/python


import main_2stage_only
#import main_dynamic_solver
import main_myAlg
#import main_2stage_only_1st
import class_info
import copy
import random
import numpy as np
# some fixed parameters
nComV = [1];
#nComV = [18,100,300,500,1000,3000,5000];

for nComponents in nComV:
	#nComponents = 5;		#fix component numbers	
	nStates = 11;					#number of states for components, 0 - (m-1);
	cInsp = 1;
	cS = 200;				#fix setup cost
	inspInterval = 30;		################################################changed!!!!!(1)	
	nStages = 2;			### do not change this!!
	#init component parameter
	#gamma distribution is assumed.
	
	## new parameters
	#alpha = ct^b
	gam_a = [];
	gam_b = [];
	gam_a = [1.0824,1];	#alpha
	#gam_a = []
	#now it a rate parameter, following gamma(kappa, lambda), lambda is a rate as well
	gam_b = [8.556,7.654];	#shape & rate
	S = 		[2]*nComponents;	#failure threshold, 13mm
	
	idxFrom = nComponents-1;
	idxTo = idxFrom+11;
	idxInt = idxTo - idxFrom;

	corError = 0;
	corErrorIdx = [];
	objError = 0;
	objErrorIdx = [];
	
	iterInfoLen = {};
	iterInfoMax = {};
	timeAll_1 = 0;
	timeAll_2 = 0;
	timeAll_3 = 0;
	solGap = 0;
	
	maxJ = 6;
	
	cost3 = [0]*maxJ;
	time3 = [0]*maxJ;
	maxTime3 = [0]*maxJ;#collect the maximum time of the third: algorithm 2
	
	cost2 = 0;
	timeAll2Rep = [];	#collect all info for time of the second: algorithm 1.
	
	#heuristic error for the third: algorithm 2.
	heu_all = [0]*maxJ;
	heu_error = [0]*maxJ;
	heu_0 = [0]*maxJ;	#no error best
	heu_1 = [0]*maxJ;	#no error best
	heu_err0 = [0]*maxJ;	#error best
	heu_err1 = [0]*maxJ;
	for idx in range(idxFrom, idxTo):			#control replicates, random seed
		#init system parameter
		sysInfo = class_info.system_info(nComponents, nStages, inspInterval, cS, cInsp);
		initState = [];
		cCM = [];
		cPM = [];
		initAge = [];
		for i in range(nComponents):
			baseSeedState = 10000+nComponents*100;
			baseSeedCM = baseSeedState * 10;
			baseSeedPM = baseSeedState * 100;
			baseSeedT = baseSeedState * 1000;
			baseSeedGa = baseSeedState * 2000;
			baseSeedGb = baseSeedState * 2000;

			#for initial state
			ranSeed = (i+1)*baseSeedState + idx;			#random seeds
			random.seed(ranSeed);
			tmp = random.randint(0, nStates-1); #0<= state <= nStates - 1
			initState.append(idx);
			#for CM
			ranSeed = (i+1)*baseSeedCM + idx;				#random seeds!!!!!!!!!!!!!remember to change it
			random.seed(ranSeed);
			tmp = random.uniform(10,30); 
			cCM.append(20);
			#for PM
			ranSeed = (i+1)*baseSeedPM + idx;				#random seeds!!!!!!!!!!!!!remember to change it
			random.seed(ranSeed);
			tmp = random.uniform(1, 5); 
			cPM.append(5);
			#for randomize initial age
			ranSeed = (i+1)*baseSeedT + idx;				#random seeds!!!!!!!!!!!!!remember to change it
			random.seed(ranSeed);
			tmp = random.uniform(1, 14); 
			initAge.append(tmp);
		
			#for randomize gamma alpha
			ranSeed = (i+1)*baseSeedGa + idx;				#random seeds!!!!!!!!!!!!!remember to change it
			random.seed(ranSeed);
			tmp = random.uniform(1, 5); 
			#gam_a.append(tmp);
			
			#for randomize gamma beta
			ranSeed = (i+1)*baseSeedGb + idx;				#random seeds!!!!!!!!!!!!!remember to change it
			random.seed(ranSeed);
			tmp = random.uniform(1, 5); 
			#gam_b.append(tmp);
		print (initState);	
		for i in range(nComponents):
			comInfo = class_info.component_info(i, gam_a, gam_b, nStates,\
									S[i], initState[i], initAge[i], inspInterval, cCM[i], cPM[i], cS);
			sysInfo.add_com(comInfo);

		
		#sysInfo1 = copy.deepcopy(sysInfo);
		#main_2stage_only.main(sysInfo1);

		#print (idx,nComponents);
		#print ("alg 1")
		sysInfo2 = copy.deepcopy(sysInfo);
		main_myAlg.main(sysInfo2, nComponents);	
		
		for i in sysInfo2.N1:
			print (i, initState[i]);
		
		
		'''
		for J in list(range(maxJ)):
			#print ("alg 2");
			#print (J);
			sysInfo3 = copy.deepcopy(sysInfo);
			#main_2stage_only_1st.main(sysInfo3);	
			main_myAlg.main(sysInfo3,J+1);
			cost3[J] += sysInfo3.objValue;
			time3[J] += sysInfo3.time;
			maxTime3[J] = max(maxTime3[J],sysInfo3.time);
			if max(sysInfo2.iterInfo) > J+1:
				heu_all[J] += 1;
				if round(sysInfo3.objValue,2) != round(sysInfo2.objValue,2):
					#error!
					print ("------------------");
					print ("heuristic error:");
					print (nComponents,idx);
					print ("obj compare(real,heuristic)");
					print (round(sysInfo2.objValue,2),round(sysInfo3.objValue,2));
					print ("N0 compare(real,heuristic)");
					print (sysInfo2.N0);
					print (sysInfo3.N0);
					print ("N1 compare(real,heuristic)");
					print (sysInfo2.N1);
					print (sysInfo3.N1);
					heu_error[J] += 1;
					if len(sysInfo3.Nu_0) == 0:
						heu_err1[J] += 1;
					elif len(sysInfo3.Nu_1) == 0: 
						heu_err0[J] += 1;
					else:
						print ("--------(in error)heuristic not 0 not 1:");	
						print (sysInfo3.Nu_0);
						print (sysInfo3.Nu_1);
				else:
					if len(sysInfo3.Nu_0) == 0:
						heu_1[J] += 1;
					elif len(sysInfo3.Nu_1) == 0: 
						heu_0[J] += 1;
					else:
						print ("++++(no error)heuristic not 0 not 1:");	
						print (sysInfo3.Nu_0);
						print (sysInfo3.Nu_1);					
		## result collection
		# for sysInfo1:
		if (0):
			timeAll_1 += sysInfo1.time;
	
		if (0):
			timeAll_2 += sysInfo2.time;
			timeAll2Rep.append(sysInfo2.time);
			cost2 += sysInfo2.objValue;		
			if len(sysInfo2.iterInfo) in iterInfoLen:
				iterInfoLen[len(sysInfo2.iterInfo)] += 1;
			else:
				iterInfoLen[len(sysInfo2.iterInfo)] = 1;
			print ("iterInfo");
			print (sysInfo2.iterInfo);
			if max(sysInfo2.iterInfo) in iterInfoMax:
				iterInfoMax[max(sysInfo2.iterInfo)] += 1;
			else:
				iterInfoMax[max(sysInfo2.iterInfo)] = 1;

		if (0):
			isN0Equal = (set(sysInfo2.N0) == set(sysInfo1.N0));
			isN1Equal = (set(sysInfo2.N1) == set(sysInfo1.N1));
			isObjEqual = (round(sysInfo2.objValue,2) == round(sysInfo1.objValue,2));

			if isN0Equal == False or isN1Equal == False:
				corError += 1;
				corErrorIdx.append(idx);
			if isObjEqual == False:
				objError += 1;
				objErrorIdx.append(idx);#neither 0 or 1

			if isN0Equal == False or isN1Equal == False or\
				isObjEqual == False or len(sysInfo2.Nu) > 0:
				print ("result comparison");
				print ("--------------------");
				print ("benchmark1 objValues:");
				print (sysInfo1.objValue);
				print ("algorithm objValues:");
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
	print ("++++++++++++++");
	print ("overall stats");
	print ("++++++++++++++");
	print ("From To");
	print (	idxFrom, idxTo);
	print ("Number of replicates:");
	print (idxInt);
	print ("number of components:");
	print (nComponents);

	if (0):	
		print ("-------");
		print ("CPLEX:");
		print ("-------");
		print ("avg time");
		print (round(timeAll_1/idxInt,4));
	if (0):	
		print ("-------");
		print ("Alg. 1:");
		print ("-------");
		print ("avg. time:");				
		print (round(timeAll_2/idxInt,4));
		tmp = 0;
		for key in iterInfoMax:
			tmp += key*iterInfoMax[key];
		print ("variance of time:");
		print (np.var(timeAll2Rep));
		print ("avg max j:");
		print (tmp/idxInt);
		print ("maximum time:")
		print (max(timeAll2Rep));
		print ("iterInfoMax");
		print (iterInfoMax);
		print ("max j:")
		print (max(iterInfoMax));
	
		
	if (0):
		print ("-------");
		print ("Alg. 2:");
		print ("-------");
		print ("average time:");
		for j in list(range(maxJ)):
			print ("j="+str(j+1)+":(time avg, time max, error%)");
			print ("-----");
			print (round(time3[j]/idxInt,4), maxTime3[j], abs(cost3[j]-cost2)/cost2);
			print ("heuristic all & error");
			print (heu_all[j], heu_error[j]);
			print ("error all 0 all 1");
			print (heu_err0[j], heu_err1[j]);
			print ("normal all 0 all 1");
			print (heu_0[j], heu_1[j]);
	if (0):
		print ("-----------");
		print ("Other info:");
		print ("-----------");
		print ("corollary violation, rate and index:");
		print (corError, corError/idxInt);
		print (corErrorIdx);
		print ("objective violation, rate and index:");
		print (objError, objError/idxInt);
		print (objErrorIdx);
'''
	