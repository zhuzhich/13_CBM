#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2019: Zhicheng Zhu. All right reserved.

#Info:
#script to control other main file for rolling horizion
#
#Last update: 03/09/2019
#!/usr/bin/python


import main_dynamic_solver
import main_dynamic_solver_1com
import main_myAlg
import class_info
import copy
import random
import time
from scipy.stats import f
from scipy.stats import gamma
import numpy as np

def get_objValue(sysInfo2, res):
	cost = 0;
	for t in range(len(res)):
		x = res[t][0];
		y = res[t][1];
		z = res[t][2];
		for i in range(len(x)):
			#print (sysInfo2.comInfoAll[i].cPM,sysInfo2.comInfoAll[i].cCM)
			cost += sysInfo2.comInfoAll[i].cPM*x[i];
			cost += (sysInfo2.comInfoAll[i].cCM - sysInfo2.comInfoAll[i].cPM)*y[i];
		cost += sysInfo2.comInfoAll[i].cS*z;
	return cost;

def rollHorizon(sysInfo2,rep,t):
	# rep&t is to control the random seed
	for i in range(sysInfo2.nComponents):
		if i in sysInfo2.N1:
			sysInfo2.comInfoAll[i].initState = 0;
		
		ranSeed = sysInfo2.nComponents*i + sysInfo2.nStages*t + rep;
		random.seed(ranSeed);
		prob = random.uniform(0,1);
		if (1):
			alpha_t = sysInfo2.comInfoAll[i].gammaAlpha[0]*\
						sysInfo2.inspItvl;
			beta0 = sysInfo2.comInfoAll[i].gammaBeta[0];
			beta1 = sysInfo2.comInfoAll[i].gammaBeta[1];
			tmp = f.ppf(prob, 2*alpha_t, 2*beta0);
			delDeg = tmp*alpha_t/(beta0*beta1);
			stepSize = sysInfo2.comInfoAll[i].failTsh\
						/(sysInfo2.comInfoAll[i].nStates - 1); 	#step size for normal states 
			degInc = int((delDeg+sysInfo2.comInfoAll[i].deg)/stepSize);
			sysInfo2.comInfoAll[i].deg = delDeg+sysInfo2.comInfoAll[i].deg - degInc*stepSize;
			sysInfo2.comInfoAll[i].initState = \
				min(sysInfo2.comInfoAll[i].nStates - 1, \
					sysInfo2.comInfoAll[i].initState + degInc);
		if (0):
			alpha_t = sysInfo2.comInfoAll[i].gammaAlpha*\
						sysInfo2.inspItvl;
			scale_p = sysInfo2.comInfoAll[i].gammaBeta			
			delDeg = gamma.ppf(prob,alpha_t,scale=scale_p);
			#print (delDeg);
			stepSize = sysInfo2.comInfoAll[i].failTsh\
						/(sysInfo2.comInfoAll[i].nStates - 1); 	#step size for normal states 
			degInc = int((delDeg+sysInfo2.comInfoAll[i].deg)/stepSize);
			sysInfo2.comInfoAll[i].deg = delDeg+sysInfo2.comInfoAll[i].deg - degInc*stepSize;
			sysInfo2.comInfoAll[i].initState = \
				min(sysInfo2.comInfoAll[i].nStates - 1, \
					sysInfo2.comInfoAll[i].initState + degInc);
		sysInfo2.comInfoAll[i].update_currentToFail();
	sysInfo2.N0 = [];
	sysInfo2.N1 = [];
	sysInfo2.Nu = [];
	sysInfo2.time = 0;
	sysInfo2.objValue = [];
	sysInfo2.iterInfo = [];

def get_lastStageSol(sysInfo2):
#x,y,z:
	x =	[];
	y = [];
	z = 0;	
	for i in range(sysInfo2.nComponents):
		if sysInfo2.comInfoAll[i].initState >= \
			sysInfo2.comInfoAll[i].nStates - 1:
			#fail
			x.append(1);
			y.append(1);
			z = 1;
			sysInfo2.N1.append(i);
		else:
			x.append(0);
			y.append(0);
			sysInfo2.N0.append(i);
	res = [];
	res.append(x);
	res.append(y);
	res.append(z);
	return res;

def get_alg1Sol(sysInfo2):
	x = [];
	y = [];
	z = 0;

	for i in range(sysInfo2.nComponents):
		#print (sysInfo2.comInfoAll[i].initState)
		if i in sysInfo2.N0:
			x.append(0);
			y.append(0);			
		else:
			x.append(1);
			z = 1;
			if sysInfo2.comInfoAll[i].initState >= \
				sysInfo2.comInfoAll[i].nStates - 1:
				#fail
				y.append(1);
			else:
				y.append(0);
	res = [];
	res.append(x);
	res.append(y);
	res.append(z);
	return res;
	
		
##########################
#start from here
##########################
# some fixed parameters
nComponents = 17;		#fix component numbers	
nStates = 11;			#number of states for components, 0 - (m-1);
cInsp = 1;
cS = 200;				#fix setup cost
inspInterval = 30;		#don't change this
nStages = 5;			
#init component parameter
#gamma distribution is assumed.
# new parameters
#alpha = ct^b
#gam_a = 0.542;
#gam_b = 1/1.147;
gam_a = [1.0824,1];	#alpha
#gam_a = []
#now it a rate parameter, following gamma(kappa, lambda), lambda is a rate as well
gam_b = [8.556,7.654];	#shape & rate
S = [2]*nComponents;	#failure threshold, 10mm
nRep = 102;
c0All = [];
c1All = [];
t0All = [];
t1All = [];
errAll = [];
for rep in range( nRep-1, nRep):			#control replicates, random seed
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
		#for initial state
		ranSeed = (i+1)*baseSeedState + rep+nComponents;			#random seeds
		random.seed(ranSeed);
		tmp = random.randint(0, nStates-1); #0<= state <= nStates - 1
		initState.append(tmp);
		#initState.append(rep);
		#for CM
		ranSeed = (i+1)*baseSeedCM + rep+nComponents;				#random seeds!!!!!!!!!!!!!remember to change it
		random.seed(ranSeed);
		tmp = np.random.normal(60,25,1); 
		#cCM.append(tmp[0]);
		cCM.append(20);
		#for PM
		ranSeed = (i+1)*baseSeedPM + rep+nComponents;				#random seeds!!!!!!!!!!!!!remember to change it
		random.seed(ranSeed);
		tmp = np.random.normal(20, 25,1); 
		#cPM.append(tmp[0]);
		cPM.append(5);
		#for randomize initial age
		ranSeed = (i+1)*baseSeedT + rep+nComponents;				#random seeds!!!!!!!!!!!!!remember to change it
		random.seed(ranSeed);
		tmp = random.uniform(1, 14); 
		initAge.append(tmp);
	#initState = [5,4,7];	
	for i in range(nComponents):
		comInfo = class_info.component_info(i, gam_a, gam_b, nStates,\
								S[i], initState[i], initAge[i], inspInterval, cCM[i], cPM[i], cS);
		sysInfo.add_com(comInfo);

	#get the multi-stage solution.
	#[t]:[0]:x[\omega]
	#[t]:[1]:y[\omega]
	#[t]:[3]:z
	
	#sysInfo1 = copy.deepcopy(sysInfo);
	#main_dynamic_solver.main(sysInfo1);
	
	nRep1 = 1;
	cost1 = 0;
	time1 = 0;
	expCost1 = 0;
	expTime1 = 0;
	for rep1 in range(nRep1):
		timeS = time.clock();
		sysInfo2 = copy.deepcopy(sysInfo);	
		res = [];
		for t in range(nStages):
			print ("============");
			print ("t="+str(t));
			print ("states:");
			for i in range(sysInfo2.nComponents):
				print (sysInfo2.comInfoAll[i].initState);
			if t == nStages - 1:
				tmp = get_lastStageSol(sysInfo2);
			else:		
				main_myAlg.main(sysInfo2, nComponents);
				tmp = get_alg1Sol(sysInfo2);
			print ("results:");
			for i in range(sysInfo2.nComponents):
				print (int(i in sysInfo2.N1));
			res.append(tmp);
			rollHorizon(sysInfo2,rep+rep1,t);	
		'''
		timeE = time.clock();		
		cost1 += get_objValue(sysInfo2, res);
		time1 += timeE - timeS;
	expCost1 = cost1/nRep1;
	expTime1 = time1/nRep1;	
	c0All.append(sysInfo1.objValue);
	c1All.append(expCost1);
	t0All.append(sysInfo1.time);
	t1All.append(expTime1);
	errAll.append(abs(c1All[-1]-c0All[-1])/c0All[-1]);

print("+++++++++++++++");
print("stat all");
print("+++++++++++++++");
print ("rep=("+str(nRep)+","+str(nRep1)+")");
print ("nComponents="+str(nComponents));
print ("-----------------------------");
print ("data for multi-stage dynamic");
print ("-----------------------------");
print ("time all rep:");
print (t0All);
print ("cost all rep:");
print (c0All);
#print ("cost avg:");
#print (c0AllExp);
print ("-----------------------------");
print ("data for rolling two-stage");
print ("-----------------------------");
print ("time all rep:");
print (t1All);
#print ("time avg:");
#print (t1AllExp);
print ("cost all rep:");
print (c1All);
#print ("cost avg:");
#print (c1AllExp);
print ("error ");
print (errAll);
'''
