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
import itertools
import time
import class_info
import copy
import random

def get_objValue(sysInfo, N0, N1):
	objValue = 0;
	#first stage
	stage1 = 0;
	for i in range(sysInfo.nComponents):
		#if i not in N0 and i not in N1:
		#	N2.append(i);
		if i in N1:
			if sysInfo.comInfoAll[i].initState >= sysInfo.comInfoAll[i].nStates - 1:
				stage1 += sysInfo.comInfoAll[i].cCM;
			else:
				stage1 += sysInfo.comInfoAll[i].cPM;
		if i in N0:
			if sysInfo.comInfoAll[i].initState >= sysInfo.comInfoAll[i].nStates - 1:
				stage1 += sysInfo.comInfoAll[i].cCM;
				print ("find bug!!!!N0 failed!!!")
	if stage1 > 0:
		stage1 += sysInfo.cS;
	#second stage:
	stage2 = 0;
	probCS = 1;
	for i in range(sysInfo.nComponents):
		if i in N0:
			failProb = sysInfo.comInfoAll[i].currentToFail;
		elif i in N1:
			failProb = sysInfo.comInfoAll[i].newToFail
		else:
			continue;
		stage2 += failProb*sysInfo.comInfoAll[i].cCM;
		probCS = probCS*(1-failProb);
	#for objValue1
	objValue = stage1 + stage2 + (1-probCS)*sysInfo.cS;

	return objValue;

def get_delta(sysInfo,N0_opt,N1_opt,Nu,N):
	deltaRS = []; #[0] is deltaR, [1] is deltaS.
					#[0], not empty, [1] empty.
	#1. get rho
	rho = 0;
	for i in N:
		rho_i = (sysInfo.comInfoAll[i].cPM - \
				(sysInfo.comInfoAll[i].currentToFail - 
				sysInfo.comInfoAll[i].newToFail)*\
				sysInfo.comInfoAll[i].cCM)/sysInfo.cS;
		rho += rho_i;
			
	# 2 cal deltaR
	# 2.1 	get prob(N0_opt+N_u,N1_opt)
	prob = 1;
	for i in N0_opt:
		if sysInfo.comInfoAll[i].currentToFail == 1:
			print ("ERROR MY ALG 1!!!!!!",__LINE__);
		prob = prob*(1-sysInfo.comInfoAll[i].currentToFail);
	for i in Nu:
		if sysInfo.comInfoAll[i].currentToFail == 1:
			print ("ERROR MY ALG 1!!!!!!",__LINE__);
		prob = prob*(1-sysInfo.comInfoAll[i].currentToFail);	
	for i in N1_opt:
		prob = prob*(1-sysInfo.comInfoAll[i].newToFail);
		
	# 2.2 get rN
	rN = 1;
	for i in N:
		if sysInfo.comInfoAll[i].currentToFail == 1:
			print ("ERROR MY ALG 1!!!!!!",__LINE__);
		rN = rN*(1-sysInfo.comInfoAll[i].newToFail)\
		/(1-sysInfo.comInfoAll[i].currentToFail);	
	rN = rN - 1;
	

	#2.3 	get deltaR
	deltaR = [];
	if rN*prob < (0.1)**9:
		deltaR.append(float("Inf"));#not empty, never move to N1
		deltaR.append(float("Inf"));#not empty
	else:
		deltaR.append(rho/(rN*prob));#not empty
		deltaR.append((1+rho)/(rN*prob));#empty
	
	# 3 cal deltaS
	# 3.1 	get prob(N0_opt,N1_opt+N_u)
	prob = 1;
	for i in N0_opt:
		if sysInfo.comInfoAll[i].currentToFail == 1:
			print ("ERROR MY ALG 1!!!!!!",__LINE__);
		prob = prob*(1-sysInfo.comInfoAll[i].currentToFail);
	for i in Nu:
		if sysInfo.comInfoAll[i].currentToFail == 1:
			print ("ERROR MY ALG 1!!!!!!",__LINE__);
		prob = prob*(1-sysInfo.comInfoAll[i].newToFail);		
	for i in N1_opt:
		prob = prob*(1-sysInfo.comInfoAll[i].newToFail);	
	# 3.2 get sN
	sN = 1;
	for i in N:
		sN = sN*(1-sysInfo.comInfoAll[i].currentToFail)\
		/(1-sysInfo.comInfoAll[i].newToFail);	
	sN = 1 - sN;
	
	#3.3 	get deltaR
	deltaS = [];
	if sN*prob < (0.1)**9:
		deltaS.append(float("Inf"));	#always move to N0
		deltaS.append(float("Inf"));
	else:
		deltaS.append(rho/(sN*prob));#not empty
		deltaS.append((1+rho)/(sN*prob));#empty

	deltaRS.append(deltaR);
	deltaRS.append(deltaS);
	return deltaRS;

def heuristic_search(sysInfo,N0_opt,N1_opt,Nu):
	#try Nu=0 and Nu = 1:
	N0_p = Nu;
	N1_p = [];
	C_opt = get_objValue(sysInfo, N0_opt+N0_p, N1_opt+N1_p);
	N0_tmp = [];
	N1_tmp = Nu;
	C = get_objValue(sysInfo, N0_opt+N0_tmp, N1_opt+N1_tmp); 
	if C < C_opt:
		N0_p = N0_tmp;
		N1_p = N1_tmp;
		C_opt = C;
	m = 1;
	M = 100;
	while (m <= M):
		m += 1;
		N0_tmp = [];
		N1_tmp = [];
		for i in Nu:
			if random.uniform(0,1)<0.5:
				N0_tmp.append(i);
			else:
				N1_tmp.append(i);
		C = get_objValue(sysInfo, N0_opt+N0_tmp, N1_opt+N1_tmp); 
		if C < C_opt:
			N0_p = N0_tmp;
			N1_p = N1_tmp;
			C_opt = C;
	out = [];
	out.append(N0_p);
	out.append(N1_p);
	return out;
	
#
#
#J = J^th order
def main(sysInfo, J):
	start_time = time.clock();
	#get Sj
	N0_opt = [];
	N1_opt = [];
	N_all = list(range(sysInfo.nComponents));
	Nu = list(range(sysInfo.nComponents));# Nu[j]: undetermine set after j=1,2,...,n orders
	iterInfo = [];
	j = 1;
	for i in N_all:
		if sysInfo.comInfoAll[i].initState >=  sysInfo.comInfoAll[i].nStates - 1:
			N1_opt.append(i);
			N = [i];
			Nu = list(set(Nu) - set(N));
	while (len(Nu) != 0 and j <= J ):
		#print (j);
		#print (Nu)
		iterInfo.append(j);
		setSj = [];
		for i in itertools.combinations(Nu, j):
			setSj.append(list(i));
		j_prime = len(Nu);
		N0 = [];
		N1 = [];
	
		for N in setSj:
			deltaRS = get_delta(sysInfo,N0_opt,N1_opt,Nu,N);
			deltaR = deltaRS[0];
			deltaS = deltaRS[1];
			
			#prop 2	
			if len(N1_opt) != 0 and \
				deltaR[0] < 1:
					N1 = N + N1;
			elif len(N1_opt) == 0 and \
				deltaR[1] < 1:
					N1 = N + N1
			#prop 3
			elif len(N1_opt+Nu)>len(N) and \
				deltaS[0] >= 1:
				N0 = N + N0;
			elif len(N1_opt+Nu)==len(N) and \
				deltaS[1] >= 1:
				N0 = N + N0;	
		Nu = list(set(Nu)-set(N0)-set(N1));
		N0_opt = N0_opt + N0;
		N1_opt = N1_opt + N1;
		if j_prime > len(Nu):
			j = 1;
		else:
			j += 1;
	
	N0_p = [];
	N1_p = [];
	#print ("outttt")
	#print (Nu);
	if len(Nu) > 0:
		ptn = heuristic_search(sysInfo,N0_opt,N1_opt,Nu);
		N0_p = ptn[0];
		N1_p = ptn[1];
		

	sysInfo.N0 = N0_opt+N0_p;
	sysInfo.N1 = N1_opt+N1_p;
	sysInfo.Nu =  Nu;
	sysInfo.Nu_0 = N0_p;
	sysInfo.Nu_1 = N1_p;	
	sysInfo.iterInfo = iterInfo;
	sysInfo.objValue = get_objValue(sysInfo, sysInfo.N0, sysInfo.N1);
	end_time = time.clock();
	sysInfo.time = end_time-start_time;
'''
	if  len(Nu) == 0:
		sysInfo.objValue = get_objValue(sysInfo, sysInfo.N0, sysInfo.N1);
	else:#heuristic
		#only check 0 or 1 case:
		#all 0
		N0 = copy.deepcopy( Nu);
		N1 = [];
		c_opt = get_objValue(sysInfo, N0_opt+N0, sysInfo.N1);
		#all 1
		c_tmp = get_objValue(sysInfo, N0_opt, N1_opt+ Nu);
		if c_tmp < c_opt:
			c_opt = c_tmp;
			N0 = [];
			N1 = copy.deepcopy( Nu);
		sysInfo.Nu_0 = N0;
		sysInfo.Nu_1 = N1;
		sysInfo.objValue = c_opt;
'''		






















