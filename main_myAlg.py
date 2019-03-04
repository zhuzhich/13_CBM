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

#Notes:
#1. In corollary, N1 cannot start from \emptyset. Othwise, the rN is
#   inf when in failure state, and it will make all component in N 
#	go to N1.
#2. We need to prove if rN < p(N0, N1), there exists a subset of N, N'
#	such that rN' < p(N0, N1). 
#	Othwise, the Algorithm 1 needs to be rewrite.
#
def get_rTildeN(sysInfo, N):
	out = [];
	if len(N) == 1:
		if sysInfo.comInfoAll[N[0]].initState == 0:
			out.append(float("Inf"));
			out.append(float("Inf"));
			return out;
		elif sysInfo.comInfoAll[N[0]].initState >= \
			sysInfo.comInfoAll[N[0]].nStates - 1:
			out.append(float("-Inf"));
			out.append(float("-Inf"));
			return out;
	rho = 0;
	rN = 1;
	for i in N:
		rho_i = (sysInfo.comInfoAll[i].cPM - \
				(sysInfo.comInfoAll[i].currentToFail - 
				sysInfo.comInfoAll[i].newToFail)*\
				sysInfo.comInfoAll[i].cCM)/sysInfo.cS;
		rho += rho_i;
		rN = rN*(1 - sysInfo.comInfoAll[i].newToFail)/ \
			(1 - sysInfo.comInfoAll[i].currentToFail);
	out.append(rho/(rN - 1));#not empty
	out.append((rho+1)/(rN - 1)); #empty
	return out;

def get_sTildeN(sysInfo, N):
	out = [];
	if len(N) == 1:
		if sysInfo.comInfoAll[N[0]].initState == 0:
			out.append(float("Inf"));
			out.append(float("Inf"));
			return out;
		elif sysInfo.comInfoAll[N[0]].initState >= \
			sysInfo.comInfoAll[N[0]].nStates:
			out.append(float("-Inf"));
			out.append(float("-Inf"));
			return out;
	rho = 0;
	sN = 1;
	for i in N:
		rho_i = (sysInfo.comInfoAll[i].cPM - \
				(sysInfo.comInfoAll[i].currentToFail - 
				sysInfo.comInfoAll[i].newToFail)*\
				sysInfo.comInfoAll[i].cCM)/sysInfo.cS;
		rho += rho_i;
		sN = sN*(1 - sysInfo.comInfoAll[i].currentToFail)/ \
			(1 - sysInfo.comInfoAll[i].newToFail);
	out.append(rho/(1 - sN));#not empty
	out.append((rho+1)/(1 - sN)); #empty
	return out;	
		



def get_objValue(sysInfo, N0, N1):
	objValue = 0;

	N2 = [];
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
	if len(N1) > 0:
		stage1 += sysInfo.cS;
	#second stage:
	stage2 = 0;
	probCS = 1;
	for i in range(sysInfo.nComponents):
		if i in N0:
			failProb = sysInfo.comInfoAll[i].transProb(\
						sysInfo.comInfoAll[i].initState,\
						sysInfo.comInfoAll[i].nStates - 1,\
						sysInfo.inspItvl);
		elif i in N1:
			failProb = sysInfo.comInfoAll[i].transProb(\
						0,\
						sysInfo.comInfoAll[i].nStates - 1,\
						sysInfo.inspItvl);
		else:
			continue;
		stage2 += failProb*sysInfo.comInfoAll[i].cCM;
		probCS = probCS*(1-failProb);
	#for objValue1
	objValue = stage1 + stage2 + (1-probCS)*sysInfo.cS;

	return objValue;

def probN0N1(sysInfo, N0, N1):
	out = 1;
	for i in N0:
		out = out*(1-sysInfo.comInfoAll[i].transProb(\
				sysInfo.comInfoAll[i].initState, \
				sysInfo.comInfoAll[i].nStates - 1, \
				sysInfo.comInfoAll[i].inspItvl));
	for i in N1:
		out = out*(1-sysInfo.comInfoAll[i].transProb(\
			0, \
			sysInfo.comInfoAll[i].nStates - 1, \
			sysInfo.comInfoAll[i].inspItvl));	
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
	while (len(Nu) != 0 and j <= J):
		iterInfo.append(j);
		setSj = [];
		for i in itertools.combinations( Nu, j):
			setSj.append(list(i));
		N_prime = [];
		while ( Nu != N_prime):
			N_prime  = copy.deepcopy( Nu);
			for N in setSj:
				rTildeN = get_rTildeN(sysInfo,N);
				sTildeN = get_sTildeN(sysInfo,N);
				probR = probN0N1(sysInfo, list(set(N_all)-set(N1_opt)),N1_opt);
				probS = probN0N1(sysInfo, N0_opt,list(set(N_all)-set(N0_opt)))
				#corollary 1
				if len(N1_opt) != 0 and \
					rTildeN[0] < probR and\
					set(N) <= set( Nu):
						N1_opt = N + N1_opt;
						setSj.remove(N);
						Nu = list(set( Nu) - set(N));
				elif len(N1_opt) == 0 and \
					rTildeN[1] < probR and\
					set(N) <= set( Nu):
						N1_opt = N + N1_opt;
						setSj.remove(N);
						Nu = list(set( Nu) - set(N));
				#corollary 2
				elif len( Nu)!= j and \
					sTildeN[0] > probS and\
					set(N) <= set( Nu):
						N0_opt = N + N0_opt;
						setSj.remove(N);
						Nu = list(set( Nu) - set(N));
				elif len( Nu) == j and \
					sTildeN[1] > probS and\
					set(N) <= set( Nu):
						N0_opt = N + N0_opt;
						setSj.remove(N);
						Nu = list(set( Nu) - set(N));
		J = len(Nu);
		if j+1 <= J:
			j += 1;
		else:
			j = 1;
	sysInfo.N0 = N0_opt;
	sysInfo.N1 = N1_opt;
	sysInfo.Nu =  Nu;
	sysInfo.iterInfo = iterInfo;
	sysInfo.objValue = get_objValue(sysInfo, sysInfo.N0, sysInfo.N1);

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
	end_time = time.clock();
	sysInfo.time = end_time-start_time;





















