#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2018: Zhicheng Zhu. All right reserved.

#Info:
#main file to solve multi-stage DEF of CBM model by using linearization and solver
#first-order approximation.
#
#Last update: 10/18/2018
#!/usr/bin/python

from __future__ import print_function
import sys
import cplex
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
		
def get_child_nodes(node, sysInfo):
	#find/generate child nodes 
	#(t, 0), (t, 1), ..., (t, m**tn)
	m = sysInfo.comInfoAll[0].nStates;
	n = sysInfo.nComponents;
	numOutcome = m**n;
	start = node*numOutcome;
	childNodes = list(range(start, start + numOutcome));
	#we only return next stage nodes index.
	return childNodes;

def node_2_outcome(node, sysInfo):
	#translate a node to outcome:
	m = sysInfo.comInfoAll[0].nStates;
	n = sysInfo.nComponents;
	numOutcome = m**n;
	outcome = node % numOutcome;
	return outcome;
#######################################
#1. initialization, START FROM HERE!!!.
#######################################
#init system parameter
start_time = time.clock();

#init system parameter
nComponents = 2;
nStates = 4;					#number of states for components, 0 - (m-1);
nStages = 6;
initState = [3,2];

inspInterval = 10;
cS = 20;				#setup cost
cInsp = 1;
sysInfo = system_info(nComponents, nStages, inspInterval, cS, cInsp);
nOrder = 0;				#order of approximation
#init component parameter
#gamma distribution is assumed.
gam_a =   	[1]*nComponents;
gam_b =   	[5]*nComponents;
S = 		[60]*nComponents;	#failure threshold
cCM = [20]*nComponents;
cPM = [5]*nComponents;

for i in range(nComponents):
	comInfo = component_info(i, gam_a[i], gam_b[i], nStates,\
							S[i], initState[i], cCM[i], cPM[i]);
	sysInfo.add_com(comInfo);


########################################
#2. build multi-stage DEF model and run 
########################################
start_time = time.clock();

#2.1 
#    (1) get scenario combinations
omega = [];
for i in itertools.product(list(range(nStates)), repeat = sysInfo.nComponents):
	omega.append(list(i));
	
'''	
#no set j when order = 0;
#	(2) get subsets which cardinality = j
setS = [];	#start from j = 2
for j in range(2, sysInfo.nComponents + 1):	#sysInfo.nComponents >=2;
	setSj = [];
	for i in itertools.combinations(list(range(sysInfo.nComponents)), j):
		setSj.append(list(i));
	setS.append(setSj);
'''		

#    get coeA and coeB 
coeA = [];				#scen*n*scen
coeB = [];				#scen*n*scen
coeAInit = []; 			#n*scen, store init coeA
coeBInit = [];			#n*scen, store init coeB
for w1 in range(len(omega)):
	stateFrom = omega[w1];
	coeAW1 = [];
	coeBW1 = [];
	for i in range(sysInfo.nComponents):
		biw = [];
		aiw = [];		
		for w2 in range(len(omega)):
			comStatesTo = omega[w2];
			comIFrom = stateFrom[i];
			comITo = comStatesTo[i];
			tmp = sysInfo.comInfoAll[i].transProb(comIFrom, comITo, sysInfo.inspItvl);
			biw.append(tmp);
			aiw.append(sysInfo.comInfoAll[i].transProb(0, comITo, sysInfo.inspItvl) - tmp);
		coeAW1.append(aiw);
		coeBW1.append(biw);
	if stateFrom == initState:
		coeAInit = coeAW1;
		coeBInit = coeBW1;
	coeA.append(coeAW1);
	coeB.append(coeBW1);
	
cpx = cplex.Cplex();		#init solver
cpx.objective.set_sense(cpx.objective.sense.minimize);
	
#2.2 add decision variables

#add X
varX = [];
dictX = {};
for stageIdx in range(nStages):
	nodeNum = sysInfo.comInfoAll[0].nStates ** (stageIdx * sysInfo.nComponents);
	for node in range(nodeNum):		#nodes in current stage
		for i in range(sysInfo.nComponents):
			scripts = str(i) + str(stageIdx) + str(node);
			nameTmp = "x"+scripts;
			dictX[scripts] = nameTmp;
			varX.append(cpx.variables.get_num());
			objCoe = 0;		
			if stageIdx == 0:
				objCoe = sysInfo.comInfoAll[i].cPM;			
			cpx.variables.add(obj = [objCoe], lb = [0.0], ub=[1.0], types=["B"], names=[nameTmp]);
		
#add Y
varY = [];
dictY = {};
for stageIdx in range(nStages):
	nodeNum = sysInfo.comInfoAll[0].nStates ** (stageIdx * sysInfo.nComponents);
	for node in range(nodeNum):		#nodes in current stage
		for i in range(sysInfo.nComponents):
			scripts = str(i)+str(stageIdx)+str(node);
			nameTmp = "y" + scripts;
			dictY[scripts] = nameTmp;
			varY.append(cpx.variables.get_num());
			objCoe = 0;		
			if stageIdx == 0:
				objCoe = sysInfo.comInfoAll[i].cCM - sysInfo.comInfoAll[i].cPM;					
			cpx.variables.add(obj = [objCoe], lb = [0.0], ub=[1.0], types=["B"], names=[nameTmp]);
							
#add Z
varZ = [];
dictZ = {};
for stageIdx in range(nStages):
	nodeNum = sysInfo.comInfoAll[0].nStates ** (stageIdx * sysInfo.nComponents);
	for node in range(nodeNum):		#nodes in current stage
		scripts = str(stageIdx) + str(node);
		nameTmp = "z" + scripts;
		dictZ[scripts] = nameTmp;
		varZ.append(cpx.variables.get_num());
		objCoe = 0;
		if stageIdx == 0: 
			objCoe = sysInfo.cS;
		cpx.variables.add(obj = [objCoe], lb = [0.0], ub=[1.0], types=["B"], names=[nameTmp]);	


#add Theta
varTheta = [];
dictTheta = {};
for stageIdx in range(1, nStages):
	nodeNum = sysInfo.comInfoAll[0].nStates ** (stageIdx * sysInfo.nComponents);
	for node in range(nodeNum):		#nodes in current stage
		coeTmp = 0;
		if stageIdx == 1:
			coeTmp = 1;
			for i in range(sysInfo.nComponents):
				coeTmp = coeTmp * coeBInit[i][node];
		#print ("ThetacoeTmp=" + str(coeTmp));
		scripts = str(stageIdx) + str(node);
		nameTmp = "th" + scripts;
		dictTheta[scripts] = nameTmp;
		varTheta.append(cpx.variables.get_num());
		cpx.variables.add(obj = [coeTmp], lb = [0.0], ub=[cplex.infinity], types=["C"], names=[nameTmp]);			
#add V
varV = [];
dictV= {};
for stageIdx in range(nStages - 1):
	nodeNum = sysInfo.comInfoAll[0].nStates ** (stageIdx * sysInfo.nComponents);
	for curNode in range(nodeNum):
		childNodes = get_child_nodes(curNode, sysInfo);
		for chNode in childNodes:
			for i in range(sysInfo.nComponents):
				#v corresponds to cardinality set when cardinality j = 1.
				if stageIdx != 0:
					coeTmp = 0;
				else:
					coeTmp = coeAInit[i][chNode];
					for r in range(sysInfo.nComponents):
						if r != i:
							coeTmp = coeTmp * coeBInit[r][chNode];
				#print ("VcoeTmp=" + str(coeTmp));
				scripts = str(i) + str(stageIdx) + str(curNode) + str(chNode);
				nameTmp = "v" + scripts;
				dictV[scripts] = nameTmp;
				varV.append(cpx.variables.get_num());
				#continuous variable
				cpx.variables.add(obj = [coeTmp], lb = [0.0], ub=[cplex.infinity], types=["C"], names=[nameTmp]);
'''					
#no  W & U when order = 1
#add W
varW = [];
dictW = {};
for stageIdx in range(nStages - 1):
	nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);
	for curNode in range(nodeNum):
		childNodes = get_child_nodes(curNode, sysInfo);
		for chNode in childNodes:
			for j in range(2, sysInfo.nComponents+1):
				#cardinality starts from 2 to n.
				setSj = setS[j-2];
				for k in range(len(setSj)):
					if stageIdx != 0:
						coeTmp = 0;
					else:
						setSjk = setSj[k];
						coeTmp = 1;
						for i in range(sysInfo.nComponents):
							if i in setSjk:
								coeTmp = coeTmp*coeAInit[i][chNode];
							else:
								coeTmp = coeTmp*coeBInit[i][chNode];		
					#print ("WcoeTmp=" + str(coeTmp));
					scripts = str(j) + str(k) + str(stageIdx) + str(curNode) + str(chNode);
					nameTmp = "w" + scripts;
					dictW[scripts] = nameTmp;
					varW.append(cpx.variables.get_num());
					#continuous variable
					cpx.variables.add(obj = [coeTmp], lb = [0.0],  ub=[cplex.infinity], types=["C"], names=[nameTmp]);

#add U: auxilary variable that used in w
varU = [];
dictU = {};
for stageIdx in range(nStages - 1):
	nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);
	for node in range(nodeNum):		#nodes in current stage
		for j in range(2, sysInfo.nComponents+1):
			#cardinality starts from 2 to n.
			setSj = setS[j-2];
			for k in range(len(setSj)):
				scripts = str(j) + str(k) + str(stageIdx) + str(node);
				nameTmp = "u" + scripts;
				dictU[scripts] = nameTmp;
				varU.append(cpx.variables.get_num());
				cpx.variables.add(obj = [0], lb = [0.0], ub=[1.0], types=["B"], names=[nameTmp]);				
				
'''				
## 2.2 add constraints
# 1
for stageIdx in range(nStages):
	nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);	
	for node in range(nodeNum):
		coefNameZ = dictZ[str(stageIdx) + str(node)];
		for i in range(sysInfo.nComponents):
			coefNameX = dictX[str(i) + str(stageIdx) + str(node)];
			cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([coefNameX, coefNameZ], [1, -1])], senses=["L"], range_values=[0.0], rhs=[0]);		

# 2 & 3
for stageIdx in range(nStages):
	nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);
	curOutcome = 0;			#distinct outcome index.
	for node in range(nodeNum):
		coefValueVec = [];
		coefNameVec = [];
		if stageIdx == 0:
			curStates = initState;
		else:
			curStates = omega[curOutcome];
		curOutcome += 1;
		if curOutcome == len(omega):
			curOutcome = 0;
		for i in range(sysInfo.nComponents):
			# 2
			curStatesI = curStates[i];
			coefNameY = dictY[str(i) + str(stageIdx) + str(node)];
			coefValueY = curStatesI;
			cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([coefNameY],[-coefValueY])], senses=["L"], range_values=[0.0], rhs=[sysInfo.comInfoAll[i].nStates-2-curStatesI]);
			# 3
			nameIdxScriptX = str(i) + str(stageIdx) + str(node);
			coefNameX = dictX[nameIdxScriptX];
			coefValueX = -1;
			coefValueY = 1;		#value changed here for 3rd constraint
			cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([coefNameY, coefNameX],[coefValueY, coefValueX])], senses=["L"], range_values=[0.0], rhs=[0.0]);

# 4: tooooo complex:
# in 4, theta starts from stage 1 to nStages - 2.
for stageIdx in range(1, nStages - 1):
	nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);
	for node in range(nodeNum):
		# do the first part
		coefNameVec = [];
		coefValueVec = [];
		nameTmp = dictTheta[str(stageIdx) + str(node)];
		coefNameVec.append(nameTmp);
		coefValueVec.append(-1);
		for i in range(sysInfo.nComponents):
			#add x
			nameTmp = dictX[str(i) +str(stageIdx) + str(node)];
			coefNameVec.append(nameTmp);
			coefValueVec.append(sysInfo.comInfoAll[i].cPM);
			#add y
			nameTmp = dictY[str(i) +str(stageIdx) + str(node)];
			coefNameVec.append(nameTmp);
			coefValueVec.append(sysInfo.comInfoAll[i].cCM - sysInfo.comInfoAll[i].cPM);
		#add z
		nameTmp = dictZ[str(stageIdx) + str(node)];
		coefNameVec.append(nameTmp);
		coefValueVec.append(sysInfo.cS);
		#do the second part
		childNodes = get_child_nodes(node, sysInfo);
		for chNode in childNodes:
			#within the second part...
			#part 1
			nameTmp = dictTheta[str(stageIdx+1) + str(chNode)];
			stateFromIdx = node_2_outcome(node, sysInfo);
			stateFrom = omega[stateFromIdx];
			stateToIdx = node_2_outcome(chNode, sysInfo);
			stateTo = omega[stateToIdx];			
			valueTmp = 1;
			for i in range(sysInfo.nComponents):
				valueTmp = valueTmp * coeB[stateFromIdx][i][stateToIdx];
				if valueTmp == 0:
					break;		#make it faster;
			coefNameVec.append(nameTmp);
			coefValueVec.append(valueTmp);
			#print (valueTmp);
			#part 2
			for i in range(sysInfo.nComponents):
				nameTmp = dictV[str(i) + str(stageIdx) + str(node) + str(chNode)];		
				valueTmp = coeA[stateFromIdx][i][stateToIdx];
				for r in range(sysInfo.nComponents):
					if r != i:
						valueTmp = valueTmp * coeB[stateFromIdx][r][stateToIdx];	
					if valueTmp == 0:
						break;	#make it faster
				coefNameVec.append(nameTmp);
				coefValueVec.append(valueTmp);	
		cpx.linear_constraints.add(lin_expr=[cplex.SparsePair(coefNameVec,coefValueVec)], senses=["E"], range_values=[0.0], rhs=[0.0]);

'''						
			#part 3:
			for j in range(2, sysInfo.nComponents + 1):
				setSj = setS[j - 2];				#setS starts from 2
				for k in range(len(setSj)):
					nameTmp = dictW[str(j) + str(k)  + str(stageIdx) + str(node) + str(chNode)];
					valueTmp = 1;
					setSjk = setSj[k];
					for i in range(sysInfo.nComponents):
						if i in setSjk:
							valueTmp = valueTmp * coeA[stateFromIdx][i][stateToIdx];
						else:
							valueTmp = valueTmp * coeB[stateFromIdx][i][stateToIdx];
						if valueTmp == 0:
							break;	#make it faster
					coefNameVec.append(nameTmp);
					coefValueVec.append(valueTmp);	
		#theta is stage * node
'''
		
					
# 5: theta at last stage
stageIdx = nStages - 1;
nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);
for node in range(nodeNum):
	coefNameVec = [];
	coefValueVec = [];
	nameTmp = dictTheta[str(stageIdx) + str(node)];
	coefNameVec.append(nameTmp);
	coefValueVec.append(-1);
	for i in range(sysInfo.nComponents):
		#add x
		nameTmp = dictX[str(i) +str(stageIdx) + str(node)];
		coefNameVec.append(nameTmp);
		coefValueVec.append(sysInfo.comInfoAll[i].cPM);
		#add y
		nameTmp = dictY[str(i) +str(stageIdx) + str(node)];
		coefNameVec.append(nameTmp);
		coefValueVec.append(sysInfo.comInfoAll[i].cCM - sysInfo.comInfoAll[i].cPM);
	#add z
	nameTmp = dictZ[str(stageIdx) + str(node)];
	coefNameVec.append(nameTmp);
	coefValueVec.append(sysInfo.cS);
	cpx.linear_constraints.add(lin_expr=[cplex.SparsePair(coefNameVec,coefValueVec)], senses=["E"], range_values=[0.0], rhs=[0.0]);

# 6: add linearization of V:
# There are 4 parts in this section:
upperM = 10000;		#upper bound of theta 
for stageIdx in range(0, nStages - 1):
	nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);
	for node in range(nodeNum):
		childNodes = get_child_nodes(node, sysInfo);
		for i in range(sysInfo.nComponents):
			nameTmpX = dictX[str(i) + str(stageIdx) + str(node)];
			valueTmpX = -upperM;
			for chNode in childNodes:
				nameTmpV = dictV[str(i) + str(stageIdx) + str(node) + str(chNode)];
				valueTmpV = 1;
				# part 1
				cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([nameTmpX, nameTmpV],[valueTmpX, valueTmpV])], senses=["L"], range_values=[0.0], rhs=[0.0]);
				# part 2
				nameTmpTheta = dictTheta[str(stageIdx + 1) + str(chNode)];
				valueTmpTheta = -1;
				cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([nameTmpTheta, nameTmpV],[valueTmpTheta, valueTmpV])], senses=["L"], range_values=[0.0], rhs=[0.0]);
				#part 3
				cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([nameTmpV, nameTmpTheta, nameTmpX],[valueTmpV, valueTmpTheta, valueTmpX])], senses=["G"], range_values=[0.0], rhs=[valueTmpX]);
				# part 4 is added when adding variable V
'''
# 7: add linearization of W:
# There are 4 parts of W
for stageIdx in range(0, nStages - 1):
	nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);
	for node in range(nodeNum):
		childNodes = get_child_nodes(node, sysInfo);
		for chNode in childNodes:
			for j in range(2, sysInfo.nComponents + 1):
				setSj = setS[j - 2];
				for k in range(len(setSj)):
					nameTmpW = dictW[str(j) + str(k) + str(stageIdx) + str(node) + str(chNode)];
					valueTmpW = 1;
					nameTmpU = dictU[str(j) + str(k) + str(stageIdx) + str(node)];
					valueTmpU = -upperM;
					# part 1
					cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([nameTmpW, nameTmpU],[valueTmpW, valueTmpU])], senses=["L"], range_values=[0.0], rhs=[0.0]);
					# part 2
					nameTmpTheta = dictTheta[str(stageIdx + 1) + str(chNode)];
					valueTmpTheta = -1;					
					cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([nameTmpW, nameTmpTheta],[valueTmpW, valueTmpTheta])], senses=["L"], range_values=[0.0], rhs=[0.0]);
					# part 3
					cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([nameTmpW, nameTmpTheta, nameTmpU],[valueTmpW, valueTmpTheta, valueTmpU])], senses=["G"], range_values=[0.0], rhs=[valueTmpU]);	
					# part 4 is added when adding variable W

# 8: add linearization of U:
# There are 3 parts of U
for stageIdx in range(nStages - 1):
	nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);
	for node in range(nodeNum):
		for j in range(2, sysInfo.nComponents + 1):
			setSj = setS[j - 2];
			for k in range(len(setSj)):	
				setSjk = setSj[k];
				nameTmpU = dictU[str(j) + str(k) + str(stageIdx) + str(node)];
				valueTmpU = 1;	
				namePart2 = [];
				valuePart2 = [];
				namePart2.append(nameTmpU);
				valuePart2.append(valueTmpU);
				for i in setSjk:
					nameTmpX = dictX[str(i) + str(stageIdx) + str(node)];
					valueTmpX = -1;		
					#part 1:
					cpx.linear_constraints.add(lin_expr=[cplex.SparsePair([nameTmpU, nameTmpX],[valueTmpU, valueTmpX])], senses=["L"], range_values=[0.0], rhs=[0.0]);	
					#prepare for part 2:
					namePart2.append(nameTmpX);
					valuePart2.append(valueTmpX);
				#part 2
				cpx.linear_constraints.add(lin_expr=[cplex.SparsePair(namePart2, valuePart2)], senses=["G"], range_values=[0.0], rhs=[-j + 1]);	# -(j - 1)		
				# part 3 is added when adding variable U
'''	
########################################
#3. solve and result handling 
########################################					
end_time = time.clock();
time_elapsed0 = end_time - start_time;


start_time = time.clock();
				
cpx.solve();
solution = cpx.solution;
#obj value
objValues = solution.get_objective_value();

#get solutions
solutionAll = solution.get_values();

#get X
minTmp = varX[0];
maxTmp = varX[-1] + 1;
solutionX = solutionAll[minTmp:maxTmp];

#get Y
minTmp = varY[0];
maxTmp = varY[-1] + 1;
solutionY = solutionAll[minTmp:maxTmp];

#get Z
minTmp = varZ[0];
maxTmp = varZ[-1] + 1;
solutionZ = solutionAll[minTmp:maxTmp];

#get theta
minTmp = varTheta[0];
maxTmp = varTheta[-1] + 1;
solutionTheta = solutionAll[minTmp:maxTmp];
'''
#get V
minTmp = varV[0];
maxTmp = varV[-1] + 1;
solutionV = solutionAll[minTmp:maxTmp];

#get W
minTmp = varW[0];
maxTmp = varW[-1] + 1;
solutionW = solutionAll[minTmp:maxTmp];

#get U 
minTmp = varU[0];
maxTmp = varU[-1] + 1;
solutionU = solutionAll[minTmp:maxTmp];
'''
end_time = time.clock();

time_elapsed = end_time - start_time;

f = open("log3.txt", "w");
old = sys.stdout;
sys.stdout = f;

print ("\n===============================main_multi_DEF_solver_0, (m, n, t)=(%d,%d,%d)============" 
		%(nStates, sysInfo.nComponents, nStages));
print ("loading time is %f"  %time_elapsed0);
print ("calculation time is %f"  %time_elapsed);
print ("objValues:");
print (objValues);

countX = 0;
countY = 0;
countZ = 0;
countV = 0;
countW = 0;
countU = 0;
countTheta = 0;
for stageIdx in range(nStages):
	nodeNum = sysInfo.comInfoAll[0].nStates**(stageIdx*sysInfo.nComponents);
	for node in range(nodeNum):
		print ("=======(stage, scen) = (%d, %d)========" %(stageIdx,node));
		#get X Y Z theta
		solX = [];
		solY = [];
		solZ = solutionZ[countZ];
		countZ += 1;
		solTheta = [];
		if stageIdx != 0:
			solTheta = solutionTheta[countTheta];
			countTheta += 1;
		for i in range(sysInfo.nComponents):
			solX.append(solutionX[countX]);
			countX += 1;
			solY.append(solutionY[countY]);
			countY += 1;
		print ("solutionX:");
		print (solX);
		print ("solutionY:");
		print (solY);
		print ("solutionZ:");
		print (solZ);
		print ("solutionTheta:");
		print (solTheta);
		'''
		#get U
		if stageIdx == nStages - 1:				#last stage, no U V W
			continue;				
		solU = [];
		for j in range(2, sysInfo.nComponents + 1):
			setSj = setS[j - 2];
			for k in range(len(setSj)):
				solU.append(solutionU[countU]);
				countU += 1;
		print ("solutionU:");
		print (solU);
		#get v and w
		childNodes = get_child_nodes(node, sysInfo);
		solV = [];
		solW = [];					
		for chNode in childNodes:
			#get V
			solVTmp = [];			
			for i in range(sysInfo.nComponents):
				solVTmp.append(solutionV[countV]);
				countV += 1;
			solV.append(solVTmp);
			#get W
			solWTmp = [];
			for j in range(2, sysInfo.nComponents + 1):
				setSj = setS[j - 2];
				for k in range(len(setSj)):
					solWTmp.append(solutionW[countW]);
					countW += 1;
			solW.append(solWTmp);
		print ("solutionV:");
		print (solV);	
		print ("solutionW:");
		print (solW);		
		print ("===================\n");
		'''
		
'''
print ("=======coeA======");
print (coeA);
print ("=======coeB======");
print (coeB);
print ("=======coeU======");
print (coeU);
print ("=======coeX======");
print (coeX);
print ("=======costTerm======");
print (consTerm);
'''
## 4. end of file 
sys.stdout = old;
f.close();		
		
		
		