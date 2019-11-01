from scipy.stats import gamma
from scipy.stats import f
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
		if (1):
		## random effect model.
			if stateTo >= self.nStates - 1:
				deltaDeg = self.failTsh - degFrom;
				tmp = self.gammaAlpha[0]*inspItvl**self.gammaAlpha[1];
				x = deltaDeg*self.gammaBeta[0]*self.gammaBeta[1]/(tmp);
				prob = 1-f.cdf(x,2*tmp,2*self.gammaBeta[0]);
			else:
				tmp = self.gammaAlpha[0]*inspItvl**self.gammaAlpha[1];
				deltaDeg1 = degToU - degFrom;
				x1 = deltaDeg1*self.gammaBeta[0]*self.gammaBeta[1]/(tmp);
				prob1 = f.cdf(x1, 2*tmp,2*self.gammaBeta[0]);
				deltaDeg2 = degToL - degFrom;
				x2 = deltaDeg2*self.gammaBeta[0]*self.gammaBeta[1]/(tmp);
				prob2 = f.cdf(x2, 2*tmp,2*self.gammaBeta[0]);
				prob = prob1 - prob2;	
			return prob;
		if (0):
			#no random effect model
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
	
	def update_currentToFail(self):
		self.currentToFail = self.transProb(self.initState, self.nStates-1, self.inspItvl);
		
	def __init__(self, idx, gam_a, gam_b, states, S,  \
					initState, initAge, inspInterval, cCM, cPM, cS):
		self.idx = idx;
		self.gammaAlpha = gam_a;
		self.gammaBeta = gam_b;
		self.nStates = states;		# 0 ... nStates - 1. nStates - 1 is failure states.
		self.failTsh = S;	#failure threshold
		self.initState = initState;
		self.initAge = initAge;
		self.inspItvl = inspInterval;
		self.deg = 0;
		#self.crtState = initState;
		#self.crtDgLvRange = self.state2lv();
		self.cCM = cCM;
		self.cPM = cPM;
		self.cS = cS;
		self.currentToFail = self.transProb(self.initState, self.nStates-1, self.inspItvl);
		self.newToFail =  self.transProb(0, self.nStates-1, self.inspItvl);

		
		
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
		#output
		self.N0 = [];
		self.N1 = [];
		self.Nu = [];
		self.Nu_0 = [];
		self.Nu_1 = [];
		self.time = 0;
		self.objValue = 0;
		self.iterInfo = [];
		