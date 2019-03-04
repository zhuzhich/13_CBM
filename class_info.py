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
			
	def __init__(self, idx, gam_a, gam_b, states, S, \
					initState, inspInterval, cCM, cPM, cS):
		self.idx = idx;
		self.gammaAlpha = gam_a;
		self.gammaBeta = gam_b;
		self.nStates = states;		# 0 ... nStates - 1. nStates - 1 is failure states.
		self.failTsh = S;	#failure threshold
		self.initState = initState;
		self.inspItvl = inspInterval;
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
		self.N0 = 0;
		self.N1 = 0;
		self.Nu = 0;
		self.Nu_0 = [];
		self.Nu_1 = [];
		self.time = 0;
		self.objValue = [];
		self.iterInfo = 0;
		