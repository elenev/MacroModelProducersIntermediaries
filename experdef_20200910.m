% mean-reverting AR(1) log(TFP)
zA= 0.0; 
sig2_A=0.023^2; 
ZA=exp(zA);
rho_A=0.4; 
nu=1;
mu_ZA=exp(zA+.5*sig2_A/(1-rho_A^2));
% deterministic trend labor-augmenting technological progress
g=0.0; 
mu_G=exp(g); 

% Omega transition matrix
Omprob_recess=[0.91, 0.09;
                0.20, 0.80];
Omprob_recess=Omprob_recess(:);    
Omprob_boom=[0.91, 0.09;
                0.20, 0.80];
Omprob_boom=Omprob_boom(:);      
crisis_eligibility='below_trend'; % either 'below_trend' or 'recession'

% borrower
alpha=0.71;
delta=0.937;
betaB=0.94;%0.925;
%betaB=0.94;%0.925;
nuB=1; % Has to be 1 for the effort model to work
sigmaB=1;
deltaK=0.0825;
mu_om=1;
sig2_om=[0.10,0.18].^2;
zeta=[.60,.60]; 
eta=0.2;
outputDWL=1;
gammaB=1-0.622;
PB=1-0.766;
psi=2.0; 
pibar=0.004; 
Phi=0.40;
theta=0.582;
F=theta/(1-delta);
Lscale= 0.7216618597; 
%take_Lscale_as_given='from_expdef';      % values are 'no', 'from_expdef', 'from_guess'
take_Lscale_as_given='from_expdef';
%warning('Calibration mode. take_Lscale_as_given = no');
CBS=false;
FLnodef=false;
phi0=.078;
phi1=0;
n0=1.4;

% intermediary
xi=0.93;
sigma_eps=0.019; 
rho=0.0;
sigmaI=7.0; 
zetaI=0.332;
etaI=0.362;
phiI=0.068;
eqfct_intern=true;
fixediv=false;
sigIexp=0.0;

% saver  
betaS=0.9815;
sigmaS=1;
nuS=1/sigmaS; 
kappa=0.00084;  
chi1adj=0.14;
chiLinear=0;
fracStarget=[0.12,0.12];
deltaT=0.92;

% government
gammaG = 0.172;                 % Government public good expenditures
T = 0.028;                      % Transfers
bT = -20;                       % Sensitivity of transfer spending to productivity growth    
bTau = 4.5  ;                      % Sensitivity of labor taxes to productivity growth
bGamma = -2.0;                   % Sensitivity of government consumption to productivity growth
BG_guess=0.5;
tauK=0.20;%0.217;                     % Corporate tax target of 3.41% of GDP
tauPi=tauK;
tauD=0.132;%0.15;  
%taushield=0.200;                % corp tax rate = tauK / (1 - alpha - taushield);
%alphawedge=0.66/alpha;          % labor share of output = alphawedge * alpha
shieldI = 1;					% Strength of intermediary deposit tax shield (1 = full shield, 0 no shield)

tau0 = 0.295;% 0.285;                     % Labor income tax target of 17.30% of GDP
BG1 = 0.1;
BG2 = 1.2;
exponent=3;
tau_min = 1e-3;
tau_max = 0.6;
BG_min = -0.4;
BG_max = 1.5;
BG_open = 0;
BG_foreign=0;
BG_stv=0.7;
bBG=0.8;
taxrule='Smooth'; %alternative: 'Profligacy/Austerity'

w = 0.00; 
xivec = [xi,xi]; 

mktConsI=1; % intermediary constraint in market value

% set params structure
params=struct;
params.betaB=betaB;
params.betaS=betaS;
params.sigmaS=sigmaS;
params.sigmaB=sigmaB;
params.sigmaI=sigmaI;
params.sigIexp=sigIexp;
params.nuB=nuB;
params.nuS=nuS;
params.gammaB=gammaB;
params.psi=psi;
params.alpha=alpha;
params.deltaK=deltaK;
params.kappa=kappa;
params.tauK=tauK;
params.tauPi=tauPi;
params.tauD=tauD;
params.Lbar=ones(2,1);
params.delta=delta;
params.zeta=zeta;
params.eta=eta;
params.zetaI=zetaI;
params.etaI=etaI;
params.g=g;
params.mu_G=mu_G;
params.zA=zA;
params.ZA=ZA;
params.sig2_A=sig2_A;
params.rho_A=rho_A;
params.mu_ZA=mu_ZA;
params.Omprob_recess=Omprob_recess;
params.Omprob_boom=Omprob_boom;
params.crisis_eligibility=crisis_eligibility;
params.PB=PB;
params.mu_om=mu_om;
params.sig2_om=sig2_om;
params.sigma_eps=sigma_eps;
params.rho=rho;
params.tau0=tau0;
params.BG1 = BG1;
params.BG2 = BG2;
params.exponent=exponent;
params.tau_min = tau_min;
params.tau_max = tau_max;
params.BG_min = BG_min;
params.BG_max = BG_max;
params.BG_open=BG_open;
params.BG_foreign=BG_foreign;
params.T=T;
params.bT=bT;
params.bTau=bTau;
params.bGamma=bGamma;
params.gammaG=gammaG;
params.BG_guess=BG_guess;
params.theta=theta;
params.F=F;
params.Phi=Phi;
params.pibar=pibar;
params.Lscale=Lscale;
params.take_Lscale_as_given=take_Lscale_as_given;
params.outputDWL=outputDWL;
params.w=w;
params.xivec=xivec;
params.CBS=CBS;
params.FLnodef=FLnodef;
params.mktConsI=mktConsI;
params.shieldI=shieldI;
params.chi1adj=chi1adj;
params.phi0=phi0;
params.phi1=phi1;
params.n0=n0;
params.chiLinear=chiLinear;
params.BG_stv=BG_stv;
params.bBG=bBG;
params.taxrule = taxrule;
params.phiI=phiI;
params.eqfct_intern=eqfct_intern;
params.fixediv=fixediv;
params.fracStarget=fracStarget;
params.deltaT=deltaT;

% Define grids

% Benchmark Grid
startgrid.label = 'ini0';
startgrid.KBpts=[1.65,1.75, 1.84, 1.98, 2.05, 2.10, 2.26, 2.45];
startgrid.LBpts=linspace(0.15,0.65,8);
startgrid.WIpts=[-0.03,-0.01,0,0.03,0.05,0.1,0.15,0.25,0.3,0.4];
startgrid.BGpts=[0,0.1833,0.4667,0.7500,1.0333,1.3167,1.4000,1.7000];

semifinegrid.label = '';
semifinegrid.KBpts=[1.6, 1.75, 1.84, 1.98, 2.05, 2.10, 2.26, 2.4];
%semifinegrid.LBpts=[0.15,0.225,0.30,0.325,0.35,0.36,0.37,0.38,0.39,0.40,0.425,0.45,0.5,0.65];
semifinegrid.LBpts=[0.225,0.30,0.325,0.35,0.37,0.39,0.40,0.41,0.42,0.43,0.44,0.46,0.47,0.48,0.49,0.5,0.55];
semifinegrid.WIpts=[-0.04,-0.03,-0.02,-0.01,0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.08,0.09,0.10,0.12,0.14,0.16,0.20,0.27];
semifinegrid.BGpts=[0.1833,0.4667,0.7500, 1.0333,1.3167,1.4000,1.7000];


allsaverstartgrid=startgrid;
allsaverstartgrid.KBpts=linspace(2,2.5,8);
allsaverstartgrid.LBpts=linspace(0.4,0.65,8);

allsaverfinegrid=allsaverstartgrid;
allsaverfinegrid.LBpts = [0.39,0.40,0.41,0.42,0.43,0.44,0.46,0.47,0.48,0.49,0.5,0.55,0.6];

safesafestartgrid = startgrid;
%safesafestartgrid.KBpts=linspace(1.7,2.2,8);
safesafestartgrid.LBpts=linspace(0.35,0.55,8);
safesafestartgrid.WIpts=linspace(0,0.15,8);

safesafefinegrid = safesafestartgrid;
safesafefinegrid.LBpts = [0.325,0.35,0.37,0.39,0.40,0.41,0.42,0.43,0.44,0.46,0.47,0.48,0.49,0.5,0.55];
safesafefinegrid.WIpts = [-0.01, 0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.08,0.09,0.10,0.12,0.14,0.16,0.20];

saferiskystartgrid=startgrid;
saferiskystartgrid.KBpts=linspace(1.7,2.3,8);
saferiskystartgrid.LBpts=linspace(0.35,0.55,8);

saferiskyfinegrid=saferiskystartgrid;
saferiskyfinegrid.LBpts = [0.325,0.35,0.37,0.39,0.40,0.41,0.42,0.43,0.44,0.46,0.47,0.48,0.49,0.5,0.55];
saferiskyfinegrid.WIpts = semifinegrid.WIpts;

sigI0startgrid = startgrid;
sigI0startgrid.LBpts=linspace(0.25, 0.65,8);
sigI0startgrid.WIpts=[0,0.03,0.05,0.1,0.15,0.25];

sigI0semifinegrid = semifinegrid;
sigI0semifinegrid.WIpts=[-0.06,-0.02,-0.01,0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.10,0.125,0.15];

sigI20startgrid = startgrid;
sigI20startgrid.WIpts=[0,0.03,0.05,0.1,0.15,0.25,0.3,0.4];
sigI20startgrid.BGpts=[0,0.1833,0.4667,0.7500,1.0333,1.3167,1.4000];

sigI20semifinegrid = semifinegrid;
sigI20semifinegrid.WIpts=[0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.08,0.09,0.10,0.12,0.14,0.16,0.20,0.27];
sigI20semifinegrid.BGpts=[0.1833,0.4667,0.7500, 1.0333,1.3167,1.4000];

% xi85 and xi80 grid
midxistartgrid.KBpts=[1.65,1.75, 1.84, 1.98, 2.05, 2.10, 2.26, 2.45];
midxistartgrid.LBpts=linspace(0.10,0.45,8);
midxistartgrid.WIpts=[0,0.03,0.05,0.1,0.125,0.15,0.25,0.3,0.4];
midxistartgrid.BGpts=[0,0.1833,0.4667,0.7500,1.0333,1.3167,1.4000,1.6000];

midxisemifinegrid.KBpts=[1.65,1.75, 1.84, 1.98, 2.05, 2.10, 2.26, 2.45];
midxisemifinegrid.LBpts=[0.19,0.22,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42];
midxisemifinegrid.WIpts=[0,0.02,0.05,0.06,0.07,0.08,0.085,0.09,0.095,0.10,0.105,0.11,0.115,0.12,0.13,0.14,0.15,0.17];
midxisemifinegrid.BGpts=[0.1833,0.4667,0.7500, 1.0333,1.3167,1.4000,1.6000];

% Low-xi grid
lowxistartgrid.KBpts=[1.65,1.75, 1.84, 1.98, 2.05, 2.10, 2.26, 2.45];
lowxistartgrid.LBpts=linspace(0.07,0.4,8);
lowxistartgrid.WIpts=[0,0.05,0.1,0.125,0.15,0.20,0.25,0.3,0.4];
lowxistartgrid.BGpts=[0,0.1833,0.4667,0.7500,1.0333,1.3167,1.4000,1.6000];

lowxisemifinegrid.KBpts=[1.65,1.75, 1.84, 1.98, 2.05, 2.10, 2.26, 2.45];
lowxisemifinegrid.LBpts=[0.10,0.14,0.16,0.18,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.30,0.32,0.35];
lowxisemifinegrid.WIpts=[0,0.025,0.05,0.08,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.18,0.20,0.22,0.24,0.26,0.30];
lowxisemifinegrid.BGpts=[0.1833,0.4667,0.7500, 1.0333,1.3167,1.4000,1.6000];

% High-xi grid
hixistartgrid = startgrid;
hixistartgrid.BGpts = [0.4667,0.7500,1.0333,1.3167,1.4000,1.6000,2.0,2.4];

hixisemifinegrid = semifinegrid;
hixisemifinegrid.BGpts=[0.4667,0.7500,1.0333,1.3167,1.4000,1.6000,2.0,2.4];

% Kappa grid
kappastartgrid = startgrid;
kappastartgrid.LBpts = lowxistartgrid.LBpts;
kappasemifinegrid = semifinegrid;
kappasemifinegrid.LBpts = lowxistartgrid.LBpts;


% High betaB grid
betaBstartgrid = startgrid;
betaBstartgrid.KBpts = startgrid.KBpts + 0.15;
betaBstartgrid.LBpts = startgrid.LBpts - 0.07;
 
betaBsemifinegrid = semifinegrid;
betaBsemifinegrid.KBpts = semifinegrid.KBpts + 0.15;
betaBsemifinegrid.LBpts = semifinegrid.LBpts - 0.07;

% Low betas grid
betaSstartgrid = startgrid;
betaSstartgrid.KBpts = startgrid.KBpts - 0.2;
betaSstartgrid.LBpts = startgrid.LBpts - 0.07;

betaSsemifinegrid = semifinegrid;
betaSsemifinegrid.KBpts = semifinegrid.KBpts - 0.2;
betaSsemifinegrid.LBpts = semifinegrid.LBpts - 0.07;

% Assemble grids
benchgrid = struct('startgrid',startgrid,'finegrid',semifinegrid);
allsavergrid = struct('startgrid',allsaverstartgrid,'finegrid',allsaverfinegrid);
safesafegrid = struct('startgrid',safesafestartgrid,'finegrid',safesafefinegrid);
saferiskygrid = struct('startgrid',saferiskystartgrid,'finegrid',saferiskyfinegrid);
lowxigrid = struct('startgrid',lowxistartgrid,'finegrid',lowxisemifinegrid);
midxigrid = struct('startgrid',midxistartgrid,'finegrid',midxisemifinegrid);
hixigrid = struct('startgrid',hixistartgrid,'finegrid',hixisemifinegrid);
betaBgrid = struct('startgrid',betaBstartgrid,'finegrid',betaBsemifinegrid);
betaSgrid = struct('startgrid',betaSstartgrid,'finegrid',betaSsemifinegrid);
kappagrid = struct('startgrid',kappastartgrid,'finegrid',kappasemifinegrid);
%auggrid = struct('startgrid',augstartgrid,'finegrid',augsemifinegrid);
sigI0grid = struct('startgrid',sigI0startgrid,'finegrid',sigI0semifinegrid);
sigI20grid = struct('startgrid',sigI20startgrid,'finegrid',sigI20semifinegrid);

grids = struct('benchgrid',benchgrid, ...
    'lowxigrid',lowxigrid, ...
    'midxigrid',midxigrid, ...
	'hixigrid',hixigrid, ...
	'kappagrid',kappagrid,...
    'betaBgrid',betaBgrid, ...
    'betaSgrid',betaSgrid,...
    'allsavergrid',allsavergrid,...
    'safesafegrid',safesafegrid,...
    'saferiskygrid',saferiskygrid,...
    'sigI0grid',sigI0grid,...
	'sigI20grid',sigI20grid);

% Define experiments as modifications to base
expers_macropru = {'xpbase',{'take_Lscale_as_given','no'};
                'bench',{};
				'xi70',{'xivec',[0.70,0.70];'grid','lowxigrid'};
                'xi75',{'xivec',[0.75,0.75];'grid','lowxigrid'};
                'xi80',{'xivec',[0.80,0.80];'grid','midxigrid'};
                'xi85',{'xivec',[0.85,0.85];'grid','midxigrid'};
                'xi90',{'xivec',[0.90,0.90]};
				'xi92',{'xivec',[0.92,0.92]};
				'xi94',{'xivec',[0.94,0.94]};
				'xi96',{'xivec',[0.96,0.96]};
                'xi975',{'xivec',[0.975,0.975];'grid','hixigrid'};
				'kappa01',{'kappa',0.01;'grid','kappagrid'};
                'xi9395',{'xivec',[0.93,0.95]};
				'xi81midxigrid',{'xivec',[0.81,0.81];'grid','midxigrid'};
                'xi82midxigrid',{'xivec',[0.82,0.82];'grid','midxigrid'};
                'xi83midxigrid',{'xivec',[0.83,0.83];'grid','midxigrid'};
                'xi84midxigrid',{'xivec',[0.84,0.84];'grid','midxigrid'};
                'xi86midxigrid',{'xivec',[0.86,0.86];'grid','midxigrid'};
                'xi87midxigrid',{'xivec',[0.87,0.87];'grid','midxigrid'};
                'xi88midxigrid',{'xivec',[0.88,0.88];'grid','midxigrid'};
                'xi89midxigrid',{'xivec',[0.89,0.89];'grid','midxigrid'};
                'xi91benchgrid',{'xivec',[0.91,0.91];'grid','benchgrid'};
                'xi93benchgrid',{'xivec',[0.93,0.93];'grid','benchgrid'};
                'xi95benchgrid',{'xivec',[0.95,0.95];'grid','benchgrid'};
                'xi97benchgrid',{'xivec',[0.97,0.97];'grid','hixigrid'};
				'xi9294',{'xivec',[0.92,0.94]};
				'xi9195',{'xivec',[0.91,0.95]};
				'xi9295',{'xivec',[0.92,0.95]};
				'FLbench',{'CBS',true;'chi1adj',0;'FLnodef',true;'pibar',0;'grid','allsavergrid';};
				'FLwithdef',{'CBS',true;'chi1adj',0;'grid','allsavergrid';};
                'benchsafesafe',{'FLnodef',true;'rho',100;'grid','safesafegrid'};
                'benchsafesafesigI0',{'FLnodef',true;'rho',100;'sigmaI',0;'grid','safesafegrid'};
                'benchsaferisky',{'FLnodef',true};
				'benchriskysafesigI0',{'sigmaI',0;'rho',100;'grid','sigI0grid'};
				'benchriskysafe',{'rho',100;};
                'benchriskierrisky',{'phi1',1;};
				'sigI0',{'sigmaI',0;'grid','sigI0grid'};
				'sigI8',{'sigmaI',8};
				'sigI20rho100',{'sigmaI',20;'rho',100;'grid','sigI20grid'};
				'sigI8xi91',{'sigmaI',8;'xivec',[0.91,0.91]};
				'sigI8xi95',{'sigmaI',8;'xivec',[0.95,0.95]};
				'sigI0xi91',{'sigmaI',0;'xivec',[0.91,0.91];'grid','sigI0grid'};
				'sigI0xi95',{'sigmaI',0;'xivec',[0.95,0.95];'grid','sigI0grid'};
				'benchriskysafexi91',{'rho',100;'xivec',[0.91,0.91]};
				'benchriskysafexi95',{'rho',100;'xivec',[0.95,0.95]};
				'noshieldsameratexi91',{'shieldI',0;'xivec',[0.91,0.91]};
				'noshieldsamerate',{'shieldI',0};
				'noshieldsameratexi95',{'shieldI',0;'xivec',[0.95,0.95]};
				'betaBupxi91',{'betaB',params.betaB+0.0015;'xivec',[0.91,0.91]};
				'betaBup',{'betaB',params.betaB+0.0015};
				'betaBupxi95',{'betaB',params.betaB+0.0015;'xivec',[0.95,0.95]};
				'betaSdownxi91',{'betaS',params.betaS-0.0015;'xivec',[0.91,0.91]};
				'betaSdown',{'betaS',params.betaS-0.0015};
				'betaSdownxi95',{'betaS',params.betaS-0.0015;'xivec',[0.95,0.95]};
				'benchriskierriskyxi91',{'phi1',1;'xivec',[0.91,0.91]};
				'benchriskierriskyxi95',{'phi1',1;'xivec',[0.95,0.95]};
				'altshares',{'PB',1-0.550;'gammaB',1-0.308};
				'altsharesxi91',{'PB',1-0.550;'gammaB',1-0.308;'xivec',[0.91,0.91]};
				'altsharesxi95',{'PB',1-0.550;'gammaB',1-0.308;'xivec',[0.95,0.95]};
				'nosavers',{'chiLinear',50;'chi1adj',0};
				'nosaversxi91',{'chiLinear',50;'chi1adj',0;'xivec',[0.91,0.91]};
				'nosaversxi95',{'chiLinear',50;'chi1adj',0;'xivec',[0.95,0.95]};
				'sig10',{'sigmaB',10;'sigmaS',10};
				'sig10FLbench',{'sigmaB',10;'sigmaS',10;'CBS',true;'chi1adj',0;'FLnodef',true;'pibar',0;'grid','allsavergrid';};
				'sig10FLwithdef',{'sigmaB',10;'sigmaS',10;'CBS',true;'chi1adj',0;'grid','allsavergrid';};
                'sig10benchsafesafe',{'sigmaB',10;'sigmaS',10;'FLnodef',true;'rho',100;'grid','safesafegrid'};
                'sig10benchsafesafesigI0',{'sigmaB',10;'sigmaS',10;'FLnodef',true;'rho',100;'sigmaI',0;'grid','safesafegrid'};
                'sig10benchsaferisky',{'sigmaB',10;'sigmaS',10;'FLnodef',true};
				'sig10benchriskysafesigI0',{'sigmaB',10;'sigmaS',10;'sigmaI',0;'rho',100;'grid','sigI0grid'};
				'sig10benchriskysafe',{'sigmaB',10;'sigmaS',10;'rho',100;};
                'sig10benchriskierrisky',{'sigmaB',10;'sigmaS',10;'phi1',1;};
				'sig10sigI0',{'sigmaB',10;'sigmaS',10;'sigmaI',0;'grid','sigI0grid'};
				'kappa001',{'kappa',0.001};
				'kappa003',{'kappa',0.003};
				'kappa005',{'kappa',0.005};
				'kappa005midxigrid',{'kappa',0.005;'grid','midxigrid'};
				'basesharethresh05',{'PB',0.261;'gammaB',0.425};
				'basesharethresh01',{'PB',0.289;'gammaB',0.472};
				'basesharethresh15',{'PB',0.213;'gammaB',0.342};
				'basesharethresh20',{'PB',0.195;'gammaB',0.311};
				'basesharethresh05xi91',{'PB',0.261;'gammaB',0.425;'xivec',[0.91,0.91]};
				'basesharethresh01xi91',{'PB',0.289;'gammaB',0.472;'xivec',[0.91,0.91]};
				'basesharethresh15xi91',{'PB',0.213;'gammaB',0.342;'xivec',[0.91,0.91]};
				'basesharethresh20xi91',{'PB',0.195;'gammaB',0.311;'xivec',[0.91,0.91]};
				'basesharethresh05xi95',{'PB',0.261;'gammaB',0.425;'xivec',[0.95,0.95]};
				'basesharethresh01xi95',{'PB',0.289;'gammaB',0.472;'xivec',[0.95,0.95]};
				'basesharethresh15xi95',{'PB',0.213;'gammaB',0.342;'xivec',[0.95,0.95]};
				'basesharethresh20xi95',{'PB',0.195;'gammaB',0.311;'xivec',[0.95,0.95]};
				};

% Cartesian Products of Parameters (for calibration)
calib = cell(0,1);	% sets of cross-products

% Insert calibration definitions here e.g.
% calib{1}.parameter = {value1, value2, value3}

expers_calib = cell( 500, 2);
counter=1;
for ii=1:numel(calib)
	param_names = fieldnames(calib{ii});
	N_params = numel(param_names);
	idxCalibParamsMesh = cell( N_params, 1);
	calibParamsVec = struct2cell(calib{ii});
	idxcell = cellfun(@(x)1:numel(x),calibParamsVec,'UniformOutput',false)';
	[idxCalibParamsMesh{:}] = ndgrid( idxcell{:} );
	idxCalibParamsFlat = cell2mat( ...
		cellfun(@(c)c(:),idxCalibParamsMesh,'UniformOutput',false)' ...
		);
	calibParamsFlat = cell(size(idxCalibParamsFlat));
	for kk=1:N_params
		calibParamsFlat(:,kk) = calib{ii}.(param_names{kk})(idxCalibParamsFlat(:,kk));
	end
	for jj=1:size(calibParamsFlat,1)
		calibName = sprintf('calib%02d',counter);
		calibCell = cell( size(calibParamsFlat,2) + 1, 2);
		calibCell(1:end-1,1) = fieldnames(calib{ii});
		calibCell(1:end-1,2) = calibParamsFlat(jj,:)';
		calibCell{end,1} = 'take_Lscale_as_given';
		calibCell{end,2} = 'no';
		expers_calib{counter,1} = calibName;
		expers_calib{counter,2} = calibCell;
		counter=counter+1;
	end
end
expers_calib = expers_calib(1:counter-1,:);

% Gentskow-Shapiro Parameter (Semi) Elasticities
%vary_params = {'alpha','bBG','betaB','betaS','bGamma','bT','bTau', ...
%	'fracStarget','chi1adj','deltaK','n0','Phi','pibar','phi1','psi','rho_A','sig2_A', ...
%	'sig2_om','sigma_eps','sigmaI','zeta'};
%vary_params = {'Phi','phi1','psi','pibar','deltaK','sigmaI','n0'};
vary_params = {'rho_A','sig2_A','sig2_om','psi','alpha','deltaK','pibar','zeta', ...
	'Phi','sigma_eps', 'sigmaI','fracStarget','chi1adj','betaB','betaS'};
adjust_Lscale = true;
eps_default = 1e-3;
% Modify interval size for specific parameters
eps_struct = struct; 
eps_struct.sig2_A=1e-6;
eps_struct.sig2_om=1e-6;
% Compute forward instead of central differences for some parameters 
% (e.g. if baseline value is at the lower bound)
forwardDifferenceOnly = {}; 

% Compute full elasticities (i.e. percent parameter changes)
fullElasticities = {};

% Vary all elements of a vector-valued parameter at once
linked = {'fracStarget','zeta'};

expers_gs=cell(0,6);
for ii=1:numel(vary_params)
	currParam = vary_params{ii};
	currValue = params.(currParam);
	isContinuous = all( ~isa(currValue,'logical') & ~isa(currValue,'char') );
	N_paramDim = length( params.(currParam) );
	if N_paramDim==1 && isContinuous
		dimLabel = {''};
	elseif N_paramDim==2 && isContinuous
		dimLabel = {'lo','hi'};
	else
		error('Cannot do GS sensitivity for this parameter.');
	end
	
	try
		currEps = eps_struct.(currParam);
	catch
		currEps = eps_default;
	end
	
	for jj=1:N_paramDim
		newval_U = currValue;
		newval_D = currValue;
		if ismember(currParam,linked)
			idx = 1:N_paramDim;
		else
			idx = jj;
		end
		if ismember(currParam,fullElasticities)
			newval_U(idx)=newval_U(idx)*exp(currEps);
			newval_D(idx)=newval_D(idx)*exp(-currEps);
		else
			newval_U(idx)=newval_U(idx)+currEps;
			newval_D(idx)=newval_D(idx)-currEps;
		end
		
		expername_U = sprintf('GS%s%sU',strrep(currParam,'_',''),dimLabel{jj});
		experdevs_U = {currParam,newval_U};
		expername_D = sprintf('GS%s%sD',strrep(currParam,'_',''),dimLabel{jj});
		experdevs_D = {currParam,newval_D};
		if adjust_Lscale
			experdevs_U = [ experdevs_U; {'take_Lscale_as_given', 'no'} ];
			experdevs_D = [ experdevs_D; {'take_Lscale_as_given', 'no'} ];
		end
		expers_gs = [ expers_gs; {expername_U,experdevs_U,currParam,jj,dimLabel{jj},'U'} ];
		if ~ismember(currParam,forwardDifferenceOnly)
			expers_gs = [ expers_gs; {expername_D,experdevs_D,currParam,jj,dimLabel{jj},'D'} ];
		end
		
		if ismember(currParam,linked)
			break;
		end
	end
end

expers = [expers_macropru; expers_calib; expers_gs(:,1:2)];

% % Relax Lscale for calibration economies (not here, comparison with julia
% code, keep take_Lscale_as_given='from_expdef')
% CALIBRATION_IDX=[2,6:size(expers,1)];
% expers(CALIBRATION_IDX,:) = ...
%     [expers(CALIBRATION_IDX,1), ...
%      cellfun(@(changes)[changes;{'take_Lscale_as_given','no'}], ...
%         expers(CALIBRATION_IDX,2),'UniformOutput',false)];

%% Create each experiment definition (fully dynamic below this line)
N_EXP=size(expers,1);
N_GRIDS = length(fieldnames(benchgrid));

% Write list of defined experiments to file
fid = fopen([mfilename,'.txt'],'w');
for i=1:N_EXP
   fprintf(fid,expers{i,1});
   fprintf(fid,'\n');
end
fclose(fid);

% Create experiment table
baserow = {'test','ini0',{startgrid.KBpts},{startgrid.LBpts},{startgrid.WIpts},{startgrid.BGpts}};
basetable = cell2table(baserow);
basetable.Properties.VariableNames = {'exper','grid','KBpts','LBpts','WIpts','BGpts'};
basetable = [basetable, struct2table(params,'AsArray',true)];

expertable = repmat(basetable,N_EXP*N_GRIDS,1);

benchgrid_cell=struct2cell(benchgrid);
for i=1:N_GRIDS
    fnames = fieldnames(benchgrid_cell{i});
    gridvectors=struct2cell(benchgrid_cell{i})';
    [~,~,order]=intersect(expertable.Properties.VariableNames(3:6),fnames, ...
            'stable');
    expertable{N_EXP*(i-1)+1 : N_EXP*i,3:6} = ...
        repmat( gridvectors(order), N_EXP, 1);
    expertable(N_EXP*(i-1)+1 : N_EXP*i,2) = ...
        {benchgrid_cell{i}.label};   
end

gridOrder = @(grid_as_cell,order)grid_as_cell(order);
for i=1:N_EXP
   expertable.exper([i,N_EXP+i]) = expers(i,1);
   changes = expers{i,2};
   for j=1:size(changes,1)
      varname=changes{j,1};
      if strcmp(varname,'grid')
          newgrid = struct2cell(grids.(changes{j,2}));
          fnames = fieldnames(newgrid{1});
          [~,~,order]=intersect(expertable.Properties.VariableNames(3:6),fnames, ...
            'stable');
        for k=1:N_GRIDS
          expertable{i+N_EXP*(k-1),3:6} = ...
              gridOrder(struct2cell(newgrid{k}),order)';
        end
      else
          if ischar(changes{j,2})
            expertable.(varname)([i,N_EXP+i],:)=repmat(changes(j,2),N_GRIDS,1);
          else
            expertable.(varname)([i,N_EXP+i],:)=repmat(changes{j,2},N_GRIDS,1);    
          end
      end
   end
end

% Name each row
row_id = cell(2*N_EXP,1);
idx_suffix = find(~strcmp(expertable.grid,''));
idx_nosuffix = find(strcmp(expertable.grid,''));
row_id(idx_suffix) = arrayfun(@(idx){[expertable.exper{idx},'_',expertable.grid{idx}]},idx_suffix,'UniformOutput',false);
row_id(idx_nosuffix) = arrayfun(@(idx)expertable.exper(idx),idx_nosuffix,'UniformOutput',false);
expertable.Properties.RowNames = [row_id{:}]';
%expertable.Properties.RowNames(1) = {'xpbase'};

% Create a struct of structs
allexpers=cell2struct(expertable.Properties.RowNames,expertable.Properties.RowNames);
for i=1:size(expertable,1)
    s = table2struct(expertable(i,3:6));
    s.params = table2struct(expertable(i,7:end));
    name = expertable.Properties.RowNames(i);
    name = name{:};
    allexpers.(name) = s;
end
