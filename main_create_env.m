if usejava('desktop')
   clear; 
end

% ===========================================
% program control
% ===========================================

% name of file with experiment definitions
if ~exist('experdef_file','var')
    experdef_file='experdef_20200910';
end
% name of experiment
if ~exist('expername','var')
    expername='xpbase_ini0'; %bench_ini0
end

% possible values 'no_guess', 'guess'
if ~exist('guess_mode','var')
    guess_mode='no_guess';
end
% path to file with initial guess; not used if guess_mode='no_guess'
if ~exist('guess_path','var')
    guess_path='res_20200215_GMSbankDWLetabank3_i50';
end

% path and file name for output file
if ~exist('outfname','var')
    outfname=['env_',expername,'.mat'];
end

% approximation mode
approxMode='linear';

% approximation mode
if ~exist('approxMode','var')
    approxMode='linear';
end

% Use analytic Jacobian
% approximation mode
if ~exist('useJacobian','var')
    useJacobian=true;
end

% only compute steady state (without writing initial model object)
if ~exist('ststonly','var')
    ststonly=0;
end

% ===================================================
% non-stochastic steady state
% ===================================================

% load experiment definitions
run(experdef_file);
disp(experdef_file);
expdef=allexpers.(expername);
params=expdef.params;

% If some parameters are dynamically overwritten during a batch job, read
% them in
controlFile='control.mat';
try
    load(controlFile);
catch
end

if expdef.params.CBS
    modelclassname='FLProductionModel';
else
    modelclassname='ProductionModel';
end
ststfun=str2func([modelclassname,'.compStSt']);
instfun=str2func(modelclassname);
guessfun=str2func([modelclassname,'.assignGuess']);
jacfun=str2func([modelclassname,'.constructJacobian']);

% compute additional parameters
% expdef.params.mu_G=exp(expdef.params.g+0.5*expdef.params.sig2_x);
jensen_corr_T=exp(expdef.params.bT^2/2*expdef.params.sig2_A);
expdef.params.T = expdef.params.T / jensen_corr_T;
jensen_corr_tau=exp(expdef.params.bTau^2/2*expdef.params.sig2_A);
expdef.params.tau0 = expdef.params.tau0 / jensen_corr_tau;
jensen_corr_gammaG=exp(expdef.params.bGamma^2/2*expdef.params.sig2_A);
expdef.params.gammaG = expdef.params.gammaG / jensen_corr_gammaG;

expdef.params.gammaS=1-expdef.params.gammaB;
expdef.params.PS=1-expdef.params.PB;

expdef.params.Pbar=[expdef.params.PB,expdef.params.PS]';

k1=(expdef.params.tau0-expdef.params.tau_min)/(expdef.params.BG1-expdef.params.BG_min)^expdef.params.exponent;
k2=(expdef.params.tau_max-expdef.params.tau0)/(expdef.params.BG_max-expdef.params.BG2)^expdef.params.exponent;

% tax rule for stationary gov debt
if strcmp(expdef.params.taxrule,'Profligacy/Austerity')
	expdef.params.tau = @(BG) min(expdef.params.tau_max, max( expdef.params.tau_min, ...
                            expdef.params.tau0 ...
                            + ( -(BG<expdef.params.BG1) .* k1 .* (expdef.params.BG1-BG).^expdef.params.exponent ) ...
                            +    (BG>expdef.params.BG2) .* k2 .* (BG-expdef.params.BG2).^expdef.params.exponent ) );

	expdef.params.dTau = @(BG) (BG<BG1 & BG>BG_min) .* k1 .* exponent .* (BG1-BG).^ (exponent-1) ...
	    + (BG>BG2 & BG<BG_max) .* k2 .* exponent .* (BG-BG2).^(exponent-1);
elseif strcmp(expdef.params.taxrule,'Smooth')
	expdef.params.tau = @(BG) min(params.tau_max, max( params.tau_min, real( params.tau0 * ...
	 ( (BG - params.BG_min) ./ (params.BG_stv - params.BG_min) ).^params.bBG ) ) );
	expdef.params.dTau = @(BG) (BG > params.BG_min) .* params.tau0 * params.bBG/(params.BG_stv-params.BG_min) .* ...
		((BG -  params.BG_min) ./ (params.BG_stv - params.BG_min) ).^(params.bBG-1) + (BG <= params.BG_min) * 0;
end

% Set bankDWL based on values of zetaI, zeta, and fracStarget
expdef.params.bankDWL = 1 - expdef.params.fracStarget(1);
                        
% compute steady state values
% do this once for each state of sigma_omega
options=optimset('Display','iter','TolX',1e-10,'TolFun',1e-10,'MaxIter',100);


% gvec={[ -2.2303    0.8717   -0.9831   -0.5776    2.3073  log(0.01/(1-0.01)) log(0.8) 0.1993 ], ...
% 	[ -3.1987    0.8529   -1.0713   -0.6658    2.3347  0.0039 log(0.8) -0.8376 ]};

% gvec={
%     [-2.4352    1.0052   -0.9083   -0.5029    2.3913    0.0018   -2.6808   -2.0370   ], ...
% 	[ -3.1987    0.8529   -1.0713   -0.6658    2.3347 .01 .07 0.0039   ]};

if expdef.params.CBS
    gvec={
        [-1.9965    0.9473   -1.0560   -0.4807    2.3113     0.2, 0.1],...
        [-2.9271    0.9517   -1.1341   -0.5588    2.3856     0.2, 0.1]
        };
else
    gvec={
        [-1.9965    0.9473   -1.0560   -0.4807    2.3113    0.0007   -2.6019   -0.2, 0.1, -2.25],...
        [-2.9271    0.9517   -1.1341   -0.5588    2.3856   -0.0005   -3.0096   -0.2, 0.1, -2.25]
        };
    if expdef.params.chiLinear>0
        gvec{1} = gvec{1}(1:end-1);
        gvec{2} = gvec{2}(1:end-1);
    end
    if expdef.params.FLnodef
       gvec{1}(8)=0.2;
       gvec{2}(8)=0.2;       
    end
end


if strcmp(expdef.params.take_Lscale_as_given,'no')
    gvec{1}=[gvec{1},-0.3266];
	gvec{2}=[gvec{2},-0.3183];
elseif strcmp(expdef.params.take_Lscale_as_given,'from_guess') && strcmp(guess_mode,'guess')
    guessobj=load(guess_path,'mobj');
    expdef.params.Lscale=guessobj.mobj.Params.Lscale;
    clear 'guessobj';
end

stv=cell(2,1);
solout = cell(2,1);
for i=1:2
	fh_compStSt=@(x)ststfun(x,expdef.params,i,0);
	[solvec,~,exfl]=fsolve(fh_compStSt,gvec{i},options);
	if exfl<1
		disp('!! Problem computing steady state');
	end
    [~,stvtmp]=ststfun(solvec,expdef.params,i,1);
	stv{i}=stvtmp;
    solout{i}=solvec;
    %gvec{i+1}=solvec;
end

expdef.params.Lscale=stv{1}.statsout.Lscale;
expdef.params.AStarget=cellfun(@(s)s.Add.AStarget,stv)';
if ~expdef.params.CBS
    expdef.params.WIbar=stv{1}.State.WI; % for dividend target with fixediv=true
end

if ststonly
	
	firmLeverage = @(i)stv{i}.Sol.qB * stv{i}.Add.AB / stv{i}.State.KB;
	capitalGDP = @(i)stv{i}.State.KB;
	
	fprintf('Firm leverage: %0.3f -- %0.3f \n', firmLeverage(1), firmLeverage(2));
	fprintf('Capital / GDP: %0.3f -- %0.3f \n', capitalGDP(1), capitalGDP(2));
    
    if ~expdef.params.CBS
        bankDefault = @(i)1-stv{i}.Add.F_eps;
        saverBondShare = @(i)stv{i}.Sol.AS_g / stv{1}.Add.AB_g;
        fprintf('Bank Default:  %0.3f -- %0.3f \n', bankDefault(1), bankDefault(2));
        fprintf('S Bond Share:  %0.3f -- %0.3f \n', saverBondShare(1), saverBondShare(2));
    end
    
	try
		prev_ss=load(['ssonly_',outfname]);
		rescell = prev_ss.rescell;
		rescell{end+1} = stv;
	catch
		rescell{1} = stv;
	end
	save(['ssonly_',outfname],'rescell');
    return;
end


% ===================================================
% Parameters of stochastic model
% ===================================================

% for mean-reverting AR(1) (log)TFP, denoted x here
N_x=5;
try
    sig2_x=expdef.params.sig2_A;
catch
    sig2_x=0.039^2;
end
mu_x=expdef.params.zA;
mu_ZA=expdef.params.mu_ZA;
rho_x=expdef.params.rho_A;
Skew_x=0;
[Xprob,X] = model.DSGEModel.rouwen(rho_x,mu_x,sqrt(sig2_x),Skew_x,N_x);
X=exp(X); % exponentiated to get TFP from log(TFP)
Xprob=Xprob';

% for shock to std.dev. of borrowers house values
threshold=0;
N_om=2;
Om=expdef.params.sig2_om';
% for countercyclical intermediary capital requirements
N_xi=1; 

try
    Omprob_recess=reshape(expdef.params.Omprob_recess,2,2);
    Omprob_boom=reshape(expdef.params.Omprob_boom,2,2);

    % total Markov transition for exogenous states
    % first find recession states 

    if strcmp(expdef.params.crisis_eligibility,'recession')
        threshold=0;
    elseif strcmp(expdef.params.crisis_eligibility,'below_trend')
        threshold=zA;
    end
catch
    Omprob_recess=[0.8, 0.2;
            0.3, 0.7];
    Omprob_boom=[0.9, 0.1;
            0.5, 0.5];
end
recess=find(X<1+threshold,1,'last');
boom=N_x-recess;
omcomb=grid.StateSpaceGrid.makeCombinations([N_x,N_om,N_xi]);
omcomb(:,3)=omcomb(:,2);
% mpts_perm=[ X(omcomb(:,1)), Om(omcomb(:,2)), expdef.params.zeta(omcomb(:,2))', expdef.params.xivec(1)*ones(size(omcomb,1),1) ];
mpts_perm=[ X(omcomb(:,1)), Om(omcomb(:,2)), expdef.params.zeta(omcomb(:,2))', expdef.params.xivec(omcomb(:,3))', expdef.params.AStarget(omcomb(:,2))' ];

% transition matrix
rectrans=kron(Xprob(1:recess,:),Omprob_recess);
boomtrans=kron(Xprob(recess+1:end,:),Omprob_boom);
mtrans=[ rectrans; boomtrans ];
exnpt=size(mpts_perm,1);
mpts_all=[mpts_perm, (1:exnpt)'];
     
% simulate exog process and produce histogram
Nsim=100000;
simst=model.DSGEModel.hitm_s(mtrans',rand(Nsim,1));
simfrac=histc(simst,1:exnpt)/Nsim;
disp('-----------------------------------------');
disp('Unconditional prob of each state: ');
disp(num2str([(1:exnpt)',mpts_perm,simfrac]));
disp(['"Recession" states: ', num2str(sum(simfrac(1:N_om*recess)))]); 
omhi_ind=N_om*(1:N_x);
disp(['High uncertainty states: ', num2str(sum(simfrac(omhi_ind)))]);
% average length of spell
forecl=ismember(simst,omhi_ind);
forestend=forecl(2:end)-forecl(1:end-1);
firstend=find(forestend==-1,1,'first');
laststart=find(forestend==1,1,'last');
forestend=forestend(firstend+1:laststart-1);
forecep=sum(forestend==1);
disp(['Uncertainty episodes (%): ',num2str(forecep/Nsim)]);
forest=find(forestend==1);
foreend=find(forestend==-1);
forelen=foreend-forest;
disp(['Average length: ',num2str(mean(forelen))]);
% SD of interquartile range
sdl = expdef.params.sig2_om(1)^0.5;
sdh = expdef.params.sig2_om(2)^0.5;
IQ_omlo= iqr(makedist('Normal',expdef.params.mu_om,sdl));
IQ_omhi= iqr(makedist('Normal',expdef.params.mu_om,sdh));
IQ_sd=std((1-forecl)*IQ_omlo + forecl*IQ_omhi);
SD_sd = std((1-forecl)*sdl + forecl*sdh);
disp(['SD of IQ range: ', num2str(IQ_sd)]); 
disp(['SD of SD: ', num2str(SD_sd)]); 
% Effective tax and spending rates
ZAvec=mpts_perm(simst,1);
ZA_gr=log(ZAvec(2:end))-log(ZAvec(1:end-1));
ZA_factor=min(1, ZAvec/mu_ZA);
ptax_eff=0.91*expdef.params.alpha*expdef.params.tau0*ZA_factor.^expdef.params.bTau;
bptax=regress(log(ptax_eff(2:end)),[ZA_gr, ones(Nsim-1,1)]);
trans_eff=expdef.params.T*ZA_factor.^expdef.params.bT;
btrans=regress(log(trans_eff(2:end)),[ZA_gr, ones(Nsim-1,1)]);
gamma_eff=expdef.params.gammaG*ZA_factor.^expdef.params.bGamma;
bgamma=regress(log(gamma_eff(2:end)),[ZA_gr, ones(Nsim-1,1)]);

disp(['Eff. personal tax/Y (m, beta): ', num2str(mean(ptax_eff)), ', ', num2str(bptax(1))]); 
disp(['Eff.  transfers/Y (m, beta): ', num2str(mean(trans_eff)), ', ', num2str(btrans(1))]); 
disp(['Eff. discr. spending/Y (m,beta): ', num2str(mean(gamma_eff)), ', ', num2str(bgamma(1))]); 


% variable lists
% exogenous state variable names
exogenv=struct;
exogenv.exnames={'ZA','Om','Zeta','xi','chi0'};
exogenv.exidx=struct('ZA',1,'Om',2,'Zeta',3,'xi',4);
exogenv.exnpt=exnpt;
exogenv.pts_perm=mpts_perm;
exogenv.pts_all=mpts_all;
exogenv.mtrans=mtrans;
exogenv.normparam=[mu_x, sig2_x];

          
% ===========================================
% create model object
% ===========================================

% endogenous state variables
endogenv=struct;
if expdef.params.CBS
    endogenv.ennames={'KB','LB','BG'};
else
    endogenv.ennames={'KB','LB','WI','BG'};
end        
                    
% assign values for initial guess
[solguessvec,Vguessvec,Vnames]=guessfun(stv{1});
            
% save name lists
endogenv.solnames=fieldnames(stv{1}.Sol);
endogenv.solbase=solguessvec;
endogenv.addnames=fieldnames(stv{1}.Add);
endogenv.Vnames=Vnames;

if expdef.params.CBS
    endogenv.condnames={'expRP','stdRP', ...
        'expRPB','stdRPB','effect_RfB','effect_CovPB','effect_collatPB','expRPB_check', ...
        'expRB','stdRB','expOmA','expCBgr','expCSgr','stdCBgr','stdCSgr',...
        'expEROWB','SRB'};
else
    endogenv.condnames={'expRP','stdRP','effect_Rf','effect_CovP','effect_collatP','expRP_check', ...
        'expRPB','stdRPB','effect_RfB','effect_CovPB','effect_collatPB','expRPB_check', ...
        'expRB','stdRB','expOmA','expCBgr','expCIgr','expCSgr','stdCBgr','stdCIgr','stdCSgr',...
        'expEROWB','SRB','expEROWI', 'SRI', 'retBE'};
end


% initial definition of grids based on steady state values of state vars
if expdef.params.CBS
    unigrids={1:exogenv.exnpt,expdef.KBpts,expdef.LBpts,expdef.BGpts};
else
    unigrids={1:exogenv.exnpt,expdef.KBpts,expdef.LBpts,expdef.WIpts,expdef.BGpts};
end
basegrid=grid.TensorGrid(unigrids);
            
% build basis for approximating function
if strcmp(approxMode,'spline')
    % spline approximation
    knotvec=[2,2,2,2];
    fbnds=ones(2,4);
end

if isequal(guess_mode,'no_guess')
    inigrid=basegrid;
    solmatguess=repmat(solguessvec',[basegrid.Npt,1]);
    Vmatguess=repmat(Vguessvec',[basegrid.Npt,1]);
%    Vmatguess(:,6)=inigrid.Pointmat(:,4);
%     % separate function for each future exog. state
%     transguess=[repmat(basegrid.Pointmat(:,3),[1,exogenv.exnpt]),...
%         repmat(basegrid.Pointmat(:,4),[1,exogenv.exnpt])];
% initialize transition matrix with all SS points (one row of
% basegrid.Pointmat is a point in the SS), for every possible future
% exogenous states (exogenv.exnpt of them)
    transguess=kron(basegrid.Pointmat(:,2:end), ones(1,exogenv.exnpt));
else
    %inigrid=grid.TensorGrid(basegrid.StateBounds,[exnpt,20,20,15]);
    inigrid=basegrid;
    guessobj=load(guess_path,'mobj');
	solmatguess=guessobj.mobj.Pfct.evaluateAt(inigrid.Pointmat)';
    Vmatguess=guessobj.mobj.Vfct.evaluateAt(inigrid.Pointmat)';
    transguess=guessobj.mobj.Tfct.evaluateAt(inigrid.Pointmat)';
	nsol = length(solguessvec);
	if guessobj.mobj.NSOL < nsol
		% Guess has fewer variables to guess
		% Check if new solution variables used to be an add var
		[newvars,idx] = setdiff( fieldnames(stv{1}.Sol), guessobj.mobj.Sol_names,'stable');
		npt=size(inigrid.Pointmat,1);
		newvarmat = zeros(npt,numel(newvars));
		guessobj.mobj=guessobj.mobj.augmentParams(allexpers.bench.params);
		for ii=1:numel(newvars)
			newvar = newvars{ii};
			if ismember(newvar,guessobj.mobj.Add_names)
				% Yes, it has. Use the add var value.
				for jj=1:npt
					point = inigrid.Pointmat(jj,:);
					solvec = solmatguess(jj,:);
					[~,outstr]=ProductionModel_july.calcStateTransition(guessobj.mobj,inigrid.Pointmat(jj,:), ...
						solmatguess(jj,:),0,transguess(jj,:)',Vmatguess(jj,6));
					newvarmat(jj,ii) = log(outstr.addvars.(newvar));
				end
			else
				% No, it hasn't. Use steady state.
				newvarmat(:,ii) = stv{1}.Sol.(newvar) * ones(npt,1);
			end
		end
		solmatguess_new = zeros(npt,nsol);
		solmatguess_new(:,idx) = newvarmat;
		solmatguess_new(:,setdiff(1:nsol,idx)) = solmatguess;
		solmatguess = solmatguess_new;
	end
    
end                  
            
% build approximating function
if strcmp(approxMode,'spline')
    % create guess for solution variable functions
    Pf=grid.SplineFunction(inigrid,solmatguess,knotvec,fbnds);
    % create guess for next-period functions
    Vf=grid.SplineFunction(inigrid,Vmatguess,knotvec,fbnds);
    % and for state transition function
    Tf=grid.SplineFunction(inigrid,transguess,knotvec,fbnds);
    % and for temporary value function (for risk-taker
    % bankruptcy)
    Tempf=grid.SplineFunction(inigrid,zeros(inigrid.Npt,1),knotvec,fbnds);
else
    % build approximating function
    % create guess for solution variable functions
    Pf=grid.LinearInterpFunction(inigrid,solmatguess);
    % create guess for next-period functions
    Vf=grid.LinearInterpFunction(inigrid,Vmatguess);
    % and for state transition function
    Tf=grid.LinearInterpFunction(inigrid,transguess);
    % and for temporary value function (for risk-taker
    % bankruptcy)
    Tempf=grid.LinearInterpFunction(inigrid,zeros(inigrid.Npt,1));
end

if useJacobian
   J=jacfun(expdef.params); 
else
   J=[];
end
expdef.params.useJacobian = useJacobian;
            
% create new object
disp('Creating object...');
mobj=instfun(expdef.params,endogenv,exogenv,Vf,Pf,Tf,Tempf,basegrid,J);
disp('Created object.');

disp('-------------------------------------');
disp('Bounds of state space:');
disp(num2str(mobj.Vfct.SSGrid.StateBounds));

% ===========================================
% save initial object
% ===========================================

disp('Saving object...');
save(outfname,'stv','mobj');              
disp('Saved object.');





