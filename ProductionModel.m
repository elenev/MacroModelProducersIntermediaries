 classdef ProductionModel < model.DSGEModel
    
    properties (SetAccess=protected)
        NSTEN % number of endogenous state variables: Vfct.Ndim-1 
        NSTEX % number of exogenous state variables
        NSOL % number of solution vars
        NV % number of forecasting variables
        NADD % number of additional endogenous variables
        NCOND % number of conditional expectations of period-ahead variables  
		NZNS % number of zero net supply prices
        Sol_names % NSOLx1 cell array with names of solution variables
                   % must be in order of functions of Pfct 
        V_names % NVx1 cell array with names of forecasting variables
                   % must be in order of functions of Vfct      
        Sol_baseguess           
        En_names % NSTENx1 cell array with names of endog state variables           
        Ex_names % NSTEXx1 cell array with names of exog state variables           
        Add_names % NADDx1 cell array with names of additional vars
        Cond_names %NCONDx1 cell array with names of conditional variables  
		Zns_names % NZNSx1 cell array with names of zero net supply prices
        Params % structure array with parameters
        Exogenv % structure array with fields mtrans, pts_perm, pts_all
                % specifying exog. Markov processes
        Vfct % ApproxFunction object for iteration-relevant functions
        Pfct % ApproxFunction object for solution jump variables
        Tfct % ApproxFunction object for transition of state variable(s) (optional)
		Zfct % ApproxFunction object for prices of assets in zero-net supply
        Tempfct % ApproxFunction for intermediate computations 
        Basegrid % base grid 
        Jacobian % Function evaluating the nonlinear system Jacobian
    end    
    
    
    methods
        % constructor
        function obj=ProductionModel(params,endogenv,exogenv,vfct,pfct,tfct,tempfct,basegrid,jac,varargin)
            % call superclass constructor
            obj=obj@model.DSGEModel(params,endogenv,exogenv,vfct,pfct,tfct);
            obj.Tempfct=tempfct;
            obj.Basegrid=basegrid; 
            obj.NCOND=length(endogenv.condnames);          
            obj.Cond_names=reshape(endogenv.condnames,1,obj.NCOND);
            obj.Jacobian=jac;
			if nargin==10
				obj.Zfct=varargin{1};
				obj.Zns_names=endogenv.znsnames;
				obj.NZNS=numel(obj.Zns_names);
			else
				obj.Zns_names=cell(0);
				obj.NZNS=0;
			end
        end
        
        function [nextst,outstr]=calcStateTransition(obj,point,solvec,mode,varargin)
            params=obj.Params;
            % unpack relevant params
            gammaB=params.gammaB;
            gammaS=params.gammaS;
            sigmaI=params.sigmaI;
            kappa=params.kappa;
            mu_G=params.mu_G;
            mu_ZA=params.mu_ZA;
            tau=params.tau;
            tauK=params.tauK;
            tauPi=params.tauPi;
            tauD=params.tauD;
            alpha=params.alpha;
            deltaK=params.deltaK;
            Lbar=params.Lbar;
            delta=params.delta;
            psi=params.psi;
            eta=params.eta;
            theta=params.theta;
            mu_om=params.mu_om;
            T=params.T;   
            bT=params.bT;
            bTau=params.bTau;
            bankDWL=params.bankDWL;
            bGamma=params.bGamma;
            Pbar=params.Pbar;
            gammaG=params.gammaG;
            Lscale=params.Lscale;
            BG_foreign=params.BG_foreign;
            tauK_int=(1-theta)*tauK;
            tauPi_int=(1-theta)*tauPi;
            outputDWL=params.outputDWL;
            sigma_eps=params.sigma_eps;
            rho=params.rho;
            pibar=params.pibar;
			shieldI=params.shieldI;
			phi0=params.phi0;
			phi1=params.phi1;
            n0=params.n0;
            chiLinear=params.chiLinear;
            phiI=params.phiI;
            eqfct_intern=params.eqfct_intern;
            fixediv=params.fixediv;
			chi1adj=params.chi1adj;
            FLnodef=params.FLnodef;
            sigIexp=params.sigIexp;
            etaI=params.etaI;
            zetaI=params.zetaI;
            
            % extract state variables
            exst=point(1);
            KB=point(2);
            LB=point(3);
            WI=point(4);
            BG=point(5);
            ZA=obj.Exogenv.pts_perm(exst,1);
            sig2B=obj.Exogenv.pts_perm(exst,2);
            zeta=obj.Exogenv.pts_perm(exst,3);
            xi=obj.Exogenv.pts_perm(exst,4);
			AStarget=obj.Exogenv.pts_perm(exst,5);
            
            fullendogvars=[KB,LB,WI,BG];

            
            % compute next period's states
			MITshock_WS = 0;
            if isempty(varargin)
                State_next=obj.evaluateTrans(point);
                thisval=obj.evaluateVal(point);
                VI=thisval(6);
            else
                State_next=varargin{1};
                VI=varargin{2};
				if nargin>6
					MITshock = varargin{3};
					MITshock_WS = MITshock;
				end
            end
            F_eps=fastnormcdf(VI+rho,0,sigma_eps);
            F_eps_cond=normpdf((VI+rho)/sigma_eps);
            F_eps_minus=-sigma_eps*F_eps_cond;
            F_eps_plus=sigma_eps*F_eps_cond;
            sigmaI=sigmaI*F_eps^sigIexp;
                                    
            % extract solution variables
            qB=exp(solvec(1));
            X=solvec(2);
            q=exp(solvec(3));
%            cB=exp(solvec(4));
            equ=solvec(4);
            if sigmaI==0
                eI=solvec(5);
                eIcost_intern = 0;
                eIcost_extern = 0;
            else
                cI=exp(solvec(5));
                if eqfct_intern
                    eI = (1-cI)/sigmaI;
                    eIcost_intern = sigmaI/2 * eI^2;
                    eIcost_extern = 0;
                else
                    eI = (cI-1)/sigmaI;                
                    eIcost_intern = 0;
                    eIcost_extern = sigmaI/2 * eI^2;
                end
            end
            cS=exp(solvec(6));  
            wvec=exp(solvec(7:8)); 
            wvec=wvec(:);
			AI_g = solvec(9);
			AS_g = solvec(10);
			cB = exp(solvec(11));
            
            % capital price           
            p=1 + psi*(X/KB-(mu_G-1+deltaK));
            % production and income
            Lbar_scale=Lscale*Lbar;
            L=Lbar_scale(1)^gammaB * Lbar_scale(2)^gammaS;
            Y=mu_om*ZA*KB^(1-alpha) * L^alpha;
            Y_om=Y/mu_om;
            wagebill=Lbar_scale'*wvec;
            
            % investment and capital
            Psi=psi/2*(X/KB-(mu_G-1+deltaK))^2;
            KB_g=(1 - deltaK + X/KB)*KB;            
            
            % borrower payoff risk
			AB = LB*p*KB/qB;
            %AB=(p*KB-WB)/qB;
            if FLnodef
                omstar=0; % turn off liquidity default
                OmA=1;
                OmK=1;
                fom=0;
            else
                omstar=(AB + wagebill + pibar*KB)/Y_om;
                [OmA,OmK,fom]=obj.Mpayoff_gamma(omstar,mu_om,sig2B);                
            end
            omplus=OmK;
            omminus=(1-OmK);
     
            
            % producer budget
            Gfun = @(om1,om2,tau,zetaDWL)( (1-zetaDWL*params.outputDWL)*(1-tau)*om1*Y_om - om2*(1-tau)*wagebill ...
                            + om2*( (1-zetaDWL)*(1-(1-tau)*deltaK)*p - (1-tau)*pibar)*KB );                            
            N = Gfun(omplus,OmA,tauPi,0) - OmA*(1-tauPi*(1-theta)+delta*qB)*AB  + (1-OmA)*n0;
            eq_N = equ/N;
            if eqfct_intern
                equcost_intern = phi1/2 * eq_N^2;
                equcost_extern = 0;
            else
                equcost_intern = 0;
                equcost_extern = phi1/2 * eq_N^2;
            end
            AB_g = (p*KB_g - (1-phi0)*N - equ + equcost_intern)/qB;
            Nminus = Gfun(omminus,1-OmA,0,zeta);
            DivP = N*(phi0 - eq_N - equcost_extern) - (1-OmA)*n0;
            nstar = (Gfun(omstar,1,tauPi,0) - (1-tauPi*(1-theta)+delta*qB)*AB )/N;
            
            
            % payoff to intermediary
            MP= OmA*(1-tauPi_int) + Nminus/AB;
            YB=OmK*ZA*KB^(1-alpha)*L^alpha;                       
                        
            % government debt and taxes
            tauN=tau(BG/Y);
            DebtB=(MP+delta*OmA*qB)*AB;
            WS=DebtB+BG-WI-BG_foreign*BG+MITshock_WS; 
%            ZA_factor=ZA/mu_ZA;
            ZA_factor=min(1,ZA/mu_ZA);
            tau_eff=tauN * ZA_factor^bTau; % procyclical 
            Tvec=Pbar * T * Y * ZA_factor^bT; % countercyclical 
            public_goods = gammaG * Y * ZA_factor^bGamma; % countercyclical 

            % intermediary default
            if fixediv
                phiIWI=phiI*params.WIbar;
            else
                phiIWI=phiI*WI;
            end
            dI_eff = phiIWI - eI - eIcost_extern - F_eps_minus - (1-F_eps)*WI;
            bailout = F_eps_plus - (1-F_eps)*(WI - bankDWL*zetaI*DebtB);            

			% Rebate
			if chi1adj>0
				PsiSA = chi1adj/2 * (AS_g / AStarget - 1)^2 * AStarget;
			else
				PsiSA = 0;
			end
			PsiS = PsiSA + chiLinear * AS_g;
			PsiP = N*phi1/2*eq_N^2;
            BR_expense = Pbar * ((1-eta)*(zeta*(1-OmA)*(1-deltaK)*p*KB + zeta*outputDWL*(1-OmK)*ZA*KB^(1-alpha)*L^alpha) ...
                                   + bankDWL*(1-etaI)*(1-F_eps)*zetaI*DebtB + PsiS + PsiP + pibar*KB);            
            
            % budget constraint borrower
            %cB=(1-tau_eff)*wvec(1)*Lbar_scale(1) + Tvec(1) + BR_expense(1) + p*X + dI_eff + DivP - X - Psi*KB;
            netResourcesB = (1-tau_eff)*wvec(1)*Lbar_scale(1) + Tvec(1) + BR_expense(1) + p*X + dI_eff + DivP - X - Psi*KB;
			
            % intermediary
            rD=1/q-1;
            bI=(WI - phiIWI + eI - eIcost_intern - qB*AI_g)/(q+shieldI*tauPi*rD-kappa);            
            % saver
            bS=( (1-tau_eff)*wvec(2)*Lbar_scale(2) + WS + Tvec(2) + BR_expense(2) - cS - qB*AS_g - PsiS)/(q+tauD*rD);          
           
            % next period government debt
             
            spending = public_goods + sum(Tvec) + bailout;
            taxes = (tau_eff-tauK*OmA)*wagebill + tauK*YB + OmA*AB*(tauPi_int-tauK_int) + bI*(shieldI*tauPi*rD-kappa) + bS*tauD*rD - tauK*OmA*(deltaK*p + pibar)*KB;
            BG_g = (BG + spending - taxes)/q;            
            
            if mode>0
                % simulation, mode contains number of next period's state
                exst=mode;
                exnpt=size(obj.Exogenv.pts_all,1);
				
				% force in state space bounds
				for ii=1:4
					unigrid = obj.Tfct.SSGrid.Unigrids{1+ii}; 
					minbnd = min( unigrid );
					maxbnd = max( unigrid );
					idx = exnpt*(ii-1)+(1:exnpt);
					State_next( idx ) = max( State_next( idx ), minbnd );
					State_next( idx ) = min( State_next( idx ), maxbnd );
				end
				
                ZAnext=obj.Exogenv.pts_all(exst,1);
                cind=obj.Exogenv.pts_all(exst,end);     
                sig2B_zeta_xi_next=obj.Exogenv.pts_all(exst,2); %:end-1
                KBnext=State_next(exst,1);
                LBnext=State_next(exnpt+exst,1);
                WInext=State_next(2*exnpt+exst,1);
                BGnext=State_next(3*exnpt+exst,1);
                nextst=[cind,KBnext,LBnext,WInext,BGnext]; 
                if ~isreal(nextst)
                    disp('dang');
                end
            else
                % solution mode, compute next period's state variables for all
                % possible Markov states
                ZAnext=obj.Exogenv.pts_all(:,1);
                %exnpt=size(obj.Exogenv.pts_all,1);
                sig2B_zeta_xi_next=obj.Exogenv.pts_all(:,2:end-1);
                cind=obj.Exogenv.pts_all(:,end);   
                nst=length(ZAnext);
                KBnext=State_next(1:nst,1);
                LBnext=State_next(nst+(1:nst),1);
                WInext=State_next(2*nst+(1:nst),1);                
                BGnext=State_next(3*nst+(1:nst),1);
                % matrix of next period states
                nextst=[cind,KBnext,LBnext,WInext,BGnext];
			end            
                               
			addvars=struct;
			addvars.OmA=OmA;
			addvars.OmK=OmK;
			addvars.MP=MP;
			addvars.omstar=omstar;
			addvars.fom=fom;
			addvars.omplus=omplus;
			addvars.F_eps=F_eps;
			addvars.F_eps_minus=F_eps_minus;
			addvars.F_eps_plus=F_eps_plus;
			addvars.bI=bI;
			addvars.VI=VI;
			addvars.bS=bS;
			addvars.BG_g=BG_g;
			addvars.AB=AB;
			addvars.KB=KB;
            addvars.WI=WI;
			addvars.L=L;
			addvars.AB_g=AB_g;
			addvars.KB_g=KB_g;
			addvars.sig2B=sig2B;
			addvars.sig2B_zeta_xi_next=sig2B_zeta_xi_next;
			addvars.Y_om=Y_om;                        
			addvars.YB=YB;                        
			addvars.zeta=zeta;
			addvars.xi=xi;
			addvars.AStarget=AStarget;
			addvars.tauN=tauN;
			addvars.ZA=ZA;
			addvars.ZA_factor=ZA_factor;
            addvars.cB=cB;
            addvars.N=N;
            addvars.nstar=nstar;
            addvars.DivP=DivP;
			addvars.netResourcesB=netResourcesB;
                     
            outstr=struct;
            outstr.addvars=addvars;
            outstr.exstvec=[ZA;sig2B;zeta;xi;AStarget];
            outstr.Jvec=[omstar,fom,mu_om,OmK,BG/Y,F_eps,F_eps_cond];
            outstr.fullendogvars=fullendogvars; 
        end
        
        
        function [fx,J,V]=calcEquations(obj,exst,nextst,solvec,instr,mode,varargin)
            % allocate result
            fx=zeros(obj.NSOL,1);
            
            % unpack params
            params=obj.Params;
            betaB=params.betaB;
            betaS=params.betaS;
            mu_G=params.mu_G;
            kappa=params.kappa;
            theta=params.theta;
            tauK=params.tauK;
            tauPi=params.tauPi;
            tauD=params.tauD;
            alpha=params.alpha;
            deltaK=params.deltaK;
            delta=params.delta;
            F=params.F;
            Phi=params.Phi;
            psi=params.psi;
            pibar=params.pibar;
            mu_om=params.mu_om;
            try
               nuinv=1./[params.nuB,params.nuS];
            catch
               nuinv=1/params.nu*ones(2,1);
            end
            tauK_int=(1-theta)*tauK;
            BG_foreign=params.BG_foreign;
            gammaB=params.gammaB;
            gammaS=params.gammaS;
            Lbar_scale=params.Lbar*params.Lscale;
            sigmaI=params.sigmaI;
            rho=params.rho;
            sigma_eps=params.sigma_eps;
            mktConsI=params.mktConsI;
			shieldI=params.shieldI;
			phi1=params.phi1;
			chiLinear=params.chiLinear;
            phiI=params.phiI;
            eqfct_intern=params.eqfct_intern;
            fixediv=params.fixediv;
			chi1adj=params.chi1adj;
            sigIexp=params.sigIexp;
                        
            % extract some other state-dependent variables
            envec=instr.addvars;
            KB=envec.KB;
            AB=envec.AB;
            WI=envec.WI;
            OmK=envec.OmK;
            OmA=envec.OmA;
            fom=envec.fom;
            omstar=envec.omstar;
            AB_g=envec.AB_g;    
            KB_g=envec.KB_g;
            BG_g=envec.BG_g;
            bI=envec.bI;
            VI=envec.VI;
            bS=envec.bS;
            L=envec.L;
            ZA=envec.ZA;
			ZA_factor=envec.ZA_factor;
            Y_om=envec.Y_om;
            xi=envec.xi;
			AStarget=envec.AStarget;
            N=envec.N;
            nstar=envec.nstar;
			netResourcesB=envec.netResourcesB;
            F_eps=envec.F_eps;
            sigmaI=sigmaI*F_eps^sigIexp;
            
            
            % extract endogeous variables 
            qB=exp(solvec(1));
            X=solvec(2);
            q=exp(solvec(3));
            equ=solvec(4);
            if sigmaI==0
                cI=solvec(5);
                eI=cI;
                cI_norm=1;
                eIcost_extern = 0;
            else
                cI=exp(solvec(5));
                if eqfct_intern
                    eI = (1-cI)/sigmaI;
                    cI_norm = 1/cI;
                    eIcost_extern = 0;
                else
                    eI = (cI-1)/sigmaI;                
                    cI_norm = cI;
                    eIcost_extern = sigmaI/2 * eI^2;
                end
            end
            cS=exp(solvec(6));  
            wvec=exp(solvec(7:8));
			AI_g = solvec(9);
			AS_g = solvec(10);
			cB = exp(solvec(11));
            muR=solvec(12);
            lamI=solvec(13);
            lamS=solvec(14);
			muS=solvec(15);
            lamB=solvec(16);
            wvec=wvec(:);                

            % multiplier transformations
            muRplus=max(0,muR)^3;
            muRminus=max(0,-muR)^3;           
            lamIplus=max(0,lamI)^3;
            lamIminus=max(0,-lamI)^3;           
            lamSplus=max(0,lamS)^3;
            lamSminus=max(0,-lamS)^3;           
			muSplus=max(0,muS)^3;
            muSminus=max(0,-muS)^3; 
            lamBplus=max(0,lamB)^3;
            lamBminus=max(0,-lamB)^3;           
            
            
            % capital price
            p = 1 + psi*(X/KB-(mu_G-1+deltaK));    
			
			% saver adjustment cost
			if chi1adj > 0
				dPsiSA = chi1adj*(AS_g/AStarget-1);
			else
				dPsiSA = 0;
			end
            
            % probabilities and states to compute expectation terms
            prnext=obj.Exogenv.mtrans(exst,:);          
            ZAnext=obj.Exogenv.pts_perm(:,1);            
            
            % projection evaluation
            if nargin==7
                % Exactly 7 arguments were passed.
                % Means expectations have already been computed (fast
                % solution mode)
                expStruct = varargin{1};
            elseif nargin==6
                % Nothing passed in ((conditional moments and EE errors in
                % simulation)
                expStruct = obj.computeExpectations(exst,nextst);
            else
                % trans, vtrans, and vtrans_def were passed in (legacy)  
                expStruct = obj.computeExpectations(exst,nextst,varargin{:});
            end            

            exp_B_qB = expStruct.exp_B_qB;
            exp_B_qK = expStruct.exp_B_qK;
            exp_I_qB = expStruct.exp_I_qB;
            exp_I_q = expStruct.exp_I_q;
            exp_S_q = expStruct.exp_S_q;
			exp_S_qB = expStruct.exp_S_qB;
            CE_vec = expStruct.CE_vec;
            PayoffI = expStruct.PayoffI;
            p_next = expStruct.p_next;
            qB_next = expStruct.qB_next;
            exp_VI_next = expStruct.exp_VI_next;
             
            % compute certainty equiv. and SDFs
            U_vec=zeros(2,1);
            V_vec=zeros(2,1);
            U_norm=zeros(2,1); % to normalize SDFs back
            cvec=[cB,cS];
            beta_vec=[betaB,betaS];
            U_vec(1)=cB;
            U_vec(2)=cS;
            nuinv_vec=nuinv;
            for i=1:2
                if nuinv_vec(i)==1
                   V_vec(i)=U_vec(i)^(1-beta_vec(i)) *  CE_vec(i)^beta_vec(i); 
                else
                   V_vec(i)=( (1-beta_vec(i))*U_vec(i)^(1-nuinv_vec(i)) + beta_vec(i)* CE_vec(i)^(1-nuinv_vec(i)) )^(1/(1-nuinv_vec(i))); 
                end
                U=U_vec(i)^(1-nuinv_vec(i))/cvec(i);  
                U_norm(i)=U;
            end      
            
            % intermediary
            if eqfct_intern
                equnorm=1/(1-phi1*equ/N);
            else
                equnorm=1+phi1*equ/N;
            end
            U_normI=U_norm(1)*cI_norm;
            U_normP=U_norm(1)*equnorm;
                                   
            % borrower FOCs
            fx(1) = qB - F*lamBplus - exp_B_qB/U_normP;
            fx(2) = p - Phi*p*lamBplus - exp_B_qK/U_normP;
            MPLvec=[gammaB, gammaS]./Lbar_scale' * alpha*mu_om*ZA*L*(KB/L)^(1-alpha);
            Fscript=nstar/Y_om;
            for j=1:2
                fx(2+j)= (1-tauK)*OmK*MPLvec(j)/mu_om - (1-tauK)*OmA*wvec(j) - fom*Fscript*(wvec(j) - omstar*MPLvec(j));
            end                        
            
            % intermediary FOC
            rD=1/q-1;
            fx(5)=qB*(1 - xi*lamIplus) - muRplus - exp_I_qB/U_normI;
            if mktConsI
                fx(6)=q + shieldI*tauPi*rD - q*lamIplus - kappa - exp_I_q/U_normI;
            else
                fx(6)=q + shieldI*tauPi*rD - lamIplus - kappa - exp_I_q/U_normI;
            end
                                      
            % saver FOC
            fx(7)=q + tauD*rD -lamSplus - exp_S_q/U_norm(2);
            
            % risk taker constraints
            if mktConsI
                fx(8)=xi*qB*AI_g + q*bI - lamIminus; 
            else
                fx(8)=xi*qB*AI_g + bI - lamIminus; 
            end
            
            % saver constraint
            fx(9)=bS - lamSminus;
            
            % risk free debt market clearing
            fx(10)=BG_g - bS - bI - BG_g*BG_foreign;
            % intermediary no-shorting constraint
            fx(11)=AI_g - muRminus;
            			
			% corporate bond market
			fx(12) = qB - muSplus + dPsiSA + chiLinear  - exp_S_qB/U_norm(2);
			fx(13) = AS_g - muSminus;
			fx(14) = AB_g - AS_g - AI_g;
			
            % budget constraint for borrower
			fx(15) = cB - netResourcesB;
            
            % borrower hard leverage constraint
            fx(16)=Phi*p*KB_g - F*AB_g - lamBminus;

			
            % Transitions
%             nst=size(nextst,1);
            nst = obj.Exogenv.exnpt;
            KBtrans = KB_g * ones(nst,1)./ mu_G;
            BGtrans = BG_g * ones(nst,1)./ mu_G;
            WInext = PayoffI*AI_g + bI;
            WItrans = WInext ./ mu_G;
            LBtrans = qB_next*AB_g ./ ( p_next*KB_g ); % PayoffB;
            
            % Jacobian computation (numeric)
            if params.useJacobian && mode == 0
                fullendogvars=instr.fullendogvars;
                ZA=instr.exstvec(1);
                sig2B=instr.exstvec(2);
                zeta=instr.exstvec(3);
                xi=instr.exstvec(4);
                
                % Compute additional things needed for Jacobian:
                % d OmA / d omstar, d OmK / d omstar
                omstar = instr.Jvec(1);
                fom = instr.Jvec(2);
                mu_om = instr.Jvec(3);
                OmK = instr.Jvec(4); 
                F_eps = instr.Jvec(6);
                F_eps_cond = instr.Jvec(7);
                dOmA = -fom;
                dfom = -fom*(omstar*mu_om - mu_om^2 + sig2B)/(omstar*sig2B); 
                dOmK = -omstar*fom;
                dF_eps = -normpdf(VI+rho,0,sigma_eps);
                dF_eps_cond = dF_eps*(VI+rho);

                % tau - tau0 and dTau to handle symbolic piecewise tax rule for
                % version earlier than Matlab R2016b
                %tau0 = params.tau0;
                tauN = envec.tauN;% - tau0; 
                dTau = obj.Params.dTau(instr.Jvec(5)); 
                if range(tauN)>0
                    error('Analytical derivative cannot be computed for this tax rule');
                end
                
                % Compute Jacobian
                                
                args = [{ZA,sig2B,zeta,xi,AStarget,ZA_factor},num2cell(fullendogvars),{qB,X,q,equ,...
					cI,cS,wvec(1),wvec(2),AI_g,AS_g,cB...
                muRplus,muRminus,lamIplus,lamIminus,lamSplus,lamSminus,muSplus,muSminus,lamBplus,lamBminus,...
                OmA,fom,OmK,tauN,dOmA,dTau,dfom,dOmK,...
                F_eps, F_eps_cond, dF_eps, dF_eps_cond,...
                exp_B_qB,exp_B_qK,exp_I_qB,exp_I_q,VI,exp_S_q,exp_S_qB}];

                % Jacobian with respect to actual solution variables
				% jacfun.m must exist in the folder for this line to work
				% otherwise, use obj.Jacobian
                tmp_Jacobian = jacfun(args{:});
				%tmp_Jacobian = real(tmp_Jacobian);

                nsol = size(solvec,1);
				nfx = nsol;
				%nfx = nsol + 1;
				
                J = zeros(nfx,nsol);

                % Now use chain rule to compute Jacobian with respect to
                % solvec
				
                
                % Exponentiated variables
                if sigmaI~=0
                    expind=[1,3,5:8,11];
                    nexpind=[2,4,9,10];
                else
                    expind=[1,3,6:8,11];
                    nexpind=[2,4,5,9,10];
                end
                
                J(:,expind) = tmp_Jacobian(1:end,expind) .* repmat( exp(solvec(expind))' ,nfx,1);
                                    
                % Investment and portfolios (not exponentiated)
                J(:,nexpind) = tmp_Jacobian(1:end,nexpind);                

                % Multipliers
                start_idx = 11;
                
                for mult_idx = 1:nsol-start_idx
                    mult = solvec(start_idx + mult_idx);
                    if mult > 0
                        J(:,start_idx + mult_idx) = ...
                            tmp_Jacobian(1:end,start_idx + 2*(mult_idx-1) + 1) .* repmat( 3*mult^2 ,nfx,1);
                    else
                        J(:,start_idx + mult_idx) = ...
                            tmp_Jacobian(1:end,start_idx + 2*(mult_idx-1) + 2) .* repmat( -3*mult^2 ,nfx,1);
                    end
                end
            else
                J=[];
            end
            
            V=cell(3,1);
            % marginal value functions
            if mode==1
                % Output new values for time iteration
                % risk-taker value function evaluated at zero risk-taker
                Vnext=zeros(obj.Vfct.Nof,1);
                Vnext(1)=cB;
                Vnext(2)=cS;
                if sigmaI==0
                    Vnext(3)=1;
                else
                    Vnext(3)=cI;
                end
                Vnext(4)=V_vec(1);
                Vnext(5)=V_vec(2);
                if fixediv
                    phiIWI=phiI*params.WIbar;
                else
                    phiIWI=phiI*WI;
                end
                Vnext(6)=phiIWI - eI - eIcost_extern + exp_VI_next/U_norm(1);
                Vnext(7)=X;
                Vnext(8)=qB;
                Vnext(9)=Lbar_scale'*wvec;
                Vnext(10)=lamB;
                V{1}=Vnext;
                % state transition
                V{2}=[KBtrans; LBtrans; WItrans; BGtrans]';
                
                % check resource constraint
                % public_goods = params.gammaG * instr.addvars.Y * exp(params.g) / obj.Exogenv.pts_perm(exst,1);
                % RC = cI + cB + cS + X + psi/2*(X/KB-(mu_G-mu_om)).^2*KB + public_goods + zeta*p*KB*(mu_om-ZK) - instr.addvars.Y;
                % if abs(RC)>0.00001
                %   disp('RC violation'); 
                % end
        
            elseif mode==2
                
                % Evaluation during simulation. Output conditional
                % variables
                
                % SDFs
                SDFI = expStruct.Conditional.SDFI/U_normI; % normalize SDFs back with current utility
                SDFB = expStruct.Conditional.SDFB/U_norm(1);
                SDFS = expStruct.Conditional.SDFS/U_norm(2);
                CB_next = expStruct.Conditional.CB_next;
                CS_next = expStruct.Conditional.CS_next;
                CI_next = expStruct.Conditional.CI_next;
                OmA_next = expStruct.Conditional.OmA_next;
                OmK_next = expStruct.Conditional.OmK_next;
                MP_next = expStruct.Conditional.MP_next;
                qB_next = expStruct.Conditional.qB_next;                
                fom_next = expStruct.Conditional.fom_next;
                omstar_next = expStruct.Conditional.omstar_next;
                muom_next = expStruct.Conditional.muom_next;
                muom_next_norm = reshape(muom_next,size(OmA_next));
                Fscript_next = expStruct.Conditional.Fscript_next;
                VI_next_ret = expStruct.Conditional.VI_next_ret;
               
                SDF.SDFI = SDFI;
                SDF.SDFI_noBR = SDFI(:,1);
                SDF.SDFB = SDFB;    
                SDF.SDFS= SDFS;
                
                % Define returns to intermediaries on corporate bonds
                retP = ( MP_next+tauPi*(1-theta)*OmA_next + delta*OmA_next.*qB_next ) / qB;
                expRP = prnext * retP;
                expOmA = prnext * OmA_next;
                stdRP = sqrt( prnext * ( retP - expRP).^2 );
                CovP = prnext* ( (SDF.SDFI - prnext*SDF.SDFI) .* (retP - expRP) );
                
                % Decomposition of risk premia
                % Note: if lamI = 0, then effect_Rf = 1 and effect_collatP = 0
                % and the standard asset pricing equation R-Rf =
                % -Rf*Cov[SDF,R] is recovered
                Rf=1./q;
                effect_Rf = 1/(1 - lamIplus*Rf);
                effect_CovP = -CovP;
                effect_collatP = lamIplus*(Rf-xi);                
                expRP_check = Rf*effect_Rf*(effect_CovP + effect_collatP);        
                
                % Return on bank equity
                VI_ex_div = exp_VI_next/U_norm(1);
                retBE = prnext*VI_next_ret/VI_ex_div;
                                
                % Define returns to borrowers on corporate bonds
                retPB = OmA_next.*(1-tauK_int+delta*qB_next)/qB;
                expRPB = prnext * retPB;
                stdRPB = sqrt( prnext * ( retPB - expRPB).^2 );
                CovPB = prnext* ( (SDF.SDFB - prnext*SDF.SDFB) .* (retPB - expRPB) );
                
                % Decomposition of risk premia
                lamBS = q - prnext*SDF.SDFB;
                effect_RfB = 1/(1 - lamBS*Rf);
                effect_CovPB = -CovPB;
                effect_collatPB = 0;                
                expRPB_check = Rf*effect_RfB*(effect_CovPB + effect_collatPB);                  
                          
                % Define returns to borrowers
                MPK_next=MP_next;               
                retB = (OmA_next.*(p_next*(1-(1-tauK)*deltaK) - pibar) + (1-tauK)*OmK_next.*MPK_next./muom_next_norm - fom_next.*Fscript_next.*(pibar-omstar_next.*MPK_next./muom_next_norm) )/p;
                expRB = prnext * retB;
                stdRB = sqrt( prnext * ( retB - expRB).^2 );
                
                % Define expected consumption growth
                CBgr =  (CB_next ./ cB) .* mu_G;
                CIgr =  (CI_next ./ cI) .* mu_G;
                CSgr =  (CS_next ./ cS) .* mu_G;
                expCBgr=prnext * CBgr;
                expCIgr=prnext * CIgr;
                expCSgr=prnext * CSgr;
                
                stdCBgr=sqrt( prnext * (CBgr - expCBgr).^2 );
                stdCIgr=sqrt( prnext * (CIgr - expCIgr).^2 );
                stdCSgr=sqrt( prnext * (CSgr - expCSgr).^2 );
                
                % Expected excess return and Sharpe ratio, including and
                % excluding bankrupcty possibility
                
                % intermediaries              
                start_wealthI=WInext;
                end_wealthI=qB*AI_g+(q+shieldI*tauPi*(1/q-1)-kappa)*bI;
                ROWI=start_wealthI./end_wealthI;
                                
                expROWI=prnext*(ROWI);
                expEROWI=expROWI-(1/q);
                stdROWI=sqrt(prnext * ( ROWI - expROWI ).^2);
                SRI=expEROWI/stdROWI;
                    
                % borrowers
                start_wealthB=p_next*KB_g.*OmA_next.*(1-(1-tauK)*deltaK) + (1-tauK).*OmK_next.*ZAnext.*KB_g^(1-alpha)*L^alpha - OmA_next.*(1-tauK_int+qB_next*delta)*AB_g;
                end_wealthB=(p*KB_g-qB*AB_g);               

                ROWB=start_wealthB./end_wealthB;
                expROWB=prnext*ROWB;
                expEROWB=expROWB-(1/q);
                stdROWB=sqrt(prnext*((ROWB-expROWB).^2));
                SRB=expEROWB/stdROWB;
                
                condvars = struct('expRP',expRP, ...
                              'stdRP',stdRP, ...
                              'effect_Rf',effect_Rf, ...
                              'effect_CovP',effect_CovP, ...
                              'effect_collatP',effect_collatP, ...
                              'expRP_check',expRP_check, ...
                              'expRPB',expRPB, ...
                              'stdRPB',stdRPB, ...
                              'effect_RfB',effect_RfB, ...
                              'effect_CovPB',effect_CovPB, ...
                              'effect_collatPB',effect_collatPB, ...
                              'expRPB_check',expRPB_check, ...
                              'expRB',expRB, ...
                              'stdRB',stdRB, ...
                              'expOmA',expOmA, ...
                              'expCBgr',expCBgr, ...
                              'expCIgr',expCIgr, ...
                              'expCSgr',expCSgr, ...
                              'stdCBgr',stdCBgr, ...
                              'stdCIgr',stdCIgr, ...
                              'stdCSgr',stdCSgr, ... % EER and SR
                              'expEROWB',expEROWB, ...
                              'SRB',SRB, ...
                              'expEROWI',expEROWI, ...
                              'SRI',SRI,...
                              'retBE',retBE);                        
                
                Wtrans.KB = KBtrans;                          
                Wtrans.LB = LBtrans;
                Wtrans.WI = WItrans;
                Wtrans.BG = BGtrans;
                V = {condvars,Wtrans,SDF};
            end
        end
    
        
        function expStruct = computeExpectations( obj, exst, nextst, varargin )
            % This function computes E[ ] terms in equilibrium conditions
            
            % unpack params
            params=obj.Params;
            betaB=params.betaB;
            betaS=params.betaS;
            sigmaB=params.sigmaB;
            sigmaS=params.sigmaS;
%             xi=params.xi;
            mu_G=params.mu_G;
            kappa=params.kappa;
            theta=params.theta;
            tauK=params.tauK;
            tauPi=params.tauPi;
            tauD=params.tauD;
            alpha=params.alpha;
            deltaK=params.deltaK;
            delta=params.delta;
            psi=params.psi;
            mu_om=params.mu_om;
            Phi=params.Phi;
            try
               nuinv=1./[params.nuB,params.nuS];
            catch
               nuinv=1/params.nu*ones(2,1);
            end
            tauK_int=(1-theta)*tauK;
            tauPi_int=(1-theta)*tauPi;
            gammaB=params.gammaB;
            gammaS=params.gammaS;
            Lbar_scale=params.Lbar*params.Lscale;
            outputDWL=params.outputDWL;
            rho=params.rho;
            sigma_eps=params.sigma_eps;
            pibar=params.pibar;
			phi0=params.phi0;
			phi1=params.phi1;
            n0=params.n0;
            phiI=params.phiI;
            eqfct_intern=params.eqfct_intern;
            fixediv=params.fixediv;
            FLnodef=params.FLnodef;
            

            % extract some other state-dependent variables
            sig2B_next=obj.Exogenv.pts_perm(nextst(:,1),obj.Exogenv.exidx.Om);
            zeta_next=obj.Exogenv.pts_perm(nextst(:,1),obj.Exogenv.exidx.Zeta);              
            xi_next=obj.Exogenv.pts_perm(nextst(:,1),obj.Exogenv.exidx.xi);
             
            % probabilities and states to compute expectation terms
            prnext=obj.Exogenv.mtrans(exst,:);     
            ZAnext=obj.Exogenv.pts_perm(nextst(:,1),1);
            
            % projection evaluation
            if nargin>3
                Pol_next=varargin{1};
            else
                Pol_next=obj.evaluateVal(nextst)';
            end
            
            CB_next=Pol_next(:,1);
            CS_next=Pol_next(:,2);
            CI_next=Pol_next(:,3);
            VB_next=Pol_next(:,4);
            VS_next=Pol_next(:,5);
            VI_next=Pol_next(:,6);
            X_next=Pol_next(:,7);
            qB_next=Pol_next(:,8);
            wagebill_next=Pol_next(:,9);
            lamBplus_next=max([Pol_next(:,10),zeros(size(ZAnext))],[],2).^3;
            
            % compute intermediary default rate
            F_eps_next=fastnormcdf(VI_next+rho,0,sigma_eps); 
            f_eps_next=normpdf(VI_next+rho,0,sigma_eps); 
            F_eps_minus_next = -sigma_eps*normpdf((VI_next+rho)/sigma_eps);
                       
            
            % compute certainty equiv. and SDFs for borrower and saver
            UB_next=CB_next;
            Vnext_mat={VB_next,VS_next};
            Cnext_mat={CB_next,CS_next};
            Unext_mat={UB_next,CS_next};
            sigma_vec=[sigmaB,sigmaS];
            beta_vec=[betaB,betaS];
            nuinv_vec=nuinv;
            SDF_mat=cell(2,1);
            for i=1:2
                % certainty equivalent
                Vnext_tmp=Vnext_mat{i};
                Cnext_tmp=Cnext_mat{i};
                Unext_tmp=Unext_mat{i};
                if sigma_vec(i)==1          
                    % log
                    CE_tmp=exp(prnext * log(mu_G.*Vnext_tmp));
                else                
                    CE_tmp=(prnext * (mu_G.*Vnext_tmp).^(1-sigma_vec(i)))^(1/(1-sigma_vec(i)));
                end
                Unext_tmp=Unext_tmp.^(1-nuinv_vec(i))./Cnext_tmp;
                SDF_tmp=beta_vec(i)* mu_G.^(-sigma_vec(i))...
                    .* Unext_tmp .* (Vnext_tmp/CE_tmp).^(nuinv_vec(i)-sigma_vec(i));
                CE_vec(i)=CE_tmp;
                SDF_mat{i}=SDF_tmp;
            end
                                  
            % borrower payoff risk next period
            muom_next_norm=mu_om*ones(size(sig2B_next));
            sig2B_next=sig2B_next(:);
            muom_next=mu_om*ones(size(sig2B_next));
            KB_next=nextst(:,2);
            LB_next=nextst(:,3);
            p_next=1+psi*(X_next./KB_next - (mu_G-1+deltaK));
			AB_next=LB_next .* p_next .* KB_next ./ qB_next;
            %AB_next=(p_next.*KB_next-WB_next)./qB_next;
            L=Lbar_scale(1)^gammaB  * Lbar_scale(2)^gammaS; 
            Y_om_next=ZAnext.*KB_next.^(1-alpha) * L^alpha; % normalized by mu_om
            if FLnodef   % turn off liquidity default for firms
                omstar_next=zeros(size(sig2B_next)); 
                OmA_next=ones(size(sig2B_next));
                OmK_next=ones(size(sig2B_next));
                fom_next=zeros(size(sig2B_next)); 
            else
                omstar_next=(AB_next+wagebill_next+pibar*KB_next)./Y_om_next;
                [OmA_next,OmK_next,fom_next]=obj.Mpayoff_gamma(omstar_next,muom_next,sig2B_next);                
            end
            omplus_next=OmK_next;
            omminus_next=(1-OmK_next);

            % producer
            Gfun = @(om1,om2,tau,zetaDWL)( (1-zetaDWL*outputDWL).*(1-tau).*om1.*Y_om_next - om2.*(1-tau).*wagebill_next ...
                            + om2.*( (1-zetaDWL).*(1-(1-tau)*deltaK).*p_next - (1-tau).*pibar).*KB_next );                            
            N_next = Gfun(omplus_next,OmA_next,tauPi,0) - OmA_next.*(1-tauPi*(1-theta)+delta*qB_next).*AB_next + (1-OmA_next)*n0;
            eq_next = p_next.*KB_next - AB_next.*qB_next - (1-phi0)*N_next;
            Nminus = Gfun(omminus_next,1-OmA_next,0,zeta_next);
            Nstar_next = Gfun(omstar_next,1,tauPi,0) - (1-tauPi*(1-theta)+delta*qB_next).*AB_next;
            
            % payoff to intermediary
            MP_next= OmA_next*(1-tauPi_int) + Nminus./AB_next;                                    
%             MP_next=OmA_next*(1-tauPi_int)+( (1-deltaK)*(1-OmA_next).*(1-zeta_next).*p_next.*KB_next...
%                                             +(1-zeta_next*outputDWL).*(1-OmK_next).*Y_om_next - (1-OmA_next).*wagebill_next )./AB_next;
            
            % borrower FOCs
            if eqfct_intern
                equnorm=1./(1-phi1*eq_next./N_next);
            else
                equnorm=1+phi1*eq_next./N_next;
            end            
            SDFB=SDF_mat{1};
            SDFP=SDFB.*(phi0 + (1-phi0)*equnorm);
            Fscript_next=Nstar_next./Y_om_next;
            FOCtmp=SDFP.*(OmA_next.*(1-tauK_int+delta*qB_next) + fom_next.*Fscript_next);
            exp_B_qB = prnext*FOCtmp;
            
            MPK_next= (1-alpha).*muom_next_norm.* ZAnext .*(KB_next/L).^(-alpha);
            FOCtmp=SDFP.*(OmA_next.*(p_next.*(1-(1-tauK)*deltaK) - pibar) ...
                            + (1-tauK)*OmK_next.*MPK_next./muom_next_norm ...
                            - fom_next.*Fscript_next.*pibar ...
                            + fom_next.*omstar_next.*MPK_next.*Fscript_next./muom_next_norm );
            exp_B_qK=prnext* FOCtmp;
                         
            % risk taker FOC
            if eqfct_intern
                CI_norm=1./CI_next;
            else
                CI_norm=CI_next;
            end
            if ~fixediv
                CI_norm = phiI + (1-phiI)*CI_norm;
            end
            SDFI=SDFB.*F_eps_next.*CI_norm;
            FOCtmp=MP_next+delta*OmA_next.*qB_next;
            exp_I_qB =  prnext*(SDFI(:,1).*FOCtmp);
            exp_I_q = prnext*SDFI;
            VI_next_ret = F_eps_next.*VI_next-F_eps_minus_next - (1-F_eps_next)*rho;
            exp_VI_next = prnext*(SDFB.*VI_next_ret);
                                      
            % saver FOC
            SDFS=SDF_mat{2};
            exp_S_q = prnext*SDFS;
			FOCtmp=MP_next+delta*OmA_next.*qB_next;
			exp_S_qB = prnext*(SDFS.*(FOCtmp));
            
            PayoffI = MP_next + delta*OmA_next.*qB_next; % in WInext
            
            expStruct = struct;
			expStruct.exp_B_qB=exp_B_qB;
			expStruct.exp_B_qK=exp_B_qK;
			expStruct.exp_I_qB=exp_I_qB;
			expStruct.exp_I_q=exp_I_q;
			expStruct.exp_VI_next=exp_VI_next;
			expStruct.exp_S_q=exp_S_q;
			expStruct.exp_S_qB=exp_S_qB;
			expStruct.CE_vec=CE_vec;
			expStruct.PayoffI=PayoffI;
			expStruct.p_next=p_next;
			expStruct.qB_next=qB_next;
               
            Conditional = struct;
			Conditional.SDFB=SDFB;
			Conditional.SDFS=SDFS;
			Conditional.SDFI=SDFI;
			Conditional.CB_next=CB_next;
			Conditional.CS_next=CS_next;
			Conditional.CI_next=CI_next;
			Conditional.OmA_next=OmA_next;
			Conditional.VS_next=VS_next;
			Conditional.OmK_next=OmK_next;
			Conditional.MP_next=MP_next;
			Conditional.qB_next=qB_next;
			Conditional.fom_next=fom_next;
			Conditional.omstar_next=omstar_next;
			Conditional.Fscript_next=Fscript_next;
			Conditional.muom_next=muom_next;
			Conditional.VI_next_ret=VI_next_ret;          
               
            expStruct.Conditional = Conditional;           
        end           
        
        
        function [errmat,solmat,condmat,Wshtrans,SDFmat]=calcEEError(obj,pointmat)
            % function to compute Euler equation error at points in state
            % space given by pointmat
            nst=size(obj.Exogenv.pts_all,1);
            
            errmat=zeros(size(pointmat,1),obj.Pfct.Nof);
            solmat=zeros(size(errmat));
            condmat=zeros(size(pointmat,1),obj.NCOND);
            SDFmat=zeros(size(pointmat,1),2*nst);
            Wshtrans=zeros(size(pointmat,1),4*nst);
            
            evaluatePol = @(point)obj.evaluatePol(point);
            calcStateTransition = @(point,soltmp)obj.calcStateTransition(point,soltmp,0);
            calcEquations = @(exst,nextst,soltmp,outstr)obj.calcEquations(exst,nextst,soltmp,outstr,2);    
            % Should be parfor. Use for when debugging only
             parfor i=1:size(errmat,1)
%            for i=1:size(errmat,1)
                point=pointmat(i,:);
                if point(4)==0
                    %disp('brupt');
                end
                soltmp=evaluatePol(point);
                % transition
                [nextst,outstr]=calcStateTransition(point,soltmp);
                % equations
                [fx,~,V]=calcEquations(point(1),nextst,soltmp,outstr);                                
                qB=exp(soltmp(1));
                p=exp(soltmp(2));
                q=exp(soltmp(3));
                wvec=exp(soltmp(7:8));
				cB=exp(soltmp(11));
                AB_g=outstr.addvars.AB_g;
                normvec=[qB,p,wvec',qB,q,q,qB*AB_g,qB*AB_g,qB*AB_g, ...
					qB*AB_g,qB,qB*AB_g,qB*AB_g,cB,qB*AB_g];
                condvars=V{1};
				KBtrans=V{2}.KB;
                LBtrans=V{2}.LB;                
                WItrans=V{2}.WI;
				BGtrans=V{2}.BG;
                SDFI=V{3}.SDFI;
                SDFS=V{3}.SDFS;
                errmat(i,:)=fx'./normvec;
                solmat(i,:)=soltmp';
                condmat(i,:)=model.DSGEModel.structToVec(condvars)';
                Wshtrans(i,:) = [KBtrans', LBtrans', WItrans', BGtrans'];
                SDFmat(i,:) = [SDFI',SDFS'];
            end
            
        end
        
        % simulate model; overwrite method from superclass 
         function [simseries,varnames,errmat,Wshtrans,SDFmat]=simulate(obj,NT,NTini,inistvec,simerror,shmat)
            if length(inistvec)~=obj.Vfct.SSGrid.Ndim
                error('inistvec must be vector of length SSGrid.Ndim');
            end
            
            NTtot=NT+NTini;
            Nvartot=1+obj.NSTEX+obj.NSTEN+obj.NSOL+obj.NV+obj.NADD;
            if obj.NZNS>0
               Nvartot=Nvartot+obj.NZNS;
            end
            simseries=zeros(NTtot,Nvartot);
            
            % if shock matrix wasn't passed in, create it
            preset_path = false;
            if nargin<6
                rng(10,'twister');
                shmat=lhsdesign(NTtot,1,'criterion','correlation'); % realizations of shocks for Markov state
            else 
                preset_path=isinteger(shmat);
            end        
            
            point=inistvec;
                       
            pointmat=zeros(NTtot,length(point));
            
            for t=1:NTtot
               pointmat(t,:)=point; 
               exst=point(1);
               
                % next period's exog. state
               if preset_path
                    exnext=shmat(t);
               else
                    transprob=cumsum(obj.Exogenv.mtrans(exst,:));
                    exnext=find(transprob-shmat(t)>0,1,'first');
               end
               
               % transition to next period
               solvec=obj.evaluatePol(point);
               valvec=obj.evaluateVal(point)';
               znsvec=[];
			   if obj.NZNS>0
					znsvec=obj.Zfct.evaluateAt(point)';
			   end
               [nextst,outstr]=obj.calcStateTransition(point,solvec,exnext);
               
               addvec=model.DSGEModel.structToVec(outstr.addvars)';
               % write different categories of variables in one row
               simnext=[point(1),outstr.exstvec',point(2:end),solvec',valvec,addvec];
               if ~isempty(znsvec)
                   simnext=[simnext,znsvec];
               end
               if length(simnext)~=size(simseries,2)
                    disp('problem');
               end
               simseries(t,:)=simnext;
               point=nextst;
            end     
            
            simseries=simseries(NTini+1:end,:);
            varnames=[{'exst'}, obj.Ex_names, obj.En_names, obj.Sol_names, ...
				obj.V_names, obj.Add_names, obj.Zns_names];
            
            errmat=[];
            Wshtrans=[];%zeros(size(obj.Exogenv.mtrans(simseries(:,1),:)));
            SDFmat=[];%zeros(size(obj.Exogenv.mtrans(simseries(:,1),:)));
            if simerror
                [errmat,~,condmat,Wshtrans,SDFmat]=obj.calcEEError(pointmat);
                errmat=errmat(NTini+1:end,:);
                condmat=condmat(NTini+1:end,:);
                Wshtrans=Wshtrans(NTini+1:end,:);
                SDFmat=SDFmat(NTini+1:end,:);
                simseries=[simseries,condmat];
                varnames=[varnames,obj.Cond_names];
            else
                simseries=[simseries,zeros(NT,obj.NCOND)];                
                varnames=[varnames,strcat(obj.Cond_names,'_nan')];
            end
		 end
		 
% 		 function [startpt_vec,simnext] = computeFirstPeriod(obj,start_ini,enstatevec,shockstruct)
% 			startpt.KB=enstatevec * (1 + shockstruct.KB);
% 			startpt.LB=enstatevec * (1 + shockstruct.LB);
% 			startpt.WI=enstatevec * (1 + shockstruct.WI);
% 			startpt.BG=enstatevec * (1 + shockstruct.BG);
% 			startpt=orderfields(startpt,obj.En_names);
% 			startpt_vec=model.DSGEModel.structToVec(startpt)';
% 			startpt_vec=[start_ini,startpt_vec];
% 
% 			if strcmp(shockstruct.MITshock,'DWL')
% 				MITshock = enstatevec(2)*shockstruct.LB + ...
% 					enstatevec(3)*shockstruct.WI + enstatevec(4)*shockstruct.BG;
% 
% 				solguess=obj.evaluatePol(startpt_vec);
% 				transguess=obj.evaluateTrans(startpt_vec);
% 				valguess=obj.evaluateVal(startpt_vec);
% 				guess=[solguess;transguess([1,(exnpt+1):(3*exnpt+1)]);valguess(6)];
% 				objfun = @(x)computeMITShockState(startpt_vec,x,obj,MITshock);
% 
% 				options=optimset('Display','off','TolX',1e-15,'TolFun',1e-12,...
% 						'MaxIter',100,'MaxFunEvals',100^2,'FinDiffType','central');
% 				[sol,~,exfl] = fsolve(objfun,guess,options);
% 				if exfl<1
% 					warning('No Eqm Found in Run %f',n);
% 				end
% 				[~,simnext]=objfun(sol);
% 				simnext = [simnext,nan(1,obj.NCOND)];
% 			else
% 				simnext = [];
% 			end
% 		end
        
        function [mobj,failedPoints,dist,distT,distmat,distTmat]=polIter(mobj,MAXIT,revisitFailed,printmode,damp,tol_avg,enforceTFbounds)
            gridSt=mobj.Vfct.SSGrid.Pointmat; % use points from BaseGrid here
            NPT=mobj.Vfct.SSGrid.Npt;
            NDIM=mobj.Vfct.SSGrid.Ndim;
            exnpt=size(mobj.Exogenv.pts_perm,1);
            
            % initialize
            resmat=mobj.evaluatePol(gridSt)';
            resmat_prev=resmat;
            
            % split up matrix of points for better output
            gr_points=cell(exnpt,1);
            gr_index=cell(exnpt,2);
            for i=1:exnpt
                grinlog=(gridSt(:,1)==i);
                grind=find(grinlog);
                gr_points{i}=gridSt(grinlog,:);
                gr_index{i,1}=grinlog;
                gr_index{i,2}=grind;
            end
            
            % value function
            VF=mobj.evaluateVal(gridSt)';
            VFnext=zeros(size(VF));
            TF=mobj.evaluateTrans(gridSt)';
            TFnext=zeros(size(VF,1),mobj.Tfct.Nof);
                        
            % control flags
            iter=0;
                        
            disp(' ');
            disp('Starting main loop ...');
            disp(' ');
            while 1
                % counter
                iter=iter+1;
                
                % ===========================================
                % loop over state space
                % ===========================================
                
                % matrix for failed points
                failedPoints=[];
                % transitions
                transmat=mobj.evaluateTrans(gridSt)';
				vmat=mobj.evaluateVal(gridSt)';
				vi=vmat(:,6);
                
                failedPoints_trans_T=zeros(0,size(transmat,2));
                failedPoints_trans_I=[];
                failedPoints_trans_V=zeros(0,size(VF,2));
                
                % a rectangular grid to speed up VF interpolations solutions           
                resmat_startiter = resmat; % for dampening
                % outer loop: all exogenous states
                for ei=1:exnpt
                    tmp_grid=gr_points{ei};
                    tmp_indlog=gr_index{ei,1};
                    tmp_index=gr_index{ei,2};
                    tmp_resmat=resmat(tmp_indlog,:);
                    tmp_resmat_prev=resmat_prev(tmp_indlog,:);
                    
                    tmp_transmat=transmat(tmp_indlog,:);
					tmp_vi = vi(tmp_indlog);
                    tmp_NPT=size(tmp_index,1); % nb SS pts for exog state ei
                    
                    % Evaluate value functions at transition points
                    transpts=reshape([repmat(1:exnpt,tmp_NPT,1),tmp_transmat],tmp_NPT*exnpt,NDIM);
                    Vtrans = mobj.evaluateVal(transpts)';
                    
                    % index matching for transitions
                    refidx=kron((1:tmp_NPT)',ones(1,exnpt));
                    refidx=refidx(:);                    
                    
                    
                    disp(['State ',num2str(ei)]);

                    [tmp_resmat_new,tmp_VF,tmp_TF,tmp_failed]=mobj.solvePointList(...
						tmp_grid,tmp_resmat,tmp_transmat,tmp_vi,...
                        refidx,Vtrans,printmode,[]);
                    failedPoints=[failedPoints; tmp_index(tmp_failed)];
                    if revisitFailed
                        failedPoints_trans_T=[failedPoints_trans_T; tmp_transmat(tmp_failed,:)];
                        refidx_failed = ismember(refidx,find(tmp_failed));
                        
                        failedPoints_trans_I=[failedPoints_trans_I; tmp_index(refidx(refidx_failed))];
                        failedPoints_trans_V=[failedPoints_trans_V; Vtrans(refidx_failed,:)];
                    end
                                       
                    resmat_prev(tmp_indlog,:)=tmp_resmat;
                    resmat(tmp_indlog,:)=tmp_resmat_new;
                    VFnext(tmp_indlog,:)=tmp_VF;
                    TFnext(tmp_indlog,:)=tmp_TF;
                end                          
                
                if (revisitFailed && ~isempty(failedPoints))
                    disp( '~~~~~~~~~~~~~~~~~~~');
                    disp(['Revisiting failed points: ',num2str(length(failedPoints)),' add. points ...']);
                    % try to solve at failed points
                    [new_resmat,new_VF,new_TF,n_succ]=mobj.solvePointListFailed(gridSt,failedPoints,resmat,vi,...
                        failedPoints_trans_T,failedPoints_trans_I,failedPoints_trans_V,...
                        1,printmode,[]);
                    resmat(failedPoints,:)=new_resmat;
                    VFnext(failedPoints,:)=new_VF;
                    TFnext(failedPoints,:)=new_TF;
                    disp(['Revisiting solved ',num2str(n_succ),' points.']);
				end     
                
                
                % approximate functions for next iteration
				if length(damp) > 1
					dampV = damp(1);
					dampT = damp(2);
				else
					dampV = damp;
					dampT = damp;
                end
                
                if enforceTFbounds
                    nsten=mobj.NSTEN;
                    tvals=TFnext;
                    for is=1:nsten
                        thisbounds=mobj.Tfct.SSGrid.StateBounds(:,is+1);
                        thisvals=tvals(:,1+(is-1)*exnpt:is*exnpt);
                        thisvals(thisvals<thisbounds(1))=thisbounds(1);
                        thisvals(thisvals>thisbounds(2))=thisbounds(2);
                        tvals(:,1+(is-1)*exnpt:is*exnpt)=thisvals;
                    end
                    TFnext=tvals;
                end
                                
                
                mobj=mobj.updateVfct((1-dampV)*VFnext+dampV*VF);
                mobj=mobj.updateTfct((1-dampT)*TFnext+dampT*TF);
                
                % convergence criterion (based on points in BaseGrid)
                val_range=4:6;
                VF_val = VF(:,val_range);
                VFnext_val=VFnext(:,val_range);
				distmat=abs(VF_val-VFnext_val);
                [dist,wh]=max(distmat(:));
                [mean_dist,col]=max(mean(distmat));
				prct=99;
				[prct_dist,prct_col]=max(prctile(distmat,prct));
				distTmat=abs(TF-TFnext);
                [distT,whT]=max(distTmat(:));                
                [wh_1,wh_2]=ind2sub(size(VFnext_val),wh);
                [whT_1,whT_2]=ind2sub(size(TFnext),whT);
                disp(['-- Iteration: ',num2str(iter),', max distance: ',num2str(dist),' in ',char(mobj.V_names(val_range(1)-1+wh_2)), ...
                    ' at point ',num2str(wh_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(wh_1,:))]);
                disp(['-- Iteration: ',num2str(iter),', mean distance: ',num2str(mean_dist),' in ',char(mobj.V_names(val_range(1)-1+col))]);
                %disp(num2str(mobj.Vfct.SSGrid.Pointmat(wh_1,:)));
                disp(['-- Iteration: ',num2str(iter),', ',num2str(prct),'th percentile distance: ',num2str(prct_dist),' in ',char(mobj.V_names(val_range(1)-1+prct_col))]);
                %disp(num2str(mobj.Vfct.SSGrid.Pointmat(wh_1,:)));
                disp(['-- Iteration: ',num2str(iter),', max T distance: ',num2str(distT),' in col ',num2str(whT_2), ...
                    ' at point ',num2str(whT_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(whT_1,:))]);
                disp(' ');
                if mean_dist<tol_avg && prct_dist<2*tol_avg
                    disp('Converged.');
                    break;
                elseif iter>=MAXIT
                    disp('Max.iter. exceeded.');
                    break;
                end
                
                % update guess (based on points in BaseGrid)
                VF=VFnext;
                TF=TFnext; 
            end
            
            % resulting policy functions
            mobj=mobj.updatePfct(resmat);       
		end
		
		function obj = priceZNS( obj, stv, TOL, MAXIT, varargin )
			% Initialize
			params = obj.Params;
			betaB = params.betaB;
			betaS = params.betaS;
			g = params.g;
			nuB = params.nuB;
			nuS = params.nuS;
			exnpt = obj.Exogenv.exnpt;
			npt = obj.Vfct.SSGrid.Npt;
			if ismember('deltaT',fieldnames(params))
				deltaT = params.deltaT;
			elseif nargin==5
				deltaT=varargin{1};
				obj.Params.deltaT=deltaT;
			else
				error('deltaT undefined. Must recreate object or specify deltaT as argument');
			end
			
			
			betaS_g = betaS * exp(-g / nuS);
			betaB_g = betaB * exp(-g / nuB);

			% Intial guess for LT bond price
			qT_ss = betaS_g / (1 - betaS_g * deltaT);

			% Initial guess for cB stream price
			qCB_ss = betaB_g / (1 - betaB_g) * stv.Add.cB;

			% Initial guess for cB stream price
			qCS_ss = betaB_g / (1 - betaS_g) * stv.Sol.cS;
			pricenames = {'qT','qCB','qCS'};

			Zfct = grid.LinearInterpFunction( obj.Vfct.SSGrid, ...
				repmat( [qT_ss, qCB_ss, qCS_ss], npt, 1 ) );
			
			% Compute SDFs (once)
			SDFS = zeros(npt,exnpt);
			SDFB = zeros(npt,exnpt);
			transtens = zeros(1+obj.NSTEN,exnpt,npt);

			pointmat = obj.Vfct.SSGrid.Pointmat;
			vmat = obj.Vfct.Vals;
			transmat=obj.Tfct.Vals;
			nsten=obj.NSTEN;

			parfor i=1:npt
			%for i=1:npt
				pt=pointmat(i,:);
				trans=transmat(i,:);%trans=mobj.Tfct.Vals(i,:);
				vvec=vmat(i,:);
				cB=vvec(1);
				cS=vvec(2);
				transpts=reshape([(1:exnpt),trans],exnpt,nsten+1);%transpts=reshape([(1:mobj.Exogenv.exnpt),trans],mobj.Exogenv.exnpt,mobj.NSTEN+1);
				vtrans=obj.evaluateVal(transpts)';%vtrans=mobj.evaluateVal(transpts)';
				expStruct=obj.computeExpectations(pt(1),transpts,vtrans);

				% Pull SDF of a particular agent (here, saver)
				SDFS(i,:)=expStruct.Conditional.SDFS(:,1)' / cS^(-1/nuS);
				SDFB(i,:)=expStruct.Conditional.SDFB(:,1)' / cB^(-1/nuB);

				% Store transition states for interpolating future payoffs in the
				% iteration part later
				transtens(:,:,i) = transpts';
			end

			transmat = reshape(permute(transtens,[2,3,1]),[],nsten+1,1);

			% Iterate on payoffs
			iter=1;
			mtrans = obj.Exogenv.mtrans;

			while true
				resmat = Zfct.Vals;
				futurePrices = Zfct.evaluateAt( transmat )';
				parfor i=1:npt
				%for i=1:npt
					prnext = mtrans(pointmat(i,1),:);
					SDFBtmp = SDFB(i,:)';
					SDFStmp = SDFS(i,:)';
					%futurePrice = Pricefct.evaluateAt( transtens(:,:,i)' )';
					futurePrice = futurePrices(exnpt*(i-1)+(1:exnpt),:);

					% LT Bond
					payoff = 1 + deltaT * futurePrice(:,1);
					qT = prnext * ( SDFStmp .* payoff );

					% Consumptions
					vvec = vmat(i,:);

					% cB
					cB = vvec(1);
					payoff = cB + futurePrice(:,2);
					qCB = prnext * ( SDFBtmp .* payoff );

					% cS
					cS = vvec(2);
					payoff = cS + futurePrice(:,3);
					qCS = prnext * ( SDFStmp .* payoff );

					resmat(i,:) = [qT, qCB, qCS];
				end

				[dist,wh]=max(abs(resmat(:)-Zfct.Vals(:)));     
				[wh_1,wh_2]=ind2sub(size(resmat),wh);
				disp(['-- Iteration: ',num2str(iter),', max distance: ',num2str(dist),' in ',char(pricenames(wh_2)), ...
					' at point ',num2str(wh_1),': ',num2str(obj.Vfct.SSGrid.Pointmat(wh_1,:))]);

				Zfct=Zfct.fitTo(resmat);

				if dist <= TOL
					disp('Converged!');
					break;
				elseif iter>=MAXIT
					disp('Max.iter. exceeded.');
					break;
				end

				iter=iter+1;
			end
			
			obj.Zfct = Zfct;
			obj.Zns_names = pricenames;
			obj.NZNS = numel(pricenames);
		end
        
        function [simseries, varnames] = computeSimulationMoments(obj, simseries, varnames, SDFRmat, varargin)
%             if isempty(SDFRmat)
%                 SDFRmat=ones(size(simseries),1);
%             end
%             if isempty(VItildemat)
%                 VItildemat=ones(size(simseries));
%             end
			if nargin>4
				firstrow = varargin{1};
				simseries = [firstrow; simseries];
			end

            NT_sim = size(simseries,1);
            
            % make HashMap with mapping of names to indices
            indexmap=java.util.HashMap;
            for i=1:length(varnames)
                if isempty(indexmap.get(varnames{i}))
                    indexmap.put(varnames{i},i);
                end
            end
            
            %--------------------------------------------------------------------------
            % transformations of raw model output
            %--------------------------------------------------------------------------

            % list of indices
            loglist=model.HelperCollection.makeListFromNames(indexmap,{'qB','q','p','cI','cS','wB','wI','wS','cB'});
            multlist=model.HelperCollection.makeListFromNames(indexmap,{'lamR','lamS','muR','muS','lamB'});    

            % conversion of log-values
            simseries(:,loglist)=exp(simseries(:,loglist));
            % conversion of multipliers
            simseries(:,multlist)=max(simseries(:,multlist),0).^(1/3);


            % Function to calculate log growth rates
%             G=simseries(:,indexmap.get('G'));
%             loggr = @(var) diff(log(var)) + log(G(2:end));
            G=obj.Params.mu_G;
            loggr = @(var) diff(log(var)) + log(G);
            
            params = obj.Params;
            xi = simseries(:,indexmap.get('xi'));

            % ---------------------------------------------------------------------
            % state vars (also levels)
            % ---------------------------------------------------------------------
            ZA = simseries(:,indexmap.get('ZA')); %G = simseries(:,indexmap.get('G'));
            KB = simseries(:,indexmap.get('KB'));
            KB_g = simseries(:,indexmap.get('KB_g'));            
            WI= simseries(:,indexmap.get('WI'));
            BG=simseries(:,indexmap.get('BG'));

            % ---------------------------------------------------------------------
            % output, production, and investment
            % ---------------------------------------------------------------------
            deltaB=params.delta;
            deltaK=params.deltaK;
            qB=simseries(:,indexmap.get('qB'));
            MP=simseries(:,indexmap.get('MP'));
            OmA=simseries(:,indexmap.get('OmA'));
            OmK=simseries(:,indexmap.get('OmK'));
            X=simseries(:,indexmap.get('X'));
            wagebill=simseries(:,indexmap.get('wbill'));
            
            YB=simseries(:,indexmap.get('YB'));
            Y_om=simseries(:,indexmap.get('Y_om')); % == output if mu_om==1
            XbyK= X ./ KB;
            p=1 + params.psi*(X./KB-(exp(params.g)-1+deltaK));
            L=params.Lscale* params.Lbar(1)^params.gammaB * params.Lbar(2)^params.gammaS;
            MPK=(1-params.alpha)*params.mu_om*(KB/L).^(-params.alpha);
            MPL= params.alpha*params.mu_om*(KB/L).^(1-params.alpha);
            Kret = (p(2:end)*(1-deltaK) + MPK(2:end)/params.mu_om) ./ p(1:end-1) ./ G - 1;
            Kret = log(1+Kret);
%             AB=simseries(:,indexmap.get('AB'));
%             AB_g=simseries(:,indexmap.get('AB_g'));
			AI_g=simseries(:,indexmap.get('AI_g'));
			AS_g=simseries(:,indexmap.get('AS_g'));
            Psi=obj.Params.psi/2*(XbyK-(exp(params.g)-1+deltaK)).^2;
            DivP = simseries(:,indexmap.get('DivP'));
			eq_N = simseries(:,indexmap.get('equ')) ./ simseries(:,indexmap.get('N'));
			dP_rate = params.phi0 - eq_N - params.eqfct_intern*params.phi1/2*eq_N.^2;
			equ_total = simseries(:,indexmap.get('equ')) + (1-OmA)*params.n0;
            
%             div = (1-params.tauK)*(YB - OmA.*wagebill) - X - Psi.*KB ...
%                 - p.*KB_g + p.*(KB.*OmA*(1-(1-params.tauK)*deltaK) + X)  ...
%                 - (1-(1-params.theta)*params.tauK + params.delta*qB).*OmA.*AB + qB.*AB_g;
               
            
            % ---------------------------------------------------------------------
            % wealth distribution
            % ---------------------------------------------------------------------

            % Total corporate debt
            corpDebt = (MP+deltaB*OmA.*qB).*simseries(:,indexmap.get('AB'));
            WSm= corpDebt+BG-WI-params.BG_foreign./G;
            WBm = p.* KB - corpDebt;
			WB = (1 - simseries(:,indexmap.get('LB'))) .* p .* KB;
            PD = WBm./DivP;

            % ---------------------------------------------------------------------
            % interest rates and returns
            % ---------------------------------------------------------------------
            rB=log(1./qB+deltaB);
            %disp('xi:'); disp(xi)
            rB_lamI = log(1./(qB.*(1-simseries(:,indexmap.get('lamR')).*xi))+deltaB);
            Cspr = rB_lamI - rB;
            q=simseries(:,indexmap.get('q'));
            rD=1./q-1;
            Lspr=rB-rD;
			
			if ismember('qT',varnames)
				qT=simseries(:,indexmap.get('qT'));
				rT=log(1./qT+params.deltaT);
				LsprT=rB-rT;
				Tspr=rT-rD;
			else
				rT=nan(size(rD));
				LsprT=nan(size(Lspr));
				Tspr=nan(size(Lspr));
			end
			
            rD_lamI=[];%nan(size(simseries,1));
            if isempty(SDFRmat)==0
                q_lamI = sum(obj.Exogenv.mtrans(simseries(:,1),:).*SDFRmat,2);
                rD_lamI = 1./q_lamI-1;
            end
            % bond excess return
            theta=params.theta;
            tauPi=params.tauPi;
			shieldI=params.shieldI;
            Lret=(MP(2:end) + tauPi*(1-theta)*OmA(2:end) +deltaB*OmA(2:end).*qB(2:end))./qB(1:end-1)-rD(1:end-1)-1;
            Lret_tax=(MP(2:end) +deltaB*OmA(2:end).*qB(2:end))./qB(1:end-1)-rD(1:end-1)*(1-shieldI*tauPi)-1;
            Lret = log(1+Lret);
            Lret_tax = log(1+Lret_tax);
            
            % producer equity excess return
            %equ = simseries(:,indexmap.get('equ'));
            N = simseries(:,indexmap.get('N'));
            Eret = log( ( N(2:end) + DivP(2:end) ) ./ N(1:end-1)) + log(G) - rD(1:end-1);
            

            % ---------------------------------------------------------------------
            % corporate debt
            % ---------------------------------------------------------------------
            % total nominal debt and borrower LTV
            AB = simseries(:,indexmap.get('AB'));        
            AB_g = simseries(:,indexmap.get('AB_g'));
            bkdebt = params.F.*AB_g;
            assets = p.*KB_g;
            LTV=bkdebt./(p.*KB_g);
            tauK=params.tauK;
            tauD=params.tauD;
            mktLTV=qB.* AB_g./(p.*KB_g);
            mktNWB= p.*KB_g-qB.*AB_g;         
            % net loss rate on mortgages
            Drate=1-OmA;

            zeta=simseries(:,indexmap.get('Zeta'));
            LGD = 1 - ( (1-deltaK)*Drate.*(1-zeta).*p.*KB + (1-zeta*params.outputDWL).*(1-OmK).*ZA.*KB.^(1-params.alpha).*L.^params.alpha - Drate.*wagebill )./AB ...
                         ./ (Drate .* (1 + qB * params.delta));
            LGD_equity = p .* KB ./ (AB .* (1 + qB * params.delta) );
            % Losses from defaults 
            Lrate = LGD .* Drate;
            % growth in debt, at market prices
            mdebt=qB.*AB_g;            
            
            % borrower constraint
            bind_lamB=(simseries(:,indexmap.get('lamB'))>0);    
            Bbkdte=params.F*AB_g./N;
            
            % ---------------------------------------------------------------------
            % Intermediary debt, leverage, and payout
            % ---------------------------------------------------------------------
            bI = simseries(:,indexmap.get('bI'));
%             Imktlev=-q.*bI./mdebt;
%             Ibklev=-bI./mdebt;
            bIstart= WSm-BG;
            bI = -bI;
            mktNWI=qB.*AI_g-bI;
            
            %dbar = simseries(:,indexmap.get('dbar'));
            F_eps=simseries(:,indexmap.get('F_eps'));
            sigmaI=params.sigmaI*F_eps.^params.sigIexp;
            if params.sigmaI>0
                cI=simseries(:,indexmap.get('cI'));
                if params.eqfct_intern
                    eI=(1-cI)./sigmaI ;
                    dIcost_intern = sigmaI.*eI.^2/2;
                    dIcost_extern = 0;
                else
                    eI=(cI-1)./sigmaI ;
                    dIcost_intern = 0;
                    dIcost_extern = sigmaI.*eI.^2/2;
                end
            else
                eI=log(simseries(:,indexmap.get('cI')));
                cI= ones(size(simseries,1),1);                
                dIcost_intern = zeros(size(simseries,1),1);
                dIcost_extern = zeros(size(simseries,1),1);
            end
            if params.fixediv
                phiIWI=params.phiI*params.WIbar*ones(size(WI));
            else
                phiIWI=params.phiI*WI;                
            end
            dI_cost = dIcost_intern + dIcost_extern;
            F_eps_minus=simseries(:,indexmap.get('F_eps_minus'));
            F_eps_plus=simseries(:,indexmap.get('F_eps_plus'));
            brupt=1-F_eps;
            newequity=(1-F_eps).*WI;
            dI_eff=phiIWI - eI - dIcost_extern - F_eps_minus - newequity;
            bailout=F_eps_plus - (1-F_eps).*(WI - zeta.*corpDebt);
            payrate = phiIWI - eI - dIcost_extern - F_eps_minus; 
            WI_end = WI - phiIWI + eI - dIcost_intern;
            dI_rate = payrate./WI_end;
            eI_rate = eI./WI_end;
            eI_rate_L = eI(2:end)./WI_end(1:end-1);
            dI_rate_L = payrate(2:end)./WI_end(1:end-1);
			eI_total = eI + newequity;
            
            % Return on intermediary wealth, after taxes
            end_equity = ( AI_g(1:end-1) .* (MP(2:end) + deltaB*ZA(2:end).*qB(2:end)) - bI(1:end-1));
            start_equity = ( AI_g(1:end-1) .* qB(1:end-1) - (q(1:end-1)+shieldI*tauPi*rD(1:end-1)).*bI(1:end-1));

            ROW = end_equity ./ start_equity - 1;
            %ROW = log(ROW+1);

            % Limited liability: Set to -100% when risk-takers declare bankruptcy
            ROW(start_equity < 0) = 0;
            profit = end_equity - start_equity;
            % Set profit_t = -equity_ {t-1} if risk taker declares bankruptcy at time t
            ROA = Lret + rD(1:end-1);
                        

            % intermediary constraint
            bind_lamI = (simseries(:,indexmap.get('lamR')) > 0);

            % EER on loans based on conditional expectation (only works if
            % compEEErr ran)
            expRP=simseries(:,indexmap.get('expRP'));
            if ~isempty(expRP)
                expERP=expRP - 1./q;
            else
                expERP=nan(NT_sim,1);
            end            
            
            % expected excess return on intermediary equity
            VI=simseries(:,indexmap.get('VI'));            
            totalvalue=VI+q.*bI;
            debtcost= 1./(q + shieldI*tauPi*rD -params.kappa) -1;
            acctprof = (1-theta)*OmA.*AB - rD.*bI;
            franch = VI./WI_end-1;
            %franch(WI<0)=0;
            accROE = (params.F*AB(2:end)).*Lret./WI_end(1:end-1);
            %accROE(WI(2:end)<0)=0;
            
            expRIE= simseries(:,indexmap.get('retBE'));
            if ~isempty(expRIE)
                expERIE= expRIE - 1./q;
                expRIE = expRIE - 1;
                wacc = VI./totalvalue .* expRIE + (q.*bI)./totalvalue .* debtcost;
                profbility = rB - Lrate - wacc;
            else
                expERIE=nan(NT_sim,1);
                expRIE=nan(NT_sim,1);
                wacc = nan(NT_sim,1);
                profbility = nan(NT_sim,1);
            end            
            

                        
            % ---------------------------------------------------------------------
            % Consumption and welfare
            % ---------------------------------------------------------------------
            
            DWL=params.eta*((1-deltaK)*zeta.*Drate.*p.*KB + zeta.*params.outputDWL.*(1-OmK).*ZA.*KB.^(1-params.alpha).*L.^params.alpha) ...
                    + params.etaI*params.bankDWL*params.zetaI.*(1-F_eps).*corpDebt ...
                    + dI_cost;
            
            DWL_byY = DWL ./ Y_om;
            
            
            cB=simseries(:,indexmap.get('cB'));
            cS=simseries(:,indexmap.get('cS'));
            
            C = cB + cS;
            
            VB=simseries(:,indexmap.get('VB'));
            VS=simseries(:,indexmap.get('VS'));

%             PB=params.Pbar(1);
%             PS=params.Pbar(2);
%             welfare=PB*VB+PS*VS+PI*VI;
%             welfare=VB.^PB .* VS.^PS .* VI.^PI;
            welfare = VB + VS;

            % Risk sharing
            VCB = VB ./ cB;
            VCS = VS ./ cS;
            VC = sum( repmat(params.Pbar',length(VCB),1) .* [VCB, VCS], 2);
            deltalogVCB = diff(log(VCB));
            deltalogVCS = diff(log(VCS));
            deltalogVC = diff(log(VC));

            deltalogVCX = {deltalogVCB,deltalogVCS};
            regcoefs = zeros(2,1);
            for i=1:2
            RHS = [ones(size(deltalogVC)),deltalogVC];
            regY = deltalogVCX{i};
            tmp = (RHS'*RHS)\(RHS'*regY);
            regcoefs(i) = tmp(2);
            end

            MUB = (1 - params.betaB) * VCB;
            MUS = (1 - params.betaS) * VCS;

            MUBS = log(MUB ./ MUS);

            % Adjust output and investment
            Y = Y_om;
            X = X + Psi.*KB;
            XDWL = X + DWL;
            
            % ---------------------------------------------------------------------
            % government debt
            % ---------------------------------------------------------------------
%             rfdebt=simseries(:,indexmap.get('bS')) - min(BG,0);

            factor = min(1,ZA/params.mu_ZA); %factor=G / params.mu_G;
%            factor = ZA/params.mu_ZA; %factor=G / params.mu_G;
            taxrate = params.tau(BG./Y) .* factor .^ params.bTau;
            Ltax = wagebill.*taxrate;
            Ltax0 = wagebill.*params.tau0;
            Ktax = tauK*(YB-OmA.*(wagebill-params.pibar*KB));
            Pitax = tauPi*(1-theta)*OmA.*AB;
            Dtax = tauD*rD.*simseries(:,indexmap.get('bS'));         
            grosstax =  Ltax + Ktax + Pitax + Dtax;
            Ktax_ben = tauK*OmA.*((1-theta)*AB + params.deltaK*p.*KB);
            Pitax_ben = shieldI*tauPi*rD.*bI;
            Corptax = Ktax - Ktax_ben + Pitax - Pitax_ben + Dtax;
            tax = grosstax - Ktax_ben - Pitax_ben;
            T_byY = params.T .* factor .^ params.bT;
            T = T_byY .* Y;
            gammaG_byY = params.gammaG .* factor .^ params.bGamma;
            gammaG = gammaG_byY .* Y;
            spend = T + gammaG;
            
            Gdef_prim = (gammaG + T) - (tax);
            Gdef0_prim = Gdef_prim + Ltax - Ltax0;
            
            BG_g = simseries(:,indexmap.get('BG_g'));
                                    
            % Compute capital ratios
            X_byKB = X ./ KB;
            XDWL_byKB = XDWL ./ KB;
            
            % Compute output ratios
            Y_denom = Y;                        
            Y_om_byY = Y_om./Y_denom;
            YDWL = Y - DWL;
            
            gammaG_byY = gammaG./Y_denom;
            T_byY = T./Y_denom;
            Ltax_byY = Ltax./Y_denom;
            
            lgammaG_byY = log(gammaG_byY);
            lT_byY = log(T_byY);
            lLtax_byY = log(Ltax_byY);            
            
            % Growth rates
            grVar_table = table(mdebt,p,Y_om,Y,C,X,XDWL,DWL,KB,AB,bI,mktNWB,mktNWI,...
                gammaG,spend,tax,cB,cS,cI,mdebt,bkdebt,welfare,VB,VI,VS,...
                WI);
            grVar_varnames=grVar_table.Properties.VariableNames;
            grVar_table = varfun(loggr,grVar_table);
            grVar_table.Properties.VariableNames=strcat(grVar_varnames,'_gr');
            
            % HP filter
            logY = log(Y);
            logC = log(C);
            logX = log(X);
            logZA = log(ZA);
			logG = log(gammaG);
            % initial values (matters for trend normalization when calculating stats)
            %logZ(1) = 0; 
            %logY(1) = logZ(1)+log(Y(1)); %log(mean(Y)); 
            %logC(1) = logZ(1)+log(C(1)); %log(mean(C)); 
            %logX(1) = logZ(1)+log(X(1)); %log(mean(X));             
            %for i=1:NT_sim-1
            %    logY(i+1)=logY(i)+grVar_table.Y_gr(i);
            %    logX(i+1)=logX(i)+grVar_table.X_gr(i);
            %    logC(i+1)=logC(i)+grVar_table.C_gr(i);
            %    logZ(i+1)=logZ(i)+log(G);
            %end
			wins_cutoff = log(0.01);
			logX_wins = logX .* (logX > wins_cutoff ) + wins_cutoff .* (logX <= wins_cutoff );

                       
            % ---------------------------------------------------------------------
            % Leverage and equity (detrended)
            % ---------------------------------------------------------------------
            
			trend_factor = G; 
            F=obj.Params.F;
            
            % borrower-entrepreneurs (B)
            % book equity
            Bbkasset = KB_g ./ trend_factor;
            Bbklbt = F.*AB_g ./ trend_factor;
            Bbklbt_byY = Bbklbt ./ Y_denom;
            Bbkeqt = Bbkasset-Bbklbt;
            % market equity
            Bmasset =p.*KB_g ./ trend_factor; 
            Bmlbt = qB.*AB_g ./ trend_factor; 
            Bmeqt = Bmasset-Bmlbt;
            % book leverage
            Bbklvg = Bbklbt./Bbkasset;
            % market leverage
            Bmlvg = Bmlbt./Bmasset;
            
            % financial intermediary (I)
            % book equity
            Ibkasset = F.*AI_g ./ trend_factor;
            Ibklbt = bI ./ trend_factor;
			Ibklbt_byY = Ibklbt ./ Y_denom;
            Ibkeqt = Ibkasset-Ibklbt;
            % market equity
            Imasset = qB.*AI_g ./ trend_factor;
            Imlbt = q.*(bI) ./ trend_factor;
            Imeqt = Imasset-Imlbt;
            IeIbymasset = eI ./ Imasset;
            IeIbybkasset = eI ./ Ibkasset;
            % book leverage
            Ibklvg = Ibklbt./Ibkasset;
            % market leverage
            Imlvg = Imlbt./Imasset;
            
            % government
            Gmlbt = q.*BG_g;
						
			% Replace intermediary leverage with nan when it's effectively
			% 0/0
			Ibklvg( abs(Ibklbt) < 1e-5 & abs(Ibkasset) < 1e-5 ) = nan;
			Imlvg( abs(Imlbt) < 1e-5 & abs(Imasset) < 1e-5 ) = nan;
			
			% fraction corporate debt held by savers
			fracS = AS_g ./ AB_g;
			d_fracS = diff(fracS);
			                        
            WI_detrend = WI./trend_factor;
                      
            % add to simulation and names vectors
            if isempty(rD_lamI)
                rD_lamI=nan(size(Cspr));
            end          
            
            % ratios to output
            byVar_table = table(WI,WBm,X,XDWL,KB,wagebill,assets,mdebt,bkdebt,DWL,...
                cB,cI,cS,bI,bIstart,eI,Imasset,tax,Ltax,Ltax0,Ktax,Pitax,Dtax,Ktax_ben,...
                Pitax_ben,spend,Corptax,Gdef_prim,Gdef0_prim,BG_g,C,Y,T,gammaG);
            byVar_varnames=byVar_table.Properties.VariableNames;
            byVar_table = varfun(@(var)var./Y_denom,byVar_table);
            byVar_table.Properties.VariableNames=strcat(byVar_varnames,'_byY');
            fracProfligacy = ( byVar_table{:,end-4} < params.BG1); 
            fracAusterity = (byVar_table{:,end-4} > params.BG2);
            
            
            statevec=simseries(:,1);
            simseries=[simseries(:,2:end),... state vars
                                      WB, WBm, WSm, ... wealth distribution
                                      XDWL, X_byKB, XDWL_byKB , p, PD, MPK, MPL, ... production
                                      rD, rB, Lspr, rT, LsprT, Tspr, rB_lamI, Cspr, rD_lamI, ...
                                      mdebt, LTV, mktLTV, Drate, LGD, LGD_equity, Lrate, ...
                                      DWL_byY, DivP, dP_rate, fracS, bind_lamB, ... corporate debt
                                      brupt,bind_lamI,... intermediary constraint
                                      eI, dI_eff, dI_cost, dI_rate, eI_rate, newequity, F_eps,...
									  payrate, eI_total, equ_total, ...
                                      welfare, MUBS, ... consumption and welfare
                                      bailout, expERP, expRIE, expERIE, wacc, profbility, acctprof, franch, ...
                                      fracProfligacy, fracAusterity, Y, YDWL, Y_om_byY, C, gammaG, Gmlbt, ...
                                      lgammaG_byY, lT_byY, lLtax_byY, logY, logX, logC, logZA, logG];
            simseries = [simseries, byVar_table{:,:}];
            simseries = [simseries,...
            Bbkasset,Bbklbt,Bbklbt_byY,Bbkeqt,Bmasset,Bmlbt,Bmeqt,Bbklvg,Bmlvg,Bbkdte,... % firm equity and leverage
            Ibkasset,Ibklbt,Ibklbt_byY,Ibkeqt,Imasset,Imlbt,Imeqt,Ibklvg,Imlvg,IeIbymasset,IeIbybkasset,... % intermediary equity and leverage
            WI_detrend];

            simseries=[simseries(2:end,:),Kret,Eret,Lret,Lret_tax,profit,ROA,ROW,accROE,eI_rate_L,dI_rate_L];
            simseries = [simseries,grVar_table{:,:},d_fracS];

            varnames_add=[{'WB', 'WBm','WSm',... wealth distribution
                      'XDWL','X_byKB','XDWL_byKB', 'p', 'PD', 'MPK', 'MPL',... production
                      'rD', 'rB', 'Lspr', 'rT', 'LsprT', 'Tspr', 'rB_lamI', 'Cspr', 'rD_lamI', ...
                      'mdebt', 'LTV', 'mktLTV', 'Drate', 'LGD', 'LGD_equity', 'Lrate', ...
                      'DWL_byY', 'DivP', 'dP_rate', 'fracS', 'bind_lamB', ... corporate debt
                      'brupt','bind_lamI',... intermediary constraint
                      'eI', 'dI_eff', 'dI_cost', 'dI_rate', 'eI_rate', 'newequity','F_eps',...
					  'payrate','eI_total','equ_total',...
                      'welfare', 'MUBS', ... consumption and welfare
                      'bailout', 'expERP', 'expRIE', 'expERIE', 'wacc', 'profbility', 'acctprof', 'franch', ...
                      'fracProfligacy', 'fracAusterity', 'Y', 'YDWL', 'Y_om_byY','C','gammaG','Gmlbt',...
                      'lgammaG_byY', 'lT_byY', 'lLtax_byY', 'logY', 'logX', 'logC', 'logZA', 'logG'},...
                      byVar_table.Properties.VariableNames,...
                      {'Bbkasset','Bbklbt','Bbklbt_byY','Bbkeqt','Bmasset','Bmlbt','Bmeqt','Bbklvg','Bmlvg','Bbkdte',... % firm equity and leverage
                      'Ibkasset','Ibklbt','Ibklbt_byY','Ibkeqt','Imasset','Imlbt','Imeqt','Ibklvg','Imlvg','IeIbymasset','IeIbybkasset',... % intermediary equity and leverage
                      'WI_detrend'},...
                      {'Kret','Eret','Lret','Lret_tax','profit','ROA','ROW','accROE','eI_rate_L','dI_rate_L'},...
                      grVar_table.Properties.VariableNames,'d_fracS'];
            varnames=[varnames(2:end), varnames_add];
            
        end
          
		function obj=toggleJacobian(obj)
			obj.Params.useJacobian = ~obj.Params.useJacobian;
		end
        
    end %of object methods
        
    
    %==============================================================================
    methods (Static)
        % static class-specific methods
        
        function [OmA,omplus,fom]=Mpayoff_norm(omstar,mu,sig2)
            omstar=real(omstar);
            sig=sqrt(sig2);
            OmA=1-normcdf(omstar,mu,sig);
            omplus=OmA.*mu+sig.*normpdf((omstar-mu)./sig);
            fom=normpdf(omstar,mu,sig);
        end
  
        function [OmA,omplus,fom]=Mpayoff_gamma(omstar,mu,sig2)
            chi1 = sig2./mu;
            chi0 = mu./chi1;
            omstar=real(omstar);
			
			z = omstar ./ chi1;
			OmA = 1 - gammainc(z, chi0);
			fom = omstar.^(chi0-1) .* exp(-z) ./ ( gamma(chi0) .* chi1.^chi0);
			omplus = mu.*(1-gammainc(z,chi0+1));
			
            %OmA=1-gamcdf(omstar,chi0,chi1);
            %fom=gampdf(omstar,chi0,chi1);
            %omplus=mu.*(1-gamcdf(omstar,chi0+1,chi1));
        end  
        
        function [Fx,Fxplus,fx]=quadapprox(omstar,eta2)
            eta0=1/2-eta2/3;
            fx = eta0 + eta2 * omstar.^2;
            Fx = eta0*omstar + eta2*omstar.^3/3 + eta0 + eta2/3;
            Fxplus = eta0*omstar.^2/2 + eta2*omstar.^4/4 - (eta0/2 + eta2/4);
        end
        
        
        function [J,fxout] = constructJacobian(params)
            syms ZA sig2B zeta xi AStarget ZA_factor
            syms KB LB WI BG
            
            syms qB X q equ cI cS wvec1 wvec2 AI_g AS_g cB
            syms muRplus muRminus
            syms lamIplus lamIminus
            syms lamSplus lamSminus
            syms lamBplus lamBminus
			syms muSplus muSminus
            
            syms exp_B_qB exp_B_qK exp_I_qB exp_I_q exp_S_q exp_S_qB VI
            
            params_numeric=params;
            params.tau=[];
            % unpack relevant parameters
            betaB=params.betaB;
            betaS=params.betaS;
            sigmaB=params.sigmaB;
            sigmaS=params.sigmaS;
            sigmaI=params.sigmaI;
%             xi=params.xi;
            mu_G=params.mu_G;
            mu_ZA=params.mu_ZA;
            kappa=params.kappa;
            theta=params.theta;
            tauK=params.tauK;
            psi=params.psi;
            tauPi=params.tauPi;
            tauD=params.tauD;
            alpha=params.alpha;
            deltaK=params.deltaK;
            delta=params.delta;
            Phi=params.Phi;
            F=params.F;
            mu_om=params.mu_om;
            try
               nuinv=1./[params.nuB,params.nuS];
            catch
               nuinv=1/params.nu*ones(3,1);
            end
            tauK_int=(1-theta)*tauK;
            tauPi_int=(1-theta)*tauPi;
            BG_foreign=params.BG_foreign;
            gammaB=params.gammaB;
            gammaS=params.gammaS;
            Lscale=params.Lscale;
            Lbar=params.Lbar;
            outputDWL=params.outputDWL;
            bankDWL=params.bankDWL;
            bTau=params.bTau;
            Pbar=params.Pbar;
            T=params.T;
            bT=params.bT;
            gammaG=params.gammaG;
            bGamma=params.bGamma;     
            eta=params.eta;
            sigma_eps=params.sigma_eps;
            pibar=params.pibar;
            mktConsI=params.mktConsI;
			shieldI=params.shieldI;
			phi0=params.phi0;
			phi1=params.phi1;
            n0=params.n0;
 			chiLinear=params.chiLinear;
            phiI=params.phiI;
            eqfct_intern=params.eqfct_intern;
            fixediv=params.fixediv;
			chi1adj=params.chi1adj;
            FLnodef=params.FLnodef;
            sigIexp=params.sigIexp;
			etaI=params.etaI;
			zetaI=params.zetaI;
 			           
            % capital price           
            p = 1 + psi*(X/KB-(mu_G-1+deltaK));
            % production and income
            Lbar_scale=Lbar.*Lscale;
            L=Lbar_scale(1)^gammaB * Lbar_scale(2)^gammaS;
            Y=mu_om*ZA*KB^(1-alpha) * L^alpha;
            Y_om=Y/mu_om;
            wagebill=Lbar_scale(1)*wvec1 +  Lbar_scale(2)*wvec2;
 
            % investment and capital
            Psi=psi/2*(X/KB-(mu_G-1+deltaK))^2;
            KB_g=(1 - deltaK + X/KB)*KB;
                        
            % borrower payoff risk
			AB=LB*p*KB/qB;
            %AB=(p*KB-WB)/qB;
            if FLnodef
                omstar=0;
            else                
                omstar=(AB + wagebill + pibar*KB)/Y_om;
            end
                         
            syms fOmA(omega) ffom(omega) fOmK(omega)
            OmA=fOmA(omstar);
            fom=ffom(omstar);
            OmK=fOmK(omstar);
            omplus=OmK;
            omminus=(1-OmK);
   
            % producer budget
            Gfun = @(om1,om2,tau,zetaDWL)( (1-zetaDWL*params.outputDWL)*(1-tau)*om1*Y_om - om2*(1-tau)*wagebill ...
                            + om2*( (1-zetaDWL)*(1-(1-tau)*deltaK)*p - (1-tau)*pibar)*KB );                            
            N = Gfun(omplus,OmA,tauPi,0) - OmA*(1-tauPi*(1-theta)+delta*qB)*AB  + (1-OmA)*n0;
            eq_N = equ/N;
            if eqfct_intern
                equcost_intern = phi1/2 * eq_N^2;
                equcost_extern = 0;
                equnorm = 1/(1-phi1*equ/N);
            else
                equcost_intern = 0;
                equcost_extern = phi1/2 * eq_N^2;
                equnorm = (1+phi1*equ/N);
            end            
            AB_g = (p*KB_g - (1-phi0)*N - equ + equcost_intern)/qB;
            Nminus = Gfun(omminus,1-OmA,0,zeta);
            DivP = N*(phi0 - eq_N - equcost_extern) - (1-OmA)*n0;
            nstar = (Gfun(omstar,1,tauPi,0) - (1-tauPi*(1-theta)+delta*qB)*AB )/N;
                    
            
            % payoff to intermediary
            MP= OmA*(1-tauPi_int) + Nminus/AB;
            YB=OmK*ZA*KB^(1-alpha)*L^alpha;                                               
            DebtB=(MP+delta*OmA*qB)*AB; 
            
            YnetDWL = Y;%Y - DWL;%Y
            
            % taxes
            syms tau(D_byY)
            taxrule_arg = BG/YnetDWL;
            

            % intermediary default
            syms fF_eps(V) fF_eps_cond(V)
            F_eps=fF_eps(VI);
            F_eps_cond=fF_eps_cond(VI);
            F_eps_minus=-sigma_eps*F_eps_cond;
            F_eps_plus=sigma_eps*F_eps_cond;
            sigmaI=sigmaI*F_eps^sigIexp;
            if sigmaI==0
                eI=cI;
                cI_norm=1;
                eIcost_intern = 0;
                eIcost_extern = 0;
            else
                if eqfct_intern
                    eI = (1-cI)/sigmaI;
                    eIcost_intern = sigmaI/2 * eI^2;
                    eIcost_extern = 0;
                    cI_norm = 1/cI;
                else
                    eI = (cI-1)/sigmaI;                
                    eIcost_intern = 0;
                    eIcost_extern = sigmaI/2 * eI^2;
                    cI_norm = cI;
                end
            end
            
            if fixediv
                phiIWI=phiI*params.WIbar;
            else
                phiIWI=phiI*WI;
            end            
            
            dI_eff = phiIWI - eI - eIcost_extern - F_eps_minus - (1-F_eps)*WI;
            bailout = F_eps_plus - (1-F_eps)*(WI - bankDWL*zetaI*DebtB);            
			
			if chi1adj>0
				PsiSA = chi1adj/2 * (AS_g / AStarget - 1)^2 * AStarget;
				dPsiSA = chi1adj * (AS_g / AStarget - 1);
			else
				PsiSA = 0;
				dPsiSA = 0;
			end
			PsiS = PsiSA + chiLinear * AS_g;
            PsiP = N*phi1/2*eq_N^2;
            BR_expense_all = (1-eta)*(zeta*(1-OmA)*(1-deltaK)*p*KB + zeta*outputDWL*(1-OmK)*ZA*KB^(1-alpha)*L^alpha) ...
                              + bankDWL*(1-etaI)*(1-F_eps)*zetaI*DebtB + PsiS + PsiP + pibar*KB;
            BR_expense1 = Pbar(1) * BR_expense_all;
            BR_expense2 = Pbar(2) * BR_expense_all;
            
            
            % government debt and taxes
            tauN=tau(taxrule_arg); % + tau0;
            WS=DebtB+BG-WI-BG_foreign*BG;
            tau_eff=tauN * ZA_factor^bTau;
            Tvec=Pbar .* T * YnetDWL * ZA_factor^bT;
            public_goods = gammaG * YnetDWL * ZA_factor^bGamma;
            
            % budget constraint borrower
            netResourcesB=(1-tau_eff)*wvec1*Lbar_scale(1) + Tvec(1) + BR_expense1 + p*X + dI_eff + DivP - X - Psi*KB;
                                 
            % intermediary
            rD=1/q-1;
            bI=( WI - phiIWI + eI - eIcost_intern - qB*AI_g)/(q+shieldI*tauPi*rD-kappa);            
            % saver
            bS=( (1-tau_eff)*wvec2*Lbar_scale(2) + WS + Tvec(2) + BR_expense2 - cS - qB*AS_g - PsiS)/(q+tauD*rD);          
           
            % next period government debt             
            spending = public_goods + sum(Tvec) + bailout;
            taxes = (tau_eff-tauK*OmA)*wagebill + tauK*YB + OmA*AB*(tauPi_int-tauK_int) + bI*(shieldI*tauPi*rD-kappa) + bS*tauD*rD - tauK*OmA*(deltaK*p+pibar)*KB;
            BG_g = (BG + spending - taxes)/q;           
            
            % compute (sym) current U's to normalize back SDFs
            nuinv_vec=nuinv;            
            %cvec=[cB,cS,cI];
            UB=cB;
            U_vecB=UB;
            U_vecS=cS;
            U_normB=U_vecB^(1-nuinv_vec(1))/cB;  
            U_normS=U_vecS^(1-nuinv_vec(2))/cS;              
            U_normI=U_normB*cI_norm; 
            U_normP=U_normB*equnorm;
            
            % borrower FOCs
            fx1 = qB - F*lamBplus - exp_B_qB/U_normP;
            fx2 = p - Phi*p*lamBplus - exp_B_qK/U_normP;
            MPLvec1=gammaB/Lbar_scale(1)*alpha*mu_om*ZA*L*(KB/L)^(1-alpha);
            MPLvec2=gammaS/Lbar_scale(2)*alpha*mu_om*ZA*L*(KB/L)^(1-alpha);
            Fscript=nstar/Y_om;
            fx3=(1-tauK)*OmK*MPLvec1/mu_om - (1-tauK)*OmA*wvec1 - fom*Fscript*(wvec1 - omstar*MPLvec1);
            fx4=(1-tauK)*OmK*MPLvec2/mu_om - (1-tauK)*OmA*wvec2 - fom*Fscript*(wvec2 - omstar*MPLvec2);
            
            % risk taker FOC
            rD=1/q-1;
            fx5=qB*(1 - xi*lamIplus) - muRplus - exp_I_qB/U_normI;
            if mktConsI
                fx6=q + shieldI*tauPi*rD - q*lamIplus - kappa - exp_I_q/U_normI;
            else
                fx6=q + shieldI*tauPi*rD - lamIplus - kappa - exp_I_q/U_normI;
            end
                                      
            % saver FOC
            fx7=q + tauD*rD -lamSplus - exp_S_q/U_normS;
            
            % risk taker constraints
            if mktConsI
                fx8=xi*qB*AI_g + q*bI - lamIminus; 
            else
                fx8=xi*qB*AI_g + bI - lamIminus; 
            end
            
            % saver constraint
            fx9=bS - lamSminus;
            
            % risk free debt market clearing
            fx10=BG_g - bS - bI - BG_g*BG_foreign;
            % Intermediary no-shorting
            fx11=AI_g - muRminus;
            			
			% corporate bond market
			fx12=qB - muSplus + dPsiSA + chiLinear - exp_S_qB/U_normS;
			fx13=AS_g - muSminus;
			fx14=AB_g - AI_g - AS_g;
			fx15=cB - netResourcesB;
			
			% borrower hard constraint
            fx16=Phi*p*KB_g - F*AB_g - lamBminus;
			
            
            fx=[fx1;fx2;fx3;fx4;fx5;fx6;fx7;...
                fx8;fx9;fx10;fx11;fx12;fx13;fx14;fx15;fx16];
            
            Jacobian = jacobian(fx,[qB,X,q,equ,cI,cS,wvec1,wvec2,...
				AI_g,AS_g,cB,...
                muRplus,muRminus,lamIplus,lamIminus,...
                lamSplus,lamSminus,...
				muSplus,muSminus,...
                lamBplus,lamBminus]);

            syms OmA fom OmK Tau dOmA dOmK dTau dfom dOmK F_eps F_eps_cond dF_eps dF_eps_cond
            sym_lookup_fOmA = subs(fOmA(omega),omega,omstar);
			sym_lookup_ffom = subs(ffom(omega),omega,omstar);
			sym_lookup_fOmK = subs(fOmK(omega),omega,omstar);
			sym_lookup_fF_eps = subs(fF_eps(V),V,VI);
			sym_lookup_fF_eps_cond = subs(fF_eps_cond(V),V,VI);
			sym_lookup_Tau = subs(tau(D_byY),D_byY,taxrule_arg);
			
			% Derivatives require a two-step process to force MATLAB to
			% recognize the chain rule
            % Example for D(fOmA)
%             (1) compute diff(fOmA,2*omega)/2 to get exactly D(fOmA)(2*omega),
%             (2) replace 2*omega with whatever the complicated expression for omstar is to get "D(fOmA)(bunch of stuff)"
%             (3) replace  "D(fOmA)(bunch of stuff)" with "dfOmA" in the Jacobian
			chain=2;
			sym_lookup_dOmA = diff(fOmA(chain*omega),omega)/chain;
			sym_lookup_dOmA = subs(sym_lookup_dOmA,chain*omega,omstar);
			sym_lookup_dTau = diff(tau(chain*D_byY),D_byY)/chain;
			sym_lookup_dTau = subs(sym_lookup_dTau,chain*D_byY,taxrule_arg);
			sym_lookup_dfom = diff(ffom(chain*omega),omega)/chain;
			sym_lookup_dfom = subs(sym_lookup_dfom,chain*omega,omstar);
			sym_lookup_dOmK = diff(fOmK(chain*omega),omega)/chain;
			sym_lookup_dOmK = subs(sym_lookup_dOmK,chain*omega,omstar);
			sym_lookup_dF_eps = diff(fF_eps(chain*V),V)/chain;
			sym_lookup_dF_eps = subs(sym_lookup_dF_eps,chain*V,VI);
            sym_lookup_dF_eps_cond = diff(fF_eps_cond(chain*V),V)/chain;
			sym_lookup_dF_eps_cond = subs(sym_lookup_dF_eps_cond,chain*V,VI);
			
			Jacobian=subs(Jacobian,sym_lookup_dOmA,'dOmA');
			Jacobian=subs(Jacobian,sym_lookup_dTau,'dTau');
			Jacobian=subs(Jacobian,sym_lookup_dfom,'dfom'); 
			Jacobian=subs(Jacobian,sym_lookup_dOmK,'dOmK'); 
			Jacobian=subs(Jacobian,sym_lookup_dF_eps,'dF_eps');
			Jacobian=subs(Jacobian,sym_lookup_dF_eps_cond,'dF_eps_cond');
			Jacobian=subs(Jacobian,sym_lookup_Tau,'Tau');
			Jacobian=subs(Jacobian,sym_lookup_fOmA,'OmA');
			Jacobian=subs(Jacobian,sym_lookup_ffom,'fom'); 
			Jacobian=subs(Jacobian,sym_lookup_fOmK,'OmK');  
			Jacobian=subs(Jacobian,sym_lookup_fF_eps,'F_eps');
			Jacobian=subs(Jacobian,sym_lookup_fF_eps_cond,'F_eps_cond');
			
			fx=subs(fx,sym_lookup_dOmA,'dOmA');
			fx=subs(fx,sym_lookup_dTau,'dTau');
			fx=subs(fx,sym_lookup_dfom,'dfom'); 
			fx=subs(fx,sym_lookup_dOmK,'dOmK'); 
			fx=subs(fx,sym_lookup_dF_eps,'dF_eps');
			fx=subs(fx,sym_lookup_dF_eps_cond,'dF_eps_cond');
			fx=subs(fx,sym_lookup_Tau,'Tau');
			fx=subs(fx,sym_lookup_fOmA,'OmA');
			fx=subs(fx,sym_lookup_ffom,'fom'); 
			fx=subs(fx,sym_lookup_fOmK,'OmK');  
			fx=subs(fx,sym_lookup_fF_eps,'F_eps');
			fx=subs(fx,sym_lookup_fF_eps_cond,'F_eps_cond');
			
			%Jacobian = simplify(Jacobian);
           
            order=[ZA, sig2B, zeta, xi, AStarget, ZA_factor, KB, LB, WI, BG,...
                qB, X, q, equ, cI, cS, wvec1, wvec2,AI_g,AS_g,cB,...
                muRplus,muRminus,lamIplus,lamIminus,...
                lamSplus,lamSminus,...
				muSplus,muSminus,...
                lamBplus,lamBminus,...
                OmA, fom,OmK, Tau, dOmA, dTau, dfom, dOmK,...
                F_eps, F_eps_cond, dF_eps, dF_eps_cond,...
                exp_B_qB, exp_B_qK, exp_I_qB, exp_I_q, VI, exp_S_q, exp_S_qB];
            J=matlabFunction(Jacobian,'Vars',order,'File','jacfun'); % convert symbolic function to numeric Matlab function (Jacobian)
            fxout=matlabFunction(fx,'Vars',order,'File','fxfun'); % convert symbolic function to numeric Matlab function (residuals)   
        end        
                
      function [fx,stvals]=compStSt(sol,params,omstate,print)
            betaB=params.betaB;
            betaS=params.betaS;
            gammaB=params.gammaB;
            gammaS=params.gammaS;
            try
                nuB=params.nuB;
                nuS=params.nuS;
            catch
                nuB=params.nu;
                nuS=params.nu;
            end
%             xi=params.xi;
            kappa=params.kappa;
            mu_ZA=params.mu_ZA;
            g=params.g;
            mu_G=params.mu_G;
            tau=params.tau;
            theta=params.theta;
            tauK=params.tauK;
            tauPi=params.tauPi;
            tauD=params.tauD;
            alpha=params.alpha;
            deltaK=params.deltaK;
            Lbar=params.Lbar;
            delta=params.delta;
            zeta=params.zeta(omstate);
            eta=params.eta;
            mu_om=params.mu_om;
            sig2_om=params.sig2_om(omstate);
            sig_om=sqrt(sig2_om);
            xi = params.xivec(2); 
            F=params.F;
            Phi=params.Phi;
            gammaG=params.gammaG;
            Pbar=params.Pbar;
            sigma_eps=params.sigma_eps;
            rho=params.rho;
            pibar=params.pibar;
            bankDWL=params.bankDWL;
			shieldI=params.shieldI;
            phi0=params.phi0;
            phi1=params.phi1;
            n0=params.n0;
			chiLinear=params.chiLinear;
            phiI=params.phiI;
            sigmaI=params.sigmaI;
            eqfct_intern=params.eqfct_intern;
			chi1adj=params.chi1adj;
			fracStarget=params.fracStarget(omstate);
			FLnodef=params.FLnodef;
            sigIexp=params.sigIexp;
			etaI=params.etaI;
			zetaI=params.zetaI;
			
            tauK_int=(1-theta)*tauK;
            tauPi_int=(1-theta)*tauPi;
            
            mktConsI=params.mktConsI;
                      
            % solution vars
            AB_g=exp(sol(1));
            KB_g=exp(sol(2));
            wB=exp(sol(3));
            wS=exp(sol(4));
            qB=exp(sol(5));
            eI=sol(6);
            VI=exp(sol(7));
            lamB=sol(8);
            equ=sol(9);
            
			% 5-dim sol: no AS_g, Lscale pre-set
			% 6-dim sol: AS_g and Lscale pre-set OR no AS_g and Lscale
			% 7-dim sol: AS_g and Lscale
		
			if numel(sol)==11 || (numel(sol)==10 && chiLinear > 0)
				if numel(sol)==11
					Lscale=exp(sol(11));
					LscaleIdx = 11;
				else
					Lscale=exp(sol(10));
					LscaleIdx = 10;
				end
				take_Lscale_as_given=0;
			else
				Lscale=params.Lscale;
				take_Lscale_as_given=1;
			end
			
			if chiLinear > 0
				fracS = 0;
			else
				fracS=exp(sol(10))/(1+exp(sol(10)));
            end                        
            
            % growth
            G=mu_G; %=exp(g)

            lamBplus=max(0,lamB)^3;
            lamBminus=max(0,-lamB)^3;            

            % short-term bond prices
            betaS_g=betaS*exp(-g/nuS);
            betaB_g=betaB*exp(-g/nuB);
            q=(tauD+betaS_g)/2 + sqrt( (tauD+betaS_g)^2/4 - tauD );
            %q=betaS_g;
            rD=1/q-1;
            
            % capital price
            p=1;
            
            % government budget and debt
            BG_g=params.BG_guess*(1-params.BG_foreign); % just as a guess (does not matter for prices)
            tauB=tau(BG_g);
            tauI=tauB;
            tauS=tauB;         
            
            % sum of wages
            Lbar = Lscale*Lbar;
            L = Lbar(1)^gammaB * Lbar(2)^gammaS;
            wagevec=[wB,wS];
            wagebill = wagevec*Lbar;
            
            % normalized output
            Y=mu_om*mu_ZA*(KB_g/G)^(1-alpha)*L^alpha;
            Y_om=Y/mu_om;
            
            % entrepreneur idiosyncratic shocks
            % default threshold
            omstar=(AB_g/G + wagebill + pibar*KB_g)/Y_om;
            if FLnodef
                omstar=0; % turn of liquidity default for firms
            end
            % OmA and OmK
            [OmA,OmK,f_om]=ProductionModel.Mpayoff_gamma(omstar,mu_om,sig_om^2);
            omplus=OmK;
            omminus=(1-OmK);
            
            % intermediary default
            F_eps=normcdf(VI+rho,0,sigma_eps);
            sigmaI=sigmaI*F_eps^sigIexp;
            
            % producer net worth etc
            Gfun = @(om1,om2,tau,zetaDWL)( (1-zetaDWL*params.outputDWL)*(1-tau)*om1*Y_om - om2*(1-tau)*wagebill ...
                            + om2*( (1-zetaDWL)*(1-(1-tau)*deltaK)*p - (1-tau)*pibar)*KB_g );                            
            N = Gfun(omplus,OmA,tauPi,0) - OmA*(1-tauPi*(1-theta)+delta*qB)*AB_g  + (1-OmA)*n0;
            %e = p*KB_g/N - qB*AB_g/N - (1-phi0);
            if eqfct_intern
                eIcost_intern = sigmaI/2*eI^2;
                eIcost_extern = 0;
                cI = 1 - sigmaI*eI;
                cI_norm = 1/cI;
                equcost_intern = phi1/2*equ^2;
                equcost_extern = 0; 
                equnorm = 1/(1-phi1*equ);
            else
                eIcost_intern = 0;
                eIcost_extern = sigmaI/2*eI^2;
                cI = 1 + sigmaI*eI;
                cI_norm = cI;                
                equcost_intern = 0;
                equcost_extern = phi1/2*equ^2;  
                equnorm = 1 + phi1*equ;
            end            
            Y_om_g = Y_om/N;
            Nminus = Gfun(omminus,1-OmA,0,zeta);
            DivP = N*(phi0 - equ - equcost_extern) - (1-OmA)*n0;
            nstar = (Gfun(omstar,1,tauPi,0) - (1-tauPi*(1-theta)+delta*qB)*AB_g )/N;
            
            % payoff to intermediary
            M= OmA*(1-tauPi_int) + Nminus/(AB_g/G);
            % marginal probability of default for borrowers
            % probability of default for intermediaries
			AS_g = fracS * AB_g;
			AStarget = fracStarget * AB_g;
            AI_g=AB_g - AS_g;
			
			if params.mktConsI
				BI_g=-xi*qB*AI_g/q;
			else
				BI_g = -xI*qB*AI_g;
            end
			


            
			WI=(AI_g*(M + qB*(delta*OmA))+BI_g)/G; % 
            DebtB=AI_g*(M + qB*(delta*OmA))/G;
            betaI_g=betaB_g*F_eps*(phiI/cI_norm + (1-phiI));
            F_eps_minus=-sigma_eps*normpdf((VI+rho)/sigma_eps);
            F_eps_plus=sigma_eps*normpdf((VI+rho)/sigma_eps);
            DI_eff=phiI*WI - eI - eIcost_extern - F_eps_minus-(1-F_eps)*WI;
            
			if params.mktConsI
				lamRplus = 1 - (betaI_g - tauPi*shieldI*rD)/q;
			else
				lamRplus = q - (betaI_g - tauPi*shieldI*rD);
			end
			
            PsiS = chi1adj/2 * ( AS_g / AStarget - 1)^2 * AStarget;
            PsiP = N*phi1/2*equ^2;
            BR_loss = zeta*(1-OmA)*(1-deltaK)*p*KB_g/G + zeta*params.outputDWL*(1-OmK)*mu_ZA*(KB_g/G)^(1-alpha)*L^alpha ...
                      + etaI/eta*bankDWL*zetaI*(1-F_eps)*DebtB + PsiS + pibar*KB_g + PsiP + sigmaI/2*eI^2;
            DWL = eta*BR_loss;
            BR_expense = Pbar*(1-eta)*BR_loss;
            		            
            
            % output, investment and borrower consumption
            X=(G-1 + deltaK)*KB_g/G;
            Tvec=Pbar.*params.T*Y;            
            CB = DivP + DI_eff + X*(p-1) + (1-tauB)*wB*Lbar(1) + Tvec(1) + BR_expense(1);
            LB = qB*AB_g/(p*KB_g);
                       
            % first-order conditions for borrowers
            % for AB
            scriptF = nstar/Y_om_g;
            betaP_g = betaB_g * ( phi0/equnorm + 1-phi0 );
            fx(1)= qB  - lamBplus*F - betaP_g*( OmA*(1-tauK_int+delta*qB) + f_om*scriptF );
            % for KB
            MPK=(1-alpha)*mu_om*mu_ZA*((KB_g/G)/L)^(-alpha); 
            fx(2)= p - lamBplus*Phi*p - betaP_g*( OmA*(p*(1-(1-tauK)*deltaK) - (1-tauK)*pibar) + (1-tauK)*OmK*MPK/mu_om ...
                                 - f_om*scriptF*pibar  ...
                                 + f_om*omstar*MPK*scriptF/mu_om ); 
            % for wages
            gammavec=[gammaB,gammaS];
            for i=1:2
                MPL= alpha*gammavec(i)*mu_om*mu_ZA*L/Lbar(i)*((KB_g/G)/L)^(1-alpha);
                fx(2+i)= (1-tauK)*OmK*MPL/mu_om - (1-tauK)*OmA*wagevec(i) - f_om*scriptF*(wagevec(i)-omstar*MPL); 
            end
            % producer budget
            fx(5) = equ - (p*KB_g/N - qB*AB_g/N - (1-phi0) + equcost_intern);
            
            % intermediary FOC
            fx(6)=qB*(1-xi*lamRplus) - betaI_g*(M + delta*OmA*qB);
            % intermediary budget
            fx(7)=(1-phiI)*WI + eI - eIcost_intern -( qB*AI_g + (q+shieldI*tauPi*rD)*BI_g);
            fx(8)=VI - (phiI*WI - eI - eIcost_extern - betaB_g*(F_eps_minus + (1-F_eps)*rho))/(1-betaB_g*F_eps);
            % borrower constraint (after-tax depreciation)
            fx(9)=Phi*p*KB_g - F*AB_g - lamBminus;
            
            %  debt holdings of savers
			if chiLinear == 0
				fx(10)=qB + chi1adj*(AS_g/AStarget-1) - betaS_g * (M + delta*OmA*qB);
				muS = -AS_g^(1/3);
			else
				muS = ( -(qB + chi1adj*(AS_g/AStarget-1) - betaS_g * (M + delta*OmA*qB) - chiLinear) )^(1/3);
            end
                        
            if ~take_Lscale_as_given
                fx(LscaleIdx)= Y-1;
            end

            
            
            % saver debt, consumption
            BS_g=BG_g-BI_g;
			WS=(AS_g*(M + qB*(delta*OmA))+BS_g)/G;
            CS = (1-tauS)*wS*Lbar(2) +WS - q*BS_g -qB*AS_g - PsiS -tauD*rD*BS_g + Tvec(2) + BR_expense(2) ;
                                   
            
            % compute other stats
            Lrate=1- (M+delta*qB*OmA)/(1+delta*qB);
            rB=log(1/qB+delta);
            mdebt=qB*AB_g;
            
            % check that goods market adds up
            Ltax=tauB*wB*Lbar(1)+tauS*wS*Lbar(2);
            YB=OmK*mu_ZA*(KB_g/G)^(1-alpha)*L^alpha;
            Corptax=tauK*YB  + (tauPi_int-tauK_int)*OmA*AB_g/G  + BI_g*(shieldI*tauPi*rD-kappa) ...
                   - OmA*tauK*(deltaK*p + pibar)*KB_g/G  - tauK*OmA*wagebill + tauD*rD*BS_g;
            T=Ltax + Corptax;
            bailout=F_eps_plus-(1-F_eps)*(WI - bankDWL*zetaI*DebtB);
            Gsp=gammaG*Y;
            GO=T-(1/G-q)*BG_g-sum(Tvec);
            Ycheck = CB + CS + GO + X + DWL - bailout;
            KBcheck = (1-deltaK)*KB_g/G  + X ;
            Gdef = -(T-Gsp-sum(Tvec));
            Xadj = X + DWL;
            
            % BGss that solves the GBC in steady-state
            tau_eff=tau(BG_g/Y);
            public_goods = gammaG * Y;
            spendingss = public_goods + sum(Tvec) + bailout;
            taxesss = (tau_eff-tauK*OmA)*wagebill + tauK*YB + OmA*AB_g*(tauPi_int-tauK_int) + BI_g*(shieldI*tauPi*rD-kappa) + BS_g*tauD*rD - tauK*OmA*(deltaK*p+pibar)*KB_g;
            surplus = taxesss - spendingss; 
            BGss = surplus/(1-q*mu_G);
            
            % return on risk-taker wealth
            row = (AI_g*(M+delta*OmA*qB) + BI_g ) ./ (AI_g*qB + BI_g*(q+shieldI*tauPi*rD)) - 1;            
            
            % state variables (other than KB and AB)
            WIsh=WI/((M+delta*OmA*qB)*AI_g+BG_g);
            BGsh=BG_g/((M+delta*OmA*qB)*AB_g+BG_g);
            WSsh=1-WIsh-BGsh;
            
            WI_end = (1-phiI)*WI+eI - eIcost_intern; 
            
            if print
                % print steady state values
                disp(' ');
                disp('Analytic steady state');
                disp('--- Rates ---');
                disp(['rB: ',num2str(rB)]);
                disp(['rD: ',num2str(rD)]);
                disp('--- Capital and Output ---');
                disp(['KB: ',num2str(KB_g)]);
                disp(['Y: ',num2str(Y)]);                
                disp(['X: ',num2str(X)]);                
                disp(['Xadj: ',num2str(Xadj)]);                
                disp(['p: ',num2str(p)]);
                disp(['Lscale: ',num2str(Lscale)]);
                disp(['Labor share: ',num2str(wagebill)]);
                disp(['Ycheck: ',num2str(Ycheck)]);                
                disp(['KBcheck: ',num2str(KBcheck)]);                
                disp('--- Debt ---');
                disp(['AB: ',num2str(AB_g)]);
                disp(['mark d: ',num2str(mdebt)]);
                disp(['BS: ',num2str(BS_g)]);
                disp(['BI: ',num2str(BI_g)]);
                disp(['BI/Y: ',num2str(-BI_g/Y)]);
                disp(['BI/K: ',num2str(-BI_g/KB_g)]);
                disp('--- Capital Risk ---');
                disp(['omstar: ',num2str(omstar)]);
                disp(['OmA: ',num2str(OmA)]);
                disp(['OmK: ',num2str(OmK)]);
                disp(['M: ',num2str(M)]);
                disp(['DWL: ',num2str(DWL)]);
                disp(['ROW: ',num2str(row)]);
                disp('--- Bond prices ---');
                disp(['qB: ',num2str(qB)]);
                disp(['q: ',num2str(q)]);
                disp('--- Producers ---');
                disp(['equ: ',num2str(equ)]);
                disp(['N: ',num2str(N)]);
                disp(['pibar*K: ',num2str(pibar*KB_g)]);
                disp('--- Borrower Quantities ---');
                disp(['CB: ',num2str(CB)]);
                disp(['LB: ',num2str(LB)]);
                disp(['WB: ',num2str(p*KB_g-qB*AB_g)]);
                disp(['DP: ',num2str(DivP)]);
                disp(['lamB: ',num2str(lamBplus)]);
                disp(['WI: ',num2str(WI)]);                
                disp('--- Saver Quantities ---');
                disp(['CS: ',num2str(CS)]);
                disp(['WS: ',num2str(WS)]); 
				disp(['AS_g/AB_g: ',num2str(AS_g/AB_g)]);
                disp(['WSsh: ',num2str(WSsh)]);                
                disp('--- Intermediary ---');
                disp(['Drate: ',num2str(1-F_eps)]);
                disp(['WI: ',num2str(WI)]);                
                disp(['VI: ',num2str(VI)]);                
                disp(['dI/WI: ',num2str((phiI*WI-eI-sigmaI/2*eI^2)/WI_end)]);                
                disp(['eI/WI: ',num2str(eI/WI_end)]);                
                disp(['DI_eff/WI: ',num2str(DI_eff/WI_end)]);                
                disp(['new equ: ',num2str((1-F_eps)*WI)]);                
                disp(['lamI: ',num2str(lamRplus)]);
                disp(['bailout: ',num2str(bailout)]);                
                disp('--- Government ---');
                disp(['BG: ',num2str(BG_g)]);                
                disp(['BGsh: ',num2str(BGsh)]);
                disp(['surplus: ',num2str(surplus)]);
                disp(['BGss: ',num2str(BGss)]);
                disp(['T: ',num2str(T)]);
                disp(['Labor Taxes: ',num2str(Ltax)]);
                disp(['Corp Taxes: ',num2str(Corptax)]);
                disp(['Gsp: ',num2str(Gsp)]);
                disp(['GO: ',num2str(GO)]);
                disp(['Gdef: ',num2str(Gdef)]);
            end
                      
            Sol=struct('qB',qB,...
                'X',X,...
                'q',q,...
                'equ',equ*N,...
                'cI',cI,...
                'cS',CS,...
                'wB',wB,...
                'wS',wS,... 
				'AI_g',AI_g,...
				'AS_g',AS_g,...
				'cB',CB,...
                'muR',-AI_g^(1/3),...
                'lamR',lamRplus^(1/3),...            
                'lamS',-BS_g^(1/3),...
				'muS',muS,...
                'lamB',lamBplus^(1/3));
			if params.sigmaI==0
				Sol.cI=exp(eI);
			end

            rhoS=1-1/nuS;  
            if nuS==1 || 1-betaS*exp(rhoS*g)<0                
                V=struct('CB',CB,...
                'CS',CS,...
                'CI',cI,...
                'VB',CB * exp(betaB * g / (1 - betaB) ),...
                'VS',CS * exp(betaS * g / (1 - betaS) ),...
                'VI',WI,...
                'wagebill',wagebill);                  
            else              
                V=struct('CB',CB,...
                'CS',CS,...
                'CI',cI,...
                'VB',CB * exp(betaB * g / (1 - betaB) ),...
                'VS',CS * ((1-betaS)/(1-betaS*exp(rhoS*g)))^(1/rhoS) ,...
                'VI',WI,...
                'wagebill',wagebill);                  
            end
 
            
            
            Add=struct('OmA',OmA,...
                         'OmK',OmK,...
                         'MP',M,...
                         'omstar',omstar,...
                         'fom',f_om,....
                         'omplus',omplus,...
                         'F_eps',F_eps,...
                         'F_eps_minus',F_eps_minus,...
                         'F_eps_plus',F_eps_plus,...
                         'bI',BI_g,...
                         'VI',WI,...
                         'bS',BS_g,...
                         'BG_g',BG_g,...
                         'AB',AB_g/G,...
                         'KB',KB_g/G,... 
                         'WI',WI,...
                         'L',L,...
                         'AB_g',AB_g,...
                         'KB_g',KB_g,...
                         'sig2B', sig2_om,...
                         'sig2B_zeta_xi_next',[sig2_om,zeta,xi],...
                         'Y_om',Y_om,...
                         'YB',Y,...                         
                         'zeta',zeta,...
                         'xi',xi,...
						 'AStarget',AStarget,...
                         'tauN',tauB,...
                         'ZA',1,...
						 'ZA_factor',1,...
                         'cB',CB,...
                         'N',N,...
                         'nstar',nstar,...
                         'DivP',DivP,...
						 'netResourcesB',CB);
                                
            State=struct('LB',qB*AB_g/p*KB_g,...
                         'KB',KB_g/G,...   
                         'WI',WI,...
                         'BG',BG_g/G,...
                         'Y',Y,...
                         'G',G);
                       
            statsout=struct('rB',rB,...
                'rD',rD,...
                'LTV',qB*AB_g/(p*KB_g),...
                'Lrate',Lrate,...
                'Lscale',Lscale);
            
            stvals=struct('Sol',Sol,...
                'V',V,...
                'Add',Add,...
                'State',State,...
                'statsout',statsout);                        
        end
                      
        
        function [solguessvec,Vguessvec,V_names]=assignGuess(stv)

           solguess=struct('qB',log(stv.Sol.qB),...
                'X',stv.Sol.X,...
                'q',log(stv.Sol.q),...
                'equ',stv.Sol.equ,...
                'cI',log(stv.Sol.cI),...
                'cS',log(stv.Sol.cS),...
                'wB',log(stv.Sol.wB),...
                'wS',log(stv.Sol.wS),...
				'AI_g',stv.Sol.AI_g,...
				'AS_g',stv.Sol.AS_g,...
				'cB',stv.Sol.cB,...
                'muR',stv.Sol.muR,...                
                'lamR',stv.Sol.lamR,...
                'lamS',stv.Sol.lamS,...
				'muS',stv.Sol.muS,...
                'lamB',stv.Sol.lamB);
            
            solguessvec=model.DSGEModel.structToVec(solguess);
                        
            Vguess=struct('UB',stv.V.CB,...
                'US',stv.V.CS,...
                'UI',stv.V.CI,...
                'VB',stv.V.VB,...
                'VS',stv.V.VS,...
                'VI',stv.V.VI,...
                'X',stv.Sol.X,...
                'qB',stv.Sol.qB,...
                'wbill',stv.V.wagebill,...
                'lamB',-0.5);
            
            Vguessvec=model.DSGEModel.structToVec(Vguess);
            V_names=fieldnames(Vguess);
            
        end
                
        
        function [x,fx,exit,i]=tryOtherGuesses(fhand,gvec,options)
            % list of guesses
            
            %  Toggle Lagrange multipliers as follows:
            %
            %                              
%             muR=solvec(11);
%             lamI=solvec(12);
%             lamS=solvec(13);
% 			muS=solvec(14);

%             gindex={[12,14],[12,14],[12,15],[12,14]};
%             gvals={[0.6,0.5],[-0.4,0.5],[-0.4,0.3],[0.6,-0.5]};            
% 
            gindex={12,12,[12,14],14,[12,16]};
            gvals={0.6,-0.4,[-0.4,0.3],0.3,[0.5,0.7]};            
            
            for i=1:length(gindex)
                newguess=gvec;
                newguess(gindex{i})=gvals{i};
                [x,fx,exit]=fsolve(fhand,newguess,options);
%                if exit>=0
                if exit>0
                    break;
                end                                
            end
        end
                     
        
    end % of static methods
    
    
end