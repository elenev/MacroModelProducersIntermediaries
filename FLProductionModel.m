 classdef FLProductionModel < model.DSGEModel
    
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
        function obj=FLProductionModel(params,endogenv,exogenv,vfct,pfct,tfct,tempfct,basegrid,jac,varargin)
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
            bGamma=params.bGamma;
            Pbar=params.Pbar;
            gammaG=params.gammaG;
            Lscale=params.Lscale;
            BG_foreign=params.BG_foreign;
            tauK_int=(1-theta)*tauK;
            tauPi_int=(1-theta)*tauPi;
            outputDWL=params.outputDWL;
            pibar=params.pibar;
			phi0=params.phi0;
			phi1=params.phi1;
            n0=params.n0;
            chiLinear=params.chiLinear;
            eqfct_intern=params.eqfct_intern;
			chi1adj=params.chi1adj;
            FLnodef=params.FLnodef;
            
            % extract state variables
            exst=point(1);
            KB=point(2);
            LB=point(3);
            BG=point(4);
            ZA=obj.Exogenv.pts_perm(exst,1);
            sig2B=obj.Exogenv.pts_perm(exst,2);
            zeta=obj.Exogenv.pts_perm(exst,3);
            xi=obj.Exogenv.pts_perm(exst,4);
			AStarget=obj.Exogenv.pts_perm(exst,5);
            
            fullendogvars=[KB,LB,BG];

            
            % compute next period's states
			MITshock_WS = 0;
            if isempty(varargin)
                State_next=obj.evaluateTrans(point);
            else
                State_next=varargin{1};
				if nargin>6
					MITshock = varargin{3};
					MITshock_WS = MITshock;
				end
            end
                                    
            % extract solution variables
            qB=exp(solvec(1));
            X=solvec(2);
            q=exp(solvec(3));
%            cB=exp(solvec(4));
            equ=solvec(4);
            cS=exp(solvec(5));  
            wvec=exp(solvec(6:7)); 
            wvec=wvec(:);
			cB = exp(solvec(8));
            
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
            
            
            % payoff to savers
            MP= OmA + Nminus/AB;
            YB=OmK*ZA*KB^(1-alpha)*L^alpha;                       
                        
            % government debt and taxes
            tauN=tau(BG/Y);
            DebtB=(MP+delta*OmA*qB)*AB;
            WS=DebtB+BG-BG_foreign*BG+MITshock_WS; 
%            ZA_factor=ZA/mu_ZA;
            ZA_factor=min(1,ZA/mu_ZA);
            tau_eff=tauN * ZA_factor^bTau; % procyclical 
            Tvec=Pbar * T * Y * ZA_factor^bT; % countercyclical 
            public_goods = gammaG * Y * ZA_factor^bGamma; % countercyclical 
       

			% Rebate
            AS_g = AB_g;
			if chi1adj>0
				PsiSA = chi1adj/2 * (AS_g / AStarget - 1)^2 * AStarget;
			else
				PsiSA = 0;
			end
			PsiS = PsiSA + chiLinear * AS_g;
			PsiP = N*phi1/2*eq_N^2;
            BR_expense = Pbar * ((1-eta)*(zeta*(1-OmA)*(1-deltaK)*p*KB + zeta*outputDWL*(1-OmK)*ZA*KB^(1-alpha)*L^alpha) ...
                                    + PsiS + PsiP + pibar*KB);            
            
            % budget constraint borrower
            %cB=(1-tau_eff)*wvec(1)*Lbar_scale(1) + Tvec(1) + BR_expense(1) + p*X + dI_eff + DivP - X - Psi*KB;
            netResourcesB = (1-tau_eff)*wvec(1)*Lbar_scale(1) + Tvec(1) + BR_expense(1) + p*X + DivP - X - Psi*KB;
			
            % saver
            rD=1/q-1;
            bS=( (1-tau_eff)*wvec(2)*Lbar_scale(2) + WS + Tvec(2) + BR_expense(2) - cS - qB*AS_g - PsiS)/(q+tauD*rD);          
           
            % next period government debt
             
            spending = public_goods + sum(Tvec);
            taxes = (tau_eff-tauK*OmA)*wagebill + tauK*YB - OmA*AB*tauK_int + bS*tauD*rD - tauK*OmA*(deltaK*p + pibar)*KB;
            BG_g = (BG + spending - taxes)/q;            
            
            if mode>0
                % simulation, mode contains number of next period's state
                exst=mode;
                exnpt=size(obj.Exogenv.pts_all,1);
                ZAnext=obj.Exogenv.pts_all(exst,1);
                cind=obj.Exogenv.pts_all(exst,end);     
                sig2B_zeta_xi_next=obj.Exogenv.pts_all(exst,2); %:end-1
                KBnext=State_next(exst,1);
                LBnext=State_next(exnpt+exst,1);
                BGnext=State_next(2*exnpt+exst,1);
                nextst=[cind,KBnext,LBnext,BGnext]; 
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
                BGnext=State_next(2*nst+(1:nst),1);
                % matrix of next period states
                nextst=[cind,KBnext,LBnext,BGnext];
			end            
                               
			addvars=struct;
			addvars.OmA=OmA;
			addvars.OmK=OmK;
			addvars.MP=MP;
			addvars.omstar=omstar;
			addvars.fom=fom;
			addvars.omplus=omplus;
			addvars.bS=bS;
			addvars.BG_g=BG_g;
			addvars.AB=AB;
			addvars.KB=KB;
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
            outstr.Jvec=[omstar,fom,mu_om,OmK,BG/Y];
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
            theta=params.theta;
            tauK=params.tauK;
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
			phi1=params.phi1;
			chiLinear=params.chiLinear;
            eqfct_intern=params.eqfct_intern;
			chi1adj=params.chi1adj;
            
            
            % extract endogeous variables 
            qB=exp(solvec(1));
            X=solvec(2);
            q=exp(solvec(3));
            equ=solvec(4);
            cS=exp(solvec(5));  
            wvec=exp(solvec(6:7));
			cB = exp(solvec(8));
            lamS=solvec(9);
			muS=solvec(10);
            lamB=solvec(11);
            wvec=wvec(:);                

            % multiplier transformations
            lamSplus=max(0,lamS)^3;
            lamSminus=max(0,-lamS)^3;           
			muSplus=max(0,muS)^3;
            muSminus=max(0,-muS)^3; 
            lamBplus=max(0,lamB)^3;
            lamBminus=max(0,-lamB)^3;           
            
            % extract some other state-dependent variables
            envec=instr.addvars;
            KB=envec.KB;
            AB=envec.AB;
            OmK=envec.OmK;
            OmA=envec.OmA;
            fom=envec.fom;
            omstar=envec.omstar;
            AB_g=envec.AB_g;    
            KB_g=envec.KB_g;
            BG_g=envec.BG_g;
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
            
            % capital price
            p = 1 + psi*(X/KB-(mu_G-1+deltaK));    
            
            % corporate bond market
            AS_g = AB_g;
			
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
            exp_S_q = expStruct.exp_S_q;
			exp_S_qB = expStruct.exp_S_qB;
            CE_vec = expStruct.CE_vec;
            p_next = expStruct.p_next;
            qB_next = expStruct.qB_next;
             
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
            
            % borrower issuance cost
            if eqfct_intern
                equnorm=1/(1-phi1*equ/N);
            else
                equnorm=1+phi1*equ/N;
            end
            U_normP=U_norm(1)*equnorm;
                                   
            % borrower FOCs
            fx(1) = qB - F*lamBplus - exp_B_qB/U_normP;
            fx(2) = p - Phi*p*lamBplus - exp_B_qK/U_normP;
            MPLvec=[gammaB, gammaS]./Lbar_scale' * alpha*mu_om*ZA*L*(KB/L)^(1-alpha);
            Fscript=nstar/Y_om;
            for j=1:2
                fx(2+j)= (1-tauK)*OmK*MPLvec(j)/mu_om - (1-tauK)*OmA*wvec(j) - fom*Fscript*(wvec(j) - omstar*MPLvec(j));
            end                        
                                                  
            % saver FOC
            rD=1/q-1;
            fx(5)=q + tauD*rD -lamSplus - exp_S_q/U_norm(2);                        
			fx(6) = qB - muSplus + dPsiSA + chiLinear  - exp_S_qB/U_norm(2);
            % saver constraints
            fx(7)=bS - lamSminus;
			fx(8) = AS_g - muSminus;
            
            % risk free debt market clearing
            fx(9)=BG_g - bS - BG_g*BG_foreign;
            					
            % budget constraint for borrower
			fx(10) = cB - netResourcesB;
            
            % borrower hard leverage constraint
            fx(11)=Phi*p*KB_g - F*AB_g - lamBminus;

			
            % Transitions
%             nst=size(nextst,1);
            nst = obj.Exogenv.exnpt;
            KBtrans = KB_g * ones(nst,1)./ mu_G;
            BGtrans = BG_g * ones(nst,1)./ mu_G;
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
                dOmA = -fom;
                dfom = -fom*(omstar*mu_om - mu_om^2 + sig2B)/(omstar*sig2B); 
                dOmK = -omstar*fom;

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
					cS,wvec(1),wvec(2),cB,...
                lamSplus,lamSminus,muSplus,muSminus,lamBplus,lamBminus,...
                OmA,fom,OmK,tauN,dOmA,dTau,dfom,dOmK,...                
                exp_B_qB,exp_B_qK,exp_S_q,exp_S_qB}];

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
                expind=[1,3,5:8];
                nexpind=[2,4];
                
                J(:,expind) = tmp_Jacobian(1:end,expind) .* repmat( exp(solvec(expind))' ,nfx,1);
                                    
                % Investment and portfolios (not exponentiated)
                J(:,nexpind) = tmp_Jacobian(1:end,nexpind);                

                % Multipliers
                start_idx = 8;
                
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
                Vnext(3)=V_vec(1);
                Vnext(4)=V_vec(2);
                Vnext(5)=X;
                Vnext(6)=qB;
                Vnext(7)=Lbar_scale'*wvec;
                V{1}=Vnext;
                % state transition
                V{2}=[KBtrans; LBtrans; BGtrans]';
                
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
                SDFB = expStruct.Conditional.SDFB/U_norm(1);
                SDFS = expStruct.Conditional.SDFS/U_norm(2);
                CB_next = expStruct.Conditional.CB_next;
                CS_next = expStruct.Conditional.CS_next;
                OmA_next = expStruct.Conditional.OmA_next;
                OmK_next = expStruct.Conditional.OmK_next;
                MP_next = expStruct.Conditional.MP_next;
                qB_next = expStruct.Conditional.qB_next;                
                fom_next = expStruct.Conditional.fom_next;
                omstar_next = expStruct.Conditional.omstar_next;
                muom_next = expStruct.Conditional.muom_next;
                muom_next_norm = reshape(muom_next,size(OmA_next));
                Fscript_next = expStruct.Conditional.Fscript_next;
               
                SDF.SDFB = SDFB;    
                SDF.SDFS= SDFS;
                
                % Define returns to intermediaries on corporate bonds
                retP = ( MP_next+delta*OmA_next.*qB_next ) / qB;
                expRP = prnext * retP;
                expOmA = prnext * OmA_next;
                stdRP = sqrt( prnext * ( retP - expRP).^2 );
                
                % Decomposition of risk premia
                % Note: if lamI = 0, then effect_Rf = 1 and effect_collatP = 0
                % and the standard asset pricing equation R-Rf =
                % -Rf*Cov[SDF,R] is recovered
                Rf=1./q;
                                
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
                CSgr =  (CS_next ./ cS) .* mu_G;
                expCBgr=prnext * CBgr;
                expCSgr=prnext * CSgr;
                
                stdCBgr=sqrt( prnext * (CBgr - expCBgr).^2 );
                stdCSgr=sqrt( prnext * (CSgr - expCSgr).^2 );
                               
                    
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
                              'expCSgr',expCSgr, ...
                              'stdCBgr',stdCBgr, ...
                              'stdCSgr',stdCSgr, ... % EER and SR
                              'expEROWB',expEROWB, ...
                              'SRB',SRB);                        
                
                Wtrans.KB = KBtrans;                          
                Wtrans.LB = LBtrans;
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
            theta=params.theta;
            tauK=params.tauK;
            tauPi=params.tauPi;
            alpha=params.alpha;
            deltaK=params.deltaK;
            delta=params.delta;
            psi=params.psi;
            mu_om=params.mu_om;
            try
               nuinv=1./[params.nuB,params.nuS];
            catch
               nuinv=1/params.nu*ones(2,1);
            end
            tauK_int=(1-theta)*tauK;
            gammaB=params.gammaB;
            gammaS=params.gammaS;
            Lbar_scale=params.Lbar*params.Lscale;
            outputDWL=params.outputDWL;
            pibar=params.pibar;
			phi0=params.phi0;
			phi1=params.phi1;
            n0=params.n0;
            eqfct_intern=params.eqfct_intern;
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
            VB_next=Pol_next(:,3);
            VS_next=Pol_next(:,4);
            X_next=Pol_next(:,5);
            qB_next=Pol_next(:,6);
            wagebill_next=Pol_next(:,7);
                                  
            
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
            
            % payoff to savers
            MP_next= OmA_next + Nminus./AB_next;                                    
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
                                                              
            % saver FOC
            SDFS=SDF_mat{2};
            exp_S_q = prnext*SDFS;
			FOCtmp=MP_next+delta*OmA_next.*qB_next;
			exp_S_qB = prnext*(SDFS.*(FOCtmp));
                       
            expStruct = struct;
			expStruct.exp_B_qB=exp_B_qB;
			expStruct.exp_B_qK=exp_B_qK;
			expStruct.exp_S_q=exp_S_q;
			expStruct.exp_S_qB=exp_S_qB;
			expStruct.CE_vec=CE_vec;
			expStruct.p_next=p_next;
			expStruct.qB_next=qB_next;
               
            Conditional = struct;
			Conditional.SDFB=SDFB;
			Conditional.SDFS=SDFS;
			Conditional.CB_next=CB_next;
			Conditional.CS_next=CS_next;
			Conditional.OmA_next=OmA_next;
			Conditional.VS_next=VS_next;
			Conditional.OmK_next=OmK_next;
			Conditional.MP_next=MP_next;
			Conditional.qB_next=qB_next;
			Conditional.fom_next=fom_next;
			Conditional.omstar_next=omstar_next;
			Conditional.Fscript_next=Fscript_next;
			Conditional.muom_next=muom_next;
               
            expStruct.Conditional = Conditional;           
        end           
        
        
        function [errmat,solmat,condmat,Wshtrans,SDFmat]=calcEEError(obj,pointmat)
            % function to compute Euler equation error at points in state
            % space given by pointmat
            nst=size(obj.Exogenv.pts_all,1);
            
            errmat=zeros(size(pointmat,1),obj.Pfct.Nof);
            solmat=zeros(size(errmat));
            condmat=zeros(size(pointmat,1),obj.NCOND);
            SDFmat=zeros(size(pointmat,1),nst);
            Wshtrans=zeros(size(pointmat,1),3*nst);
            
            evaluatePol = @(point)obj.evaluatePol(point);
            calcStateTransition = @(point,soltmp)obj.calcStateTransition(point,soltmp,0);
            calcEquations = @(exst,nextst,soltmp,outstr)obj.calcEquations(exst,nextst,soltmp,outstr,2);    
            % Should be parfor. Use for when debugging only
             parfor i=1:size(errmat,1)
%            for i=1:size(errmat,1)
                point=pointmat(i,:);
                soltmp=evaluatePol(point);
                % transition
                [nextst,outstr]=calcStateTransition(point,soltmp);
                % equations
                [fx,~,V]=calcEquations(point(1),nextst,soltmp,outstr);                                
                qB=exp(soltmp(1));
                p=exp(soltmp(2));
                q=exp(soltmp(3));
                wvec=exp(soltmp(6:7));
				cB=exp(soltmp(8));
                AB_g=outstr.addvars.AB_g;
                normvec=[qB,p,wvec',q,qB,qB*AB_g,AB_g,qB*AB_g,cB,qB*AB_g];
                               
                condvars=V{1};
				KBtrans=V{2}.KB;
                LBtrans=V{2}.LB;                
				BGtrans=V{2}.BG;
                SDFS=V{3}.SDFS;
                errmat(i,:)=fx'./normvec;
                solmat(i,:)=soltmp';
                condmat(i,:)=model.DSGEModel.structToVec(condvars)';
                Wshtrans(i,:) = [KBtrans', LBtrans', BGtrans'];
                SDFmat(i,:) = SDFS';
            end
            
        end
        
        % simulate model; overwrite method from superclass 
         function [simseries,varnames,errmat,Wshtrans,SDFmat]=simulate(obj,NT,NTini,inistvec,simerror,shmat)
            if length(inistvec)~=obj.Vfct.SSGrid.Ndim
                error('inistvec must be vector of length SSGrid.Ndim');
            end
            
            NTtot=NT+NTini;
            simseries=zeros(NTtot,1+obj.NSTEX+obj.NSTEN+obj.NSOL+obj.NV+obj.NADD+obj.NZNS);
            
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
                val_range=3:4;
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
        

        
        function [simseries, varnames] = computeSimulationMoments(obj, simseries, varnames, varargin)

			if nargin>3
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
            loglist=model.HelperCollection.makeListFromNames(indexmap,{'qB','q','p','cS','wB','wS','cB'});
            multlist=model.HelperCollection.makeListFromNames(indexmap,{'lamS','muS','lamB'});    

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
%            AB_g=simseries(:,indexmap.get('AB_g'));
            Psi=obj.Params.psi/2*(XbyK-(exp(params.g)-1+deltaK)).^2;
            DivP = simseries(:,indexmap.get('DivP'));
			eq_N = simseries(:,indexmap.get('equ')) ./ simseries(:,indexmap.get('N'));
			dP_rate = params.phi0 - eq_N - params.eqfct_intern*params.phi1/2*eq_N.^2;
            
            
            % ---------------------------------------------------------------------
            % wealth distribution
            % ---------------------------------------------------------------------

            % Total corporate debt
            corpDebt = (MP+deltaB*OmA.*qB).*simseries(:,indexmap.get('AB'));
            WSm= corpDebt+BG-params.BG_foreign./G;
            WBm = p.* KB - corpDebt;
			WB = (1 - simseries(:,indexmap.get('LB'))) .* p .* KB;
            PD = WBm./DivP;

            % ---------------------------------------------------------------------
            % interest rates and returns
            % ---------------------------------------------------------------------
            rB=log(1./qB+deltaB);
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
			
            % bond excess return
            theta=params.theta;
            Lret=(MP(2:end) +deltaB*OmA(2:end).*qB(2:end))./qB(1:end-1)-rD(1:end-1)-1;
            Lret = log(1+Lret);
            
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
%             bI = simseries(:,indexmap.get('bI'));
% %             Imktlev=-q.*bI./mdebt;
% %             Ibklev=-bI./mdebt;
%             bIstart= WSm-BG;
%             bI = -bI;
%             mktNWI=qB.*AI_g-bI;
%             
%             %dbar = simseries(:,indexmap.get('dbar'));
%             if params.sigmaI>0
%                 cI=simseries(:,indexmap.get('cI'));
%                 if params.eqfct_intern
%                     eI=(1-cI)/params.sigmaI ;
%                     dIcost_intern = params.sigmaI/2*eI.^2;
%                     dIcost_extern = 0;
%                 else
%                     eI=(cI-1)/params.sigmaI ;
%                     dIcost_intern = 0;
%                     dIcost_extern = params.sigmaI/2*eI.^2;
%                 end
%             else
%                 eI=log(simseries(:,indexmap.get('cI')));
%                 cI= ones(size(simseries,1),1);                
%                 dIcost_intern = zeros(size(simseries,1),1);
%                 dIcost_extern = zeros(size(simseries,1),1);
%             end
%             if params.fixediv
%                 phiIWI=params.phiI*params.WIbar*ones(size(WI));
%             else
%                 phiIWI=params.phiI*WI;                
%             end
%             dI_cost = dIcost_intern + dIcost_extern;
%             F_eps=simseries(:,indexmap.get('F_eps'));
%             F_eps_minus=simseries(:,indexmap.get('F_eps_minus'));
%             F_eps_plus=simseries(:,indexmap.get('F_eps_plus'));
%             brupt=1-F_eps;
%             newequity=(1-F_eps).*WI;
%             dI_eff=phiIWI - eI - dIcost_extern - F_eps_minus - newequity;
%             bailout=F_eps_plus - (1-F_eps).*(WI - zeta.*corpDebt);
%             payrate = phiIWI - eI - dIcost_extern - F_eps_minus; 
%             WI_end = WI - phiIWI + eI - dIcost_intern;
%             dI_rate = payrate./WI_end;
%             eI_rate = eI./WI_end;
%             
%             % Return on intermediary wealth, after taxes
%             end_equity = ( AI_g(1:end-1) .* (MP(2:end) + deltaB*ZA(2:end).*qB(2:end)) - bI(1:end-1));
%             start_equity = ( AI_g(1:end-1) .* qB(1:end-1) - (q(1:end-1)+shieldI*tauPi*rD(1:end-1)).*bI(1:end-1));
% 
%             ROW = end_equity ./ start_equity - 1;
%             %ROW = log(ROW+1);
% 
%             % Limited liability: Set to -100% when risk-takers declare bankruptcy
%             ROW(start_equity < 0) = 0;
%             profit = end_equity - start_equity;
%             % Set profit_t = -equity_ {t-1} if risk taker declares bankruptcy at time t
%             ROA = Lret + rD(1:end-1);
%                         
% 
%             % intermediary constraint
%             bind_lamI = (simseries(:,indexmap.get('lamR')) > 0);
%         
%             
%             % expected excess return on intermediary equity
%             VI=simseries(:,indexmap.get('VI'));            
%             totalvalue=VI+q.*bI;
%             debtcost= 1./(q + shieldI*tauPi*rD -params.kappa) -1;
%             acctprof = (1-theta)*OmA.*AB - rD.*bI;
%             franch = VI./WI-1;
%             franch(WI<0)=0;
%             accROE = (params.F*AB(2:end)).*Lret./WI(2:end);
%             accROE(WI(2:end)<0)=0;
%             
%             expRIE= simseries(:,indexmap.get('retBE'));
%             if ~isempty(expRIE)
%                 expERIE= expRIE - 1./q;
%                 expRIE = expRIE - 1;
%                 wacc = VI./totalvalue .* expRIE + (q.*bI)./totalvalue .* debtcost;
%                 profbility = rB - Lrate - wacc;
%             else
%                 expERIE=nan(NT_sim,1);
%                 expRIE=nan(NT_sim,1);
%                 wacc = nan(NT_sim,1);
%                 profbility = nan(NT_sim,1);
%             end            
%             

            % EER on loans based on conditional expectation (only works if
            % compEEErr ran)
            expRP=simseries(:,indexmap.get('expRP'));
            if ~isempty(expRP)
                expERP=expRP - 1./q;
            else
                expERP=nan(NT_sim,1);
            end   

                        
            % ---------------------------------------------------------------------
            % Consumption and welfare
            % ---------------------------------------------------------------------
            
            DWL=params.eta*((1-deltaK)*zeta.*Drate.*p.*KB + zeta.*params.outputDWL.*(1-OmK).*ZA.*KB.^(1-params.alpha).*L.^params.alpha);
            
            DWL_byY = DWL ./ Y_om;
            
            
            cB=simseries(:,indexmap.get('cB'));
            cS=simseries(:,indexmap.get('cS'));
            
            C = cB + cS;
            
            VB=simseries(:,indexmap.get('VB'));
            VS=simseries(:,indexmap.get('VS'));

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
            Pitax = zeros(size(Ktax));
            Dtax = tauD*rD.*simseries(:,indexmap.get('bS'));         
            grosstax =  Ltax + Ktax + Pitax + Dtax;
            Ktax_ben = tauK*OmA.*((1-theta)*AB + params.deltaK*p.*KB);
            Pitax_ben = zeros(size(Ktax_ben));
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
            byVar_table = table(WBm,X,XDWL,KB,wagebill,assets,mdebt,bkdebt,DWL,...
                cB,cS,tax,Ltax,Ltax0,Ktax,Pitax,Dtax,Ktax_ben,...
                Pitax_ben,spend,Corptax,Gdef_prim,Gdef0_prim,BG_g,C,Y,T,gammaG);
            byVar_varnames=byVar_table.Properties.VariableNames;
            byVar_table = varfun(@(var)var./Y_denom,byVar_table);
            byVar_table.Properties.VariableNames=strcat(byVar_varnames,'_byY');       
                        
            fracProfligacy = ( byVar_table{:,end-4} < params.BG1); 
            fracAusterity = (byVar_table{:,end-4} > params.BG2);
            Y_om_byY = Y_om./Y_denom;
            YDWL = Y - DWL;
            
            gammaG_byY = gammaG./Y_denom;
            T_byY = T./Y_denom;
            Ltax_byY = Ltax./Y_denom;
            
            lgammaG_byY = log(gammaG_byY);
            lT_byY = log(T_byY);
            lLtax_byY = log(Ltax_byY);            
            
            % Growth rates
            grVar_table = table(mdebt,p,Y_om,Y,C,X,XDWL,DWL,KB,AB,mktNWB,...
                gammaG,spend,tax,cB,cS,mdebt,bkdebt,welfare,VB,VS);
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
            
%             % financial intermediary (I)
%             % book equity
%             Ibkasset = F.*AI_g ./ trend_factor;
%             Ibklbt = bI ./ trend_factor;
% 			Ibklbt_byY = Ibklbt ./ Y_denom;
%             Ibkeqt = Ibkasset-Ibklbt;
%             % market equity
%             Imasset = qB.*AI_g ./ trend_factor;
%             Imlbt = q.*(bI) ./ trend_factor;
%             Imeqt = Imasset-Imlbt;
%             % book leverage
%             Ibklvg = Ibklbt./Ibkasset;
%             % market leverage
%             Imlvg = Imlbt./Imasset;
			
		                        
            
            statevec=simseries(:,1);
            simseries=[simseries(:,2:end),... state vars
                                      WB, WBm, WSm, ... wealth distribution
                                      XDWL, X_byKB, XDWL_byKB , p, PD, MPK, MPL, ... production
                                      rD, rB, Lspr, rT, LsprT, Tspr,  ...
                                      mdebt, LTV, mktLTV, Drate, LGD, LGD_equity, Lrate, ...
                                      DWL_byY, DivP, dP_rate, bind_lamB, ... corporate debt
                                      welfare, MUBS, ... consumption and welfare
                                      expERP,...
                                      fracProfligacy, fracAusterity, Y, YDWL, Y_om_byY, C, gammaG,...
                                      lgammaG_byY, lT_byY, lLtax_byY, logY, logX, logC, logZA, logG];
            simseries = [simseries, byVar_table{:,:}];
            simseries = [simseries,...
            Bbkasset,Bbklbt,Bbklbt_byY,Bbkeqt,Bmasset,Bmlbt,Bmeqt,Bbklvg,Bmlvg,Bbkdte]; % firm equity and leverage

            simseries=[simseries(2:end,:),Kret,Eret,Lret];
            simseries = [simseries,grVar_table{:,:}];

            varnames_add=[{'WB', 'WBm','WSm',... wealth distribution
                      'XDWL','X_byKB','XDWL_byKB', 'p', 'PD', 'MPK', 'MPL',... production
                      'rD', 'rB', 'Lspr', 'rT', 'LsprT', 'Tspr',   ...
                      'mdebt', 'LTV', 'mktLTV', 'Drate', 'LGD', 'LGD_equity', 'Lrate', ...
                      'DWL_byY', 'DivP', 'dP_rate', 'bind_lamB', ... corporate debt
                      'welfare', 'MUBS', ... consumption and welfare
                      'expERP', ...
                      'fracProfligacy', 'fracAusterity', 'Y', 'YDWL', 'Y_om_byY','C','gammaG',...
                      'lgammaG_byY', 'lT_byY', 'lLtax_byY', 'logY', 'logX', 'logC', 'logZA', 'logG'},...
                      byVar_table.Properties.VariableNames,...
                      {'Bbkasset','Bbklbt','Bbklbt_byY','Bbkeqt','Bmasset','Bmlbt','Bmeqt','Bbklvg','Bmlvg','Bbkdte'},...
                      {'Kret','Eret','Lret'},...
                      grVar_table.Properties.VariableNames];
            varnames=[varnames(2:end), varnames_add];
            
        end
        
%         function [res, wealth,transfers] = solveCVEV( obj, point, x, LHS )
%            KB_new = exp(x(1));
%            LB_new = x(2);
%            WI_new = x(3);
%            
%            point_new = [point(1), KB_new, LB_new, WI_new, point(5)];
%            val = obj.evaluateVal( point_new );
%            
%            if nargout>1
%                pol = obj.evaluatePol( point_new );
%                p_new = ( 1 + obj.Params.psi*(pol(2)/KB_new-(obj.Params.mu_G-1+obj.Params.deltaK)) );
%                qB_new = exp(pol(1));
%                AB_new = LB_next * (p_new * KB_new) / qB_new;
% 			   %AB_new = -(WB_new - p_new * KB_new) / qB_new;
%                [~,outstr] = obj.calcStateTransition( point_new, pol, 0);
%                envec = outstr.addvars;
%                wagebill = obj.Params.Lbar'*obj.Params.Lscale*exp(pol(7:9));
% 
%                BR = envec.zeta * ( (1-envec.OmA)*p_new*KB_new*(1-obj.Params.deltaK) ...
%                       + obj.Params.outputDWL*(1-envec.OmK)*envec.Y_om );
% 
%                borrowerWealth_new = (envec.YB - envec.OmA * wagebill) * (1-obj.Params.tauK) ...
%                    + envec.OmA * ( p_new * KB_new * (1 - (1 - obj.Params.tauK) * obj.Params.deltaK ) ...
%                    - (1 - (1-obj.Params.theta)*obj.Params.tauK + obj.Params.delta * qB_new) * AB_new) ...
%                    + obj.Params.Pbar(1) * (1-obj.Params.eta) * BR;
% 
%                intermediaryWealth_new = WI_new ...
%                    + obj.Params.Pbar(2) * (1-obj.Params.eta) * BR;
% 
%                CapitalTaxes = obj.Params.tauK * (envec.YB - envec.OmA * wagebill) ...
%                    + envec.OmA * AB_new * (1 - obj.Params.theta) * (obj.Params.tauPi - obj.Params.tauK) ...
%                    - envec.OmA * obj.Params.tauK * obj.Params.deltaK * p_new * KB_new;
%                totalWealth_new = envec.Y_om - wagebill + p_new * KB_new * (1 - obj.Params.deltaK) + point(5) ...
%                    - CapitalTaxes - obj.Params.eta * BR;
%                
%                wealth = struct('totalWealth',totalWealth_new,...
%                    'borrowerWealth',borrowerWealth_new,...
%                    'intermediaryWealth',intermediaryWealth_new);
%                res = zeros(size(x));
%                
%                if nargout>2
%                   totalWealth = LHS(4);
%                   borrowerWealth = LHS(5);
%                   intermediaryWealth = LHS(6);
%                    
%                   TB = borrowerWealth_new - borrowerWealth; 
%                   TI = intermediaryWealth_new - intermediaryWealth;
%                   TS = totalWealth_new - totalWealth - (TB + TI);
%                   transfers = struct('TB',TB,'TS',TS,'TI',TI);
%                end
%            else
%                res = val(4:6) ./ LHS(1:3) - 1;             
%            end
%            % res = [VB, VS, VI, p*KB]
% 		end
        
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
            syms KB LB BG
            
            syms qB X q equ cS wvec1 wvec2 cB
            syms lamSplus lamSminus
            syms lamBplus lamBminus
			syms muSplus muSminus
            
            syms exp_B_qB exp_B_qK exp_S_q exp_S_qB
            
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
                    
            
            % payoff to savers
            MP= OmA + Nminus/AB;
            YB=OmK*ZA*KB^(1-alpha)*L^alpha;                                               
            DebtB=(MP+delta*OmA*qB)*AB; 
            
            YnetDWL = Y;%Y - DWL;%Y
            
            % taxes
            syms tau(D_byY)
            taxrule_arg = BG/YnetDWL;
            
			AS_g = AB_g;
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
                              + PsiS + PsiP + pibar*KB;
            BR_expense1 = Pbar(1) * BR_expense_all;
            BR_expense2 = Pbar(2) * BR_expense_all;
            
            
            % government debt and taxes
            tauN=tau(taxrule_arg); % + tau0;
            WS=DebtB+BG-BG_foreign*BG;
            tau_eff=tauN * ZA_factor^bTau;
            Tvec=Pbar .* T * YnetDWL * ZA_factor^bT;
            public_goods = gammaG * YnetDWL * ZA_factor^bGamma;
            
            % budget constraint borrower
            netResourcesB=(1-tau_eff)*wvec1*Lbar_scale(1) + Tvec(1) + BR_expense1 + p*X + DivP - X - Psi*KB;
                                 
            % saver
            rD=1/q-1;
            bS=( (1-tau_eff)*wvec2*Lbar_scale(2) + WS + Tvec(2) + BR_expense2 - cS - qB*AS_g - PsiS)/(q+tauD*rD);          
           
            % next period government debt             
            spending = public_goods + sum(Tvec);
            taxes = (tau_eff-tauK*OmA)*wagebill + tauK*YB - OmA*AB*tauK_int + bS*tauD*rD - tauK*OmA*(deltaK*p+pibar)*KB;
            BG_g = (BG + spending - taxes)/q;           
            
            % compute (sym) current U's to normalize back SDFs
            nuinv_vec=nuinv;            
            %cvec=[cB,cS,cI];
            UB=cB;
            U_vecB=UB;
            U_vecS=cS;
            U_normB=U_vecB^(1-nuinv_vec(1))/cB;  
            U_normS=U_vecS^(1-nuinv_vec(2))/cS;              
            U_normP=U_normB*equnorm;
                        
            
            % borrower FOCs
            fx1 = qB - F*lamBplus - exp_B_qB/U_normP;
            fx2 = p - Phi*p*lamBplus - exp_B_qK/U_normP;
            MPLvec1=gammaB/Lbar_scale(1)*alpha*mu_om*ZA*L*(KB/L)^(1-alpha);
            MPLvec2=gammaS/Lbar_scale(2)*alpha*mu_om*ZA*L*(KB/L)^(1-alpha);
            Fscript=nstar/Y_om;
            fx3=(1-tauK)*OmK*MPLvec1/mu_om - (1-tauK)*OmA*wvec1 - fom*Fscript*(wvec1 - omstar*MPLvec1);
            fx4=(1-tauK)*OmK*MPLvec2/mu_om - (1-tauK)*OmA*wvec2 - fom*Fscript*(wvec2 - omstar*MPLvec2);
                                                  
            % saver FOC
            fx5=q + tauD*rD -lamSplus - exp_S_q/U_normS;
			fx6=qB - muSplus + dPsiSA + chiLinear - exp_S_qB/U_normS;
            % saver constraints
            fx7=bS - lamSminus;
			fx8=AS_g - muSminus;
                                    
            % risk free debt market clearing
            fx9=BG_g - bS - BG_g*BG_foreign;
            			
            % budget constraint for borrower
			fx10=cB - netResourcesB;
			
			% borrower hard constraint
            fx11=Phi*p*KB_g - F*AB_g - lamBminus;
			
            
            fx=[fx1;fx2;fx3;fx4;fx5;fx6;fx7;fx8;fx9;fx10;fx11];
            
            Jacobian = jacobian(fx,[qB,X,q,equ,cS,wvec1,wvec2,cB,...
                lamSplus,lamSminus,...
				muSplus,muSminus,...
                lamBplus,lamBminus]);

            syms OmA fom OmK Tau dOmA dOmK dTau dfom dOmK 
            sym_lookup_fOmA = subs(fOmA(omega),omega,omstar);
			sym_lookup_ffom = subs(ffom(omega),omega,omstar);
			sym_lookup_fOmK = subs(fOmK(omega),omega,omstar);
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
			
			Jacobian=subs(Jacobian,sym_lookup_dOmA,'dOmA');
			Jacobian=subs(Jacobian,sym_lookup_dTau,'dTau');
			Jacobian=subs(Jacobian,sym_lookup_dfom,'dfom'); 
			Jacobian=subs(Jacobian,sym_lookup_dOmK,'dOmK'); 
			Jacobian=subs(Jacobian,sym_lookup_Tau,'Tau');
			Jacobian=subs(Jacobian,sym_lookup_fOmA,'OmA');
			Jacobian=subs(Jacobian,sym_lookup_ffom,'fom'); 
			Jacobian=subs(Jacobian,sym_lookup_fOmK,'OmK');  
			
			fx=subs(fx,sym_lookup_dOmA,'dOmA');
			fx=subs(fx,sym_lookup_dTau,'dTau');
			fx=subs(fx,sym_lookup_dfom,'dfom'); 
			fx=subs(fx,sym_lookup_dOmK,'dOmK'); 
			fx=subs(fx,sym_lookup_Tau,'Tau');
			fx=subs(fx,sym_lookup_fOmA,'OmA');
			fx=subs(fx,sym_lookup_ffom,'fom'); 
			fx=subs(fx,sym_lookup_fOmK,'OmK');  
			
			%Jacobian = simplify(Jacobian);
           
            order=[ZA, sig2B, zeta, xi, AStarget, ZA_factor, KB, LB, BG,...
                qB, X, q, equ, cS, wvec1, wvec2,cB,...
                lamSplus,lamSminus,...
				muSplus,muSminus,...
                lamBplus,lamBminus,...
                OmA, fom,OmK, Tau, dOmA, dTau, dfom, dOmK,...
                exp_B_qB, exp_B_qK, exp_S_q, exp_S_qB];
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
            pibar=params.pibar;
            phi0=params.phi0;
            phi1=params.phi1;
            n0=params.n0;
            eqfct_intern=params.eqfct_intern;
			chi1adj=params.chi1adj;
			fracStarget=params.fracStarget(omstate);
			FLnodef=params.FLnodef;
            
            
            tauK_int=(1-theta)*tauK;
            
                      
            % solution vars
            AB_g=exp(sol(1));
            KB_g=exp(sol(2));
            wB=exp(sol(3));
            wS=exp(sol(4));
            qB=exp(sol(5));
            lamB=sol(6);
            equ=sol(7);
            
			% 5-dim sol: no AS_g, Lscale pre-set
			% 6-dim sol: AS_g and Lscale pre-set OR no AS_g and Lscale
			% 7-dim sol: AS_g and Lscale
		
			if numel(sol)==8 
				Lscale=exp(sol(8));
				LscaleIdx = 8;
				take_Lscale_as_given=0;
			else
				Lscale=params.Lscale;
				take_Lscale_as_given=1;
			end
			
            
            % growth
            G=mu_G; %=exp(g)

            lamBplus=max(0,lamB)^3;
            lamBminus=max(0,-lamB)^3;            

            % short-term bond prices
            betaS_g=betaS*exp(-g/nuS);
            betaB_g=betaB*exp(-g/nuB);
            q=(tauD+betaS_g)/2 + sqrt( (tauD+betaS_g)^2/4 - tauD );
            rD=1/q-1;
            
            % capital price
            p=1;
            
            % government budget and debt
            BG_g=params.BG_guess*(1-params.BG_foreign); % just as a guess (does not matter for prices)
            tauB=tau(BG_g);
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
            
            % producer net worth etc
            Gfun = @(om1,om2,tau,zetaDWL)( (1-zetaDWL*params.outputDWL)*(1-tau)*om1*Y_om - om2*(1-tau)*wagebill ...
                            + om2*( (1-zetaDWL)*(1-(1-tau)*deltaK)*p - (1-tau)*pibar)*KB_g );                            
            N = Gfun(omplus,OmA,tauPi,0) - OmA*(1-tauPi*(1-theta)+delta*qB)*AB_g  + (1-OmA)*n0;
            %e = p*KB_g/N - qB*AB_g/N - (1-phi0);
            if eqfct_intern
                equcost_intern = phi1/2*equ^2;
                equcost_extern = 0; 
                equnorm = 1/(1-phi1*equ);
            else
                equcost_intern = 0;
                equcost_extern = phi1/2*equ^2;  
                equnorm = 1 + phi1*equ;
            end            
            Y_om_g = Y_om/N;
            Nminus = Gfun(omminus,1-OmA,0,zeta);
            DivP = N*(phi0 - equ - equcost_extern) - (1-OmA)*n0;
            nstar = (Gfun(omstar,1,tauPi,0) - (1-tauPi*(1-theta)+delta*qB)*AB_g )/N;
            
            % payoff to savers
            M= OmA + Nminus/(AB_g/G);
            % marginal probability of default for borrowers
            % probability of default for intermediaries
			AS_g = AB_g;
			AStarget = fracStarget * AB_g;
	           
			
            PsiS = chi1adj/2 * ( AS_g / AStarget - 1)^2 * AStarget;
            PsiP = N*phi1/2*equ^2;
            BR_loss = zeta*(1-OmA)*(1-deltaK)*p*KB_g/G + zeta*params.outputDWL*(1-OmK)*mu_ZA*(KB_g/G)^(1-alpha)*L^alpha ...
                      + PsiS + pibar*KB_g + PsiP;
            DWL = eta*BR_loss;
            BR_expense = Pbar*(1-eta)*BR_loss;
            		            
            
            % output, investment and borrower consumption
            X=(G-1 + deltaK)*KB_g/G;
            Tvec=Pbar.*params.T*Y;            
            CB = DivP + X*(p-1) + (1-tauB)*wB*Lbar(1) + Tvec(1) + BR_expense(1);
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
            % borrower constraint (after-tax depreciation)
            fx(6)=Phi*p*KB_g - F*AB_g - lamBminus;
            
            %  debt holdings of savers
			fx(7)=qB + chi1adj*(AS_g/AStarget-1) - betaS_g * (M + delta*OmA*qB);
			muS = -AS_g^(1/3);
                        
            if ~take_Lscale_as_given
                fx(LscaleIdx)= Y-1;
            end           
            
            % saver debt, consumption
            BS_g=BG_g;
			WS=(AS_g*(M + qB*(delta*OmA))+BS_g)/G;
            CS = (1-tauS)*wS*Lbar(2) +WS - q*BS_g -qB*AS_g - PsiS -tauD*rD*BS_g + Tvec(2) + BR_expense(2) ;
                                   
            
            % compute other stats
            Lrate=1- (M+delta*qB*OmA)/(1+delta*qB);
            rB=log(1/qB+delta);
            mdebt=qB*AB_g;
            
            % check that goods market adds up
            Ltax=tauB*wB*Lbar(1)+tauS*wS*Lbar(2);
            YB=OmK*mu_ZA*(KB_g/G)^(1-alpha)*L^alpha;
            Corptax=tauK*YB -tauK_int*OmA*AB_g/G  ...
                   - OmA*tauK*(deltaK*p + pibar)*KB_g/G  - tauK*OmA*wagebill + tauD*rD*BS_g;
            T=Ltax + Corptax;
            Gsp=gammaG*Y;
            GO=T-(1/G-q)*BG_g-sum(Tvec);
            Ycheck = CB + CS + GO + X + DWL;
            KBcheck = (1-deltaK)*KB_g/G  + X ;
            Gdef = -(T-Gsp-sum(Tvec));
            Xadj = X + DWL;
            
            % BGss that solves the GBC in steady-state
            tau_eff=tau(BG_g/Y);
            public_goods = gammaG * Y;
            spendingss = public_goods + sum(Tvec);
            taxesss = (tau_eff-tauK*OmA)*wagebill + tauK*YB - OmA*AB_g*tauK_int + BS_g*tauD*rD - tauK*OmA*(deltaK*p+pibar)*KB_g;
            surplus = taxesss - spendingss; 
            BGss = surplus/(1-q*mu_G);
            
            
            % state variables (other than KB and AB)
            BGsh=BG_g/((M+delta*OmA*qB)*AB_g+BG_g);
            WSsh=1-BGsh;
            
            
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
                disp('--- Capital Risk ---');
                disp(['omstar: ',num2str(omstar)]);
                disp(['OmA: ',num2str(OmA)]);
                disp(['OmK: ',num2str(OmK)]);
                disp(['M: ',num2str(M)]);
                disp(['DWL: ',num2str(DWL)]);
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
                disp('--- Saver Quantities ---');
                disp(['CS: ',num2str(CS)]);
                disp(['WS: ',num2str(WS)]); 
                disp(['WSsh: ',num2str(WSsh)]);                
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
                'cS',CS,...
                'wB',wB,...
                'wS',wS,... 
				'cB',CB,...
                'lamS',-BS_g^(1/3),...
				'muS',muS,...
                'lamB',0.3);

            rhoS=1-1/nuS;  
            if nuS==1 || 1-betaS*exp(rhoS*g)<0                
                V=struct('CB',CB,...
                'CS',CS,...
                'VB',CB * exp(betaB * g / (1 - betaB) ),...
                'VS',CS * exp(betaS * g / (1 - betaS) ),...
                'wagebill',wagebill);                  
            else              
                V=struct('CB',CB,...
                'CS',CS,...
                'VB',CB * exp(betaB * g / (1 - betaB) ),...
                'VS',CS * ((1-betaS)/(1-betaS*exp(rhoS*g)))^(1/rhoS) ,...
                'wagebill',wagebill);                  
            end
 
            
            
            Add=struct('OmA',OmA,...
                         'OmK',OmK,...
                         'MP',M,...
                         'omstar',omstar,...
                         'fom',f_om,....
                         'omplus',omplus,...
                         'bS',BS_g,...
                         'BG_g',BG_g,...
                         'AB',AB_g/G,...
                         'KB',KB_g/G,... 
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
                'cS',log(stv.Sol.cS),...
                'wB',log(stv.Sol.wB),...
                'wS',log(stv.Sol.wS),...
				'cB',stv.Sol.cB,...
                'lamS',stv.Sol.lamS,...
				'muS',stv.Sol.muS,...
                'lamB',stv.Sol.lamB);
            
            solguessvec=model.DSGEModel.structToVec(solguess);
                        
            Vguess=struct('UB',stv.V.CB,...
                'US',stv.V.CS,...
                'VB',stv.V.VB,...
                'VS',stv.V.VS,...
                'X',stv.Sol.X,...
                'qB',stv.Sol.qB,...
                'wbill',stv.V.wagebill);
            
            Vguessvec=model.DSGEModel.structToVec(Vguess);
            V_names=fieldnames(Vguess);
            
        end
                
        
        function [x,fx,exit,i]=tryOtherGuesses(fhand,gvec,options)
            % list of guesses
            
            %  Toggle Lagrange multipliers as follows:
            %
            %                              
%             lamS=solvec(9);
% 			muS=solvec(10);
%             lamB=solvec(11);

 
            gindex={11};
            gvals={0.6,-0.4};            
            
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