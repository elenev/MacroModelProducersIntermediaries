function [res,simnext]  = computeMITShockState( point, x, mobj, MITshock)

NSOL=mobj.NSOL;
exnpt=mobj.Exogenv.exnpt;

solvec = x(1:NSOL);
KBtrans = x(NSOL+1);
LBtrans = x(NSOL+1+(1:exnpt));
WItrans = x(NSOL+1+exnpt+(1:exnpt));
BGtrans = x(NSOL+1+2*exnpt+1);
thisVI = x(NSOL+1+2*exnpt+2);
%point = [avgstate,KBpt,WBpt,WIpt,BGpt];

%trans=[(1:exnpt)',KBtrans*ones(exnpt,1),LBtrans,WItrans,BGtrans*ones(exnpt,1)];
trans=[KBtrans*ones(exnpt,1);LBtrans;WItrans;BGtrans*ones(exnpt,1)];
%trans=reshape(trans,exnpt,4);

[nextst,outstr]=mobj.calcStateTransition(point,solvec,0,trans,thisVI,MITshock);
[fx,~,V]=mobj.calcEquations(point(1),nextst,solvec,outstr,1);
KBtrans_check = V{2}(1);
LBtrans_check = V{2}(exnpt+(1:exnpt))';
WItrans_check = V{2}(2*exnpt+(1:exnpt))';
BGtrans_check = V{2}(3*exnpt+1);

res = zeros(size(x));
res(1:NSOL) = fx;
res(NSOL+1) = KBtrans - KBtrans_check;
res(NSOL+1+(1:exnpt)) = LBtrans - LBtrans_check;
res(NSOL+1+exnpt+(1:exnpt)) = WItrans - WItrans_check;
res(NSOL+1+2*exnpt+1) = BGtrans - BGtrans_check;
res(NSOL+1+2*exnpt+2) = thisVI - V{1}(6);

if nargout>1
	outstr.addvars.sig2B_zeta_xi_next = outstr.addvars.sig2B_zeta_xi_next(5,2);
	addvec=model.DSGEModel.structToVec(outstr.addvars)';
	valvec=mobj.evaluateVal(point)'; 
	znsvec=[];
    if mobj.NZNS>0
	 	znsvec=mobj.Zfct.evaluateAt(point)';
    end
	% write different categories of variables in one row
	simnext=[point(1),outstr.exstvec',point(2:end),solvec',valvec,addvec,znsvec];
end

end