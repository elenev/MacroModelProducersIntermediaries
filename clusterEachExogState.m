function [clusters,probs,cindex,rsq] = clusterEachExogState( statevec, simenstates, maxclust, diagnostics )
	exogStates = unique(statevec);
	exnpt = length(exogStates);
	totclust = maxclust * exnpt;
	
	cindex = zeros( size( statevec) );
	clusters = zeros( totclust, 1 + size(simenstates,2));
	probs = zeros( totclust, 1);
	
	for ii = 1:exnpt
		exogIdx = statevec==exogStates(ii);
		enstatemat=simenstates(exogIdx, :);
		cindex( exogIdx ) = (ii-1)*maxclust + ...
			clusterdata(enstatemat,'criterion','distance','maxclust',maxclust,'linkage','weighted');
		clusters( (ii-1)*maxclust + (1:maxclust), 1) = exogStates(ii);
	end
	
	for jj=1:totclust
		clusters(jj,2:end) = mean( simenstates( cindex == jj, : ), 1 );
		probs(jj) = mean( cindex == jj );
	end
	
	rsq = zeros( 1, size(simenstates,2) );
	if ~isempty(diagnostics) && diagnostics
		tab=array2table(simenstates);
		tab.statevec=categorical(statevec);
		tab.cluster=categorical(cindex);
		rsq = zeros( 1, size(simenstates,2) );
		for kk = 1:size(simenstates,2)
			mdl = fitlm(tab,sprintf('simenstates%d ~ cluster',kk));
			rsq(kk) = mdl.Rsquared.Ordinary;
		end
	end
end