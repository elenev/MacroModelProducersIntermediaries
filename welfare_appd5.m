function welfare_appd5(numClust)

experdef = '20200910';
suffix = '_s130';

econ_prefixes = {'altshares','sigI8','sigI0','','benchriskysafe', ...
	'noshieldsamerate','benchriskierrisky','betaBup','betaSdown', 'nosavers'};

bench = [];

EVs=zeros(numel(econ_prefixes),5);
EVs=array2table(EVs);
EVs.Properties.VariableNames={'pctdiff_VB','pctdiff_VS','deltaWB','deltaWS','Total'};
for i=1:numel(econ_prefixes)
	if strcmp(econ_prefixes{i},'')
		bench_econ_name = 'xi91benchgrid';
		alt_econ_name = 'xi95benchgrid';
		econ = 'bench';
	else
		bench_econ_name = [econ_prefixes{i},'xi91'];
		alt_econ_name = [econ_prefixes{i},'xi95'];
		econ = econ_prefixes{i};
	end
	bench.res = load(['res_',experdef,'_',bench_econ_name,suffix,'.mat']);
	bench.sim = load(['sim_res_',experdef,'_',bench_econ_name,suffix,'.mat']);
	alt.res = load(['res_',experdef,'_',alt_econ_name,suffix,'.mat']);
	alt.sim = load(['sim_res_',experdef,'_',alt_econ_name,suffix,'.mat']);
	
	values=struct2cell(bench.res.mobj.Params);
	paramout=cell2table(values,'RowNames',fieldnames(bench.res.mobj.Params)); 

	% Cluster
	gv = @(x)bench.sim.indexmap.get(x);
	enstatesidx = [gv('KB'),gv('LB'),gv('WI'),gv('BG')];
	[clusters,probs,~,rsq] = clusterEachExogState( bench.sim.statevec(2:end), ...
		bench.sim.simseries( :, enstatesidx), numClust, false );

	%EXOG_STATE=5;
	%OFFSET=4;
	%state = [EXOG_STATE, mean(bench.pzns.simseries(:,OFFSET+1:OFFSET+4))];
	val_bench = bench.res.mobj.evaluateVal(clusters);
	VB_bench = val_bench(4,:)';
	VS_bench = val_bench(5,:)';
	pzns_bench = bench.res.mobj.Zfct.evaluateAt(clusters);
	qCB_bench = pzns_bench(2,:)';
	qCS_bench = pzns_bench(3,:)';
	
	fprintf('Processing %s\n',econ);

	% For Each alternative economy, compute VX_bench / VX_alternative
	% diff_VX = how much higher agent's welfare is in alternative vs. bench
	% diff_VX = by how many percent do you need to increase consumption in 
	% every state of the world in benchmark economy to make her as well 
	% of as she was in the alternative
	%alt_state =  [EXOG_STATE, mean(alt.sim.simseries(:,OFFSET+1:OFFSET+4))];
	gv = @(x)alt.sim.indexmap.get(x);
	enstatesidx = [gv('KB'),gv('LB'),gv('WI'),gv('BG')];
	[alt_clusters,alt_probs] = clusterEachExogState( alt.sim.statevec(2:end), ...
		alt.sim.simseries( :, enstatesidx), numClust, false );
	val_alt = alt.res.mobj.evaluateVal(alt_clusters);
	VB_alt = val_alt(4,:)';
	VS_alt = val_alt(5,:)';
	diff_VB = VB_alt' ./ VB_bench - 1;
	diff_VS = VS_alt' ./ VS_bench - 1;
	delta_WB = diff_VB .* qCB_bench;
	delta_WS = diff_VS .* qCS_bench;
	totalDelta = delta_WB + delta_WS;

	EVs{i,1} = probs' * diff_VB * alt_probs;
	EVs{i,2} = probs' * diff_VS * alt_probs;
	EVs{i,3} = probs' * delta_WB * alt_probs;
	EVs{i,4} = probs' * delta_WS * alt_probs;
	EVs{i,5} = probs' * totalDelta * alt_probs;
	EVs.Properties.RowNames{i} = econ;

end

%for j=1:4
% writetable(EVs,'Results/statsexog_res_20200124_bench_s130.xls','Sheet',j,'WriteRowNames',true);
%end
%writetable(paramout,'Results/statsexog_res_20200124_bench_s130.xls','Sheet','params','WriteRowNames',true);
writetable(EVs,['Results/welfareappd5_',experdef,'_s130.xls'], ...
	'WriteRowNames',true,'Sheet',sprintf('NumClust %d',numClust));
end

