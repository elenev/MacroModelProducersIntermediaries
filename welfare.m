function welfare(numClust)

experdef = '20200910';

results_files = dir(['sim_res_',experdef,'_*_s130.mat']);

pzns_res = ['res_',experdef,'_bench_s130.mat'];
pzns_simres = ['sim_res_',experdef,'_bench_s130.mat'];
%bench.pzns=load(['pzns_',pzns_res]);
bench.res=load(pzns_res);
bench.sim=load(pzns_simres);

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

bench = [];

EVs=zeros(size(results_files,1),5);
EVs=array2table(EVs);
EVs.Properties.VariableNames={'pctdiff_VB','pctdiff_VS','deltaWB','deltaWS','Total'};
for i=1:size(results_files,1)
   fname = results_files(i).name;
   fprintf('Processing %s last modified at %s\n',fname,results_files(i).date);
   parts = strread(fname,'%s','delimiter','_');
   alt.sim=load(fname,'simseries','errmat','indexmap','statevec');
   alt.res=load(fname(5:end),'mobj');
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
   EVs.Properties.RowNames(i) = parts(4);
   
end

%for j=1:4
%    writetable(EVs,'Results/statsexog_res_20200124_bench_s130.xls','Sheet',j,'WriteRowNames',true);
%end
%writetable(paramout,'Results/statsexog_res_20200124_bench_s130.xls','Sheet','params','WriteRowNames',true);
writetable(EVs,['Results/welfare_',experdef,'_bench_s130.xls'], ...
	'WriteRowNames',true,'Sheet',sprintf('NumClust %d',numClust));
end

