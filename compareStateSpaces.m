% Plot Appendix Figure B.2

reload = 1;
if reload
    run('experdef_20200910.m');
    bench=load('sim_res_20200910_bench_s130.mat');
    benchgrid = benchgrid.finegrid;
    fields = fieldnames(benchgrid);
    %compare=load('sim_res_20180511_betaB94_s130.mat');
    %comparegrid = betaBgrid.finegrid;
end

%legend_entry='betaB=0.94';

idxarray = { bench.indexmap.get('KB'), ...
    bench.indexmap.get('LB'), ...
    bench.indexmap.get('WI'), ...
    bench.indexmap.get('BG') };

captions = {'K^B', 'L^P', 'N^I', 'B^G'}; 

f=figure;
for i=1:4
    idx = idxarray{i};
    h = subplot(2,2,i);
    colorList = get(h,'ColorOrder');
    benchhist=histogram(bench.simseries(:,idx));
	benchhist.LineWidth = 0.25;
    %hold on;
    %comparehist=histogram(compare.simseries(:,idx),'BinWidth',benchhist.BinWidth);
    
    benchpoints = benchgrid.(fields{i+1});
    %comparepoints = comparegrid.(fields{1+i});
    for j=1:length(benchpoints)
       line([benchpoints(j), benchpoints(j)],get(h,'YLim'),'Color',colorList(1,:)) 
    end
    %for j=1:length(comparepoints)
    %   line([comparepoints(j), comparepoints(j)],get(h,'YLim'),'Color',colorList(2,:)) 
    %end
    %hold off;
    title(captions{i});
end
%legend('Bench',legend_entry);

printfile = 'Results/stateSpaceHistograms';
fig_pos = f.PaperPosition;
f.PaperSize = [fig_pos(3) fig_pos(4)];
print('-dpdf',[printfile,'.pdf']);
print('-depsc',[printfile,'.eps']);

% 0.9315      2.2       1.4
% 0.94        2.45      1.7
% 0.95        2.65      1.8
% 0.96        2.7       1.9