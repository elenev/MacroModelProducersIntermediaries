clear;
close all;
% file with model
respath='./';
outpath='./Results/';
resfile='res_20200910_bench_s130';

% Initial Economy Config
varlist={'simseries','varnames','params'};
load(['sim_',resfile],varlist{:});

% make HashMap with mapping of names to indices
indexmap=java.util.HashMap;
for i=1:length(varnames)
    indexmap.put(varnames{i},i);
end

% variables for graphs
WI = simseries(:,indexmap.get('WI'));
VI = simseries(:,indexmap.get('VI'));
WS=simseries(:,indexmap.get('WSm'));
WB=simseries(:,indexmap.get('WBm'));
WT=WI+ WB + WS;
WIsh=simseries(:,indexmap.get('WI'))./WT;
WSsh=simseries(:,indexmap.get('WSm'))./WT;
WBsh=simseries(:,indexmap.get('WB'))./WT;
expERP=simseries(:,indexmap.get('expERP'));
expRP=simseries(:,indexmap.get('expRP'))-1;
expRP_check=simseries(:,indexmap.get('expRP_check'));
effect_CovP=simseries(:,indexmap.get('effect_CovP'));
effect_CovP_noBR=simseries(:,indexmap.get('effect_CovP_noBR'));
effect_collatP=simseries(:,indexmap.get('effect_collatP'));
effect_Rf=simseries(:,indexmap.get('effect_Rf'));
Lspr=simseries(:,indexmap.get('Lspr'));
q = simseries(:,indexmap.get('q'));
lamI = simseries(:,indexmap.get('lamR'));
Rf=1./q;
rD=Rf-1;
rB = simseries(:,indexmap.get('rB'));

%expERP_noRf=(expERP+rD)./effect_Rf - Rf.*effect_collatP - rD;
expERP_corr=(expERP+rD).*(1-params.xivec(1).*lamI.^(1/3)) - rD;

%% intermediary
% histogram counts
edgesWI = [min(WIsh),0:0.005:round(prctile(WIsh,99.0),2) ,prctile(WIsh,[99.5,100])];
plotseriesWI = getPlotSeries( edgesWI, WIsh, {Lspr,expERP,expERP_corr}, [1,0] );

printfile=[outpath,'WI_EERhist',resfile];
makePlot( plotseriesWI(1,:), plotseriesWI(end,:), ...
	mat2cell(plotseriesWI(2:3,:),ones(2,1),size(plotseriesWI,2)), ...
	'Intermediary Wealth Share', ...
	{'Credit Spread and Expected Excess Return','Frequency'}, ...
	{'r','b'},...
	{'-','--'},...
	{'Credit Spread','EER'},...
	printfile );

%% entrepreneur
% histogram counts
edgesWB = [min(WBsh), ...
	round(prctile(WBsh,1),2):0.01:round(prctile(WBsh,99),2) max(WBsh)];
plotseriesWB = getPlotSeries( edgesWB, WBsh, {Lspr,expERP,expERP_corr}, [0,0] );

printfile=[outpath,'WB_EERhist',resfile];
makePlot( plotseriesWB(1,:), plotseriesWB(end,:), ...
	mat2cell(plotseriesWB(2,:),ones(1,1),size(plotseriesWB,2)), ...
	'Borrower-Entrepreneur Wealth Share', ...
	{'Expected Excess Return','Frequency'}, ...
	{'r','b'},...
	{'--'},...
	{'EER'},...
	printfile );

%% joint histogram

printfile=[outpath,'WI_WB_hist'];
edges={edgesWI, edgesWB};
data={WIsh, WBsh};
titles={[],'WI share','WB share'};
ticks={[1:length(edges{1})], [1:length(edges{2})] };
makeHist3D(data,ones(size(WIsh)),edges,titles,ticks);
f = gcf;
D = f.PaperPosition; % Returns 1x4 vector [left right width height]
f.PaperSize = [D(3) D(4)];
f.PaperPositionMode = 'auto';
print('-dpdf',[printfile,'.pdf']);
print('-depsc',[printfile,'.eps']);

function plotseries = getPlotSeries( edges, xvar, yvars, chop )
	nbins = length(edges)-1;
	[cnt,~,bins]=histcounts(xvar,edges);

	plotseries=zeros(numel(yvars)+2,nbins);
	plotseries(1,:) = splitapply(@mean,xvar,bins);
	plotseries(2:end-1,:) = splitapply(@mean,[yvars{:}],bins)';
	plotseries(5,:) = cnt/sum(cnt);

	pl0bin1=1+chop(1);
	plmaxbin=nbins-chop(2);
	
	plotseries = plotseries(:,pl0bin1:plmaxbin);
end

function makePlot( xvar, xhist, yvars, xlbl, ylbl, colors, styles, leg, printfile )
	f=figure;
	hold on;
	
	yyaxis right;
	bar(xvar,xhist,0.6,colors{2});
	ylabel(ylbl{2});
	xlabel(xlbl);
	ax=gca;
	ax.YColor=colors{2};

	yyaxis left;
	for ii=1:numel(yvars)
		plot(xvar,yvars{ii},[colors{1},styles{ii}],'LineWidth',2.5);
	end
	ylabel(ylbl{1},'Color',colors{1});
	ax=gca;
	ax.YColor=colors{1};
	ylim( [min(cellfun(@(x)min(x),yvars)), max(cellfun(@(x)max(x),yvars)) ] );

	legend(leg{:});

	D = f.PaperPosition;
	f.PaperSize = [D(3) D(4)];
	f.PaperPositionMode = 'auto';
	print('-dpdf',[printfile,'.pdf']);
	print('-depsc',[printfile,'.eps']);
end