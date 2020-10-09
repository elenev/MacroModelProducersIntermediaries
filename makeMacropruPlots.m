% Make macropru plots
clear; clc; close all;

% Economies
experdef = '20200910';
suffix = 's130';
bench_econ = 'bench';
other_econs = {'xi75','xi80','xi81midxigrid','xi82midxigrid','xi83midxigrid',...
	'xi84midxigrid','xi85','xi86midxigrid','xi87midxigrid','xi88midxigrid',...
	'xi89midxigrid','xi90','xi91benchgrid','xi92','xi93benchgrid',...
	'xi94','xi95benchgrid','xi96','xi97benchgrid'};
econs = [{bench_econ},other_econs];

% Define series
% Each row represents an economys-specific statistic to plot and contains 4 elements
% (1) unique identifier, (2) variable, (3) caption, (4) statistic (as
% defined in statsexog_XXXXXX or cvwelfare for CV welfare)
% (1) and (2) can be identical, but (1) must be unique even when (2) is not
% e.g. when plotting mean AND std of the same variable
sercell = { 'Lrate', 'Lrate','Loan Loss Rate','mean';
			'brupt', 'brupt','Bank Brupt','mean';
			'Y', 'Y','Y','mean';
			'bI_byY', 'bI_byY','Deposits / Y','mean';
			'VB', 'VB','Borrower','mean';
			'VS', 'VS','Saver','mean';
			'cvwelfare', 'cvwelfare','Ex-Ante (CV)','mean';
			'C_gr', 'C_gr','Cons Growth','std';
			'Y_gr', 'Y_gr','GDP Growth','std';
			'X_gr', 'X_gr','Inv Growth','std';
			'DivP', 'DivP','Producers','mean';
			'dI_eff', 'dI_eff','Intermediaries','mean';};
			%'welfare', 'welfare','Pop-Weighted','mean';
			
% Which series IDs are changes relative to benchmark
chvar = {'VB';'VS'};
% Which series IDs are to be expressed in percent
pctvar = setdiff( sercell(:,1), {'Y'; 'bI_byY'; 'DivP'; 'dI_eff'; 'cvwelfare'} );

% Graph layout
% This is a cell of cells.
% Outer cell is N x 2, where N is the number of plots to create
% First column will be plotted on the left axis, second column (if not
% empty) will be plotted on the right axis
% Each element is a cell of series names to plot on that plot and axis
graph_sets = {  {'Lrate'}, {'brupt'};
				{'Y'}, {'bI_byY'};
				{'VB', 'VS'}, {};
				{'cvwelfare'}, {};
				{'C_gr'}, {'X_gr'};
				{'DivP'}, {'dI_eff'} };

titles = {'Financial Fragility'; 'Size of the Economy'; ...
	'Change Relative to Benchmark'; 'Aggregate Welfare (CEV)'; ...
	'Volatility'; 'Net Payouts'};

align_axes_at = {nan; nan; nan; 0; nan; nan};

xirange = [0.73,0.99];

% Graph parameters
lineWidth = 2;
markerSize = 17;
styles = {'.-','.--','.:'}; % At least as many as would appear on one plot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FULLY DYNAMIC BELOW THIS LINE %%%

sertable = cell2table( sercell(:,2:end) );
sertable.Properties.VariableNames = {'var','cap','stat'};
sertable.Properties.RowNames = sercell(:,1)';


%% Assemble data
graph_pairs_vec = sertable.var;
stats_vec = sertable.stat;
simidx_mean = find( strcmp(stats_vec,'mean') & ~strcmp(graph_pairs_vec,'cvwelfare') );
simidx_std = find( strcmp(stats_vec,'std') & ~strcmp(graph_pairs_vec,'cvwelfare') );
cvwelfare_idx = find( strcmp(graph_pairs_vec,'cvwelfare') );

data = zeros( numel(econs), 1 + numel(graph_pairs_vec) );
use_excel = false;
if use_excel
	for econCounter = 1:numel(econs)
		econ = econs{econCounter};
		file = ['Results/statsexog_res_',experdef,'_',econ,'_',suffix,'.xls'];
		uncond = readtable( file, 'Sheet','Sheet1', 'ReadRowNames', true);
		params = readtable( file, 'Sheet','params', 'ReadRowNames', true);
		data( econCounter, 1 ) = params.values_2('xivec');
		data( econCounter, 1 + simidx_mean ) = uncond.mean( graph_pairs_vec(simidx_mean) );
		data( econCounter, 1 + simidx_std ) = uncond.std( graph_pairs_vec(simidx_std) );
	end
else
	for econCounter = 1:numel(econs)
		econ = econs{econCounter};
		file = ['sim_res_',experdef,'_',econ,'_',suffix,'.mat'];
		warning off
		sim = load(file);
		warning on
		sim.simseries(:,sim.indexmap.get('X_gr')) = winsor( sim.simseries(:,sim.indexmap.get('X_gr')), 0.3 );
		
		data( econCounter, 1 ) = sim.params.xivec(1);
		data( econCounter, 1 + simidx_mean ) = mean( sim.simseries( :, ...
			cellfun(@(x)sim.indexmap.get(x), graph_pairs_vec(simidx_mean) ) ) );
		data( econCounter, 1 + simidx_std ) = std( sim.simseries( :, ...
			cellfun(@(x)sim.indexmap.get(x), graph_pairs_vec(simidx_std) ) ) );
	end
end

%% CV welfare
welfare = readtable( ['Results/welfare_',experdef,'_',bench_econ,'_',suffix,'.xls'], ...
	'Sheet','NumClust 10','ReadRowNames',true);
[~,ia,ib] = intersect( econs, welfare.Properties.RowNames, 'stable' );
data(ia,1+cvwelfare_idx) = welfare.Total(ib);

%% Arrange data
tabdata = array2table(data);
tabdata.Properties.VariableNames = [{'xi'}; sertable.Properties.RowNames]';
tabdata.Properties.RowNames = econs;
tabdata = sortrows(tabdata,'xi');

%% Transform data
% Replace welfare with changes from benchmark
tabdata{:,chvar} = tabdata{:,chvar} ./ tabdata{'bench',chvar} - 1;

% Percent
tabdata{:,pctvar} = 100 * tabdata{:,pctvar};
pctExtra = " (%)";
sertable.pctExtra(:) = "";
sertable.pctExtra( ismember( sertable.Properties.RowNames, pctvar) ) = pctExtra;
sertable.legendEntries = strcat( sertable.cap, sertable.pctExtra );

%% Make plots
close all;
for graphCounter = 1:size(graph_sets,1)
	figure; hold on;
	set(gcf,'DefaultLineLineWidth',lineWidth);
	
	twoAxes = ~isempty( graph_sets{ graphCounter, 2 } );
	if twoAxes
		title( titles{graphCounter} );
		yyaxis left;
	end
	
	series_left = graph_sets{ graphCounter, 1 };
	ax=gca;
	initColorOrder = ax.ColorOrderIndex;
	for seriesCounter = 1:numel( series_left )
		h=plot( tabdata.xi, tabdata.(series_left{seriesCounter}), styles{seriesCounter}, 'MarkerSize', markerSize );
		plot( tabdata{'bench','xi'}, tabdata{'bench',series_left{seriesCounter}}, ...
			'.', 'MarkerSize', 2*markerSize, ...
			'MarkerEdgeColor', h.Color, 'MarkerFaceColor', h.Color, ...
			'HandleVisibility','off');
		ax.ColorOrderIndex = max( initColorOrder, ax.ColorOrderIndex - 1);
	end
	Lax_ylim = ax.YLim;
	
	if ~twoAxes
		ylabel( titles{graphCounter} );
		series_right = {};
	else
		series_right = graph_sets{ graphCounter, 2};
		yyaxis right;
		ax=gca;
		initColorOrder = ax.ColorOrderIndex;
		for seriesCounter = 1:numel( series_right )
			h=plot( tabdata.xi, tabdata.(series_right{seriesCounter}), ...
				styles{numel(series_left)+seriesCounter}, 'MarkerSize', markerSize );
			plot( tabdata{'bench','xi'}, tabdata{'bench',series_right{seriesCounter}}, ...
				'.', 'MarkerSize', 2*markerSize, ...
				'MarkerEdgeColor', h.Color, 'MarkerFaceColor', h.Color, ...
				'HandleVisibility','off');
			ax.ColorOrderIndex = max( initColorOrder, ax.ColorOrderIndex - 1);
		end
		
		if ~isnan(align_axes_at{graphCounter})
			align = align_axes_at{graphCounter};
			
			Lratio = ( Lax_ylim(2) - align ) / ( align - Lax_ylim(1) );
			
			Rax_ylim = ax.YLim;
			Rax_ytick = ax.YTick;
			Rratio = ( Rax_ylim(2) - align ) / ( align - Rax_ylim(1) );
			
			if Rratio < 0 || Lratio < 0
				warning('Axis alignment point is outside YLims');
			else
				if Lratio > Rratio
					% Adjust right top
					ax.YLim(2) = align + Lratio * ( Rax_ylim(1) - align);
					step = ( Rax_ytick(2) - Rax_ytick(1) ) * round( range(ax.YLim) / range( Rax_ylim ) );
					ax.YTick = Rax_ylim(1) : step : ax.YLim(2);
				else
					% Adjust right bottom
					ax.YLim(1) = align - ( Rax_ylim(2) - align ) / Lratio;
					step = ( Rax_ytick(2) - Rax_ytick(1) ) * round( range(ax.YLim) / range( Rax_ylim ) );
					tmp = Rax_ylim(2) : -step : ax.YLim(1);
					ax.YTick = tmp(end:-1:1);
				end
				
				plot( xirange, [align, align], 'k:' );
			end
		end
	end	

	xlim(xirange);
	series = [series_left, series_right];
	legend( sertable{ series, 'legendEntries' }', ...
			'location', 'southoutside', 'Orientation','horizontal');
		
	% For some reason, changing font size messes with the y axes
	% Reset after
	if twoAxes
		yyaxis right;
		y2 = ylim;
		yyaxis left;
	end
	y1 = ylim;
	set(gca,'FontSize',14);	
	ylim(y1);
	if twoAxes
		yyaxis right;
		ylim(y2);
	end
	
	figfile = ['macropru',sprintf('__%s',series{:})];
	set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[4 4]);
    set(gcf,'PaperPosition',[0 0 4 4]);
    set(gcf,'PaperPositionMode','manual');
	print('-dpdf',['Figures/',figfile,'.pdf']);
	print('-depsc',['Figures/',figfile,'.eps']);
end

function out = winsor( in, p )
	cutoffs = prctile( in, [p, 100-p] );
	out = in;
	out( out < cutoffs(1) ) = cutoffs(1);
	out( out > cutoffs(end) ) = cutoffs(end);
end

 