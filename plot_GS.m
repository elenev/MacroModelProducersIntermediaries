% Make parameter elasticity plots
clear; clc; close all;

% Economies
experdef = '20200910';
suffix = 's130';
bench_econ = 'bench';

read_from_latex = true;
printfile = 'Figures/GS';
codepath = './';

%% Identify parameters
initialVars = who;
run(sprintf([codepath,'experdef_%s.m'],experdef));
clearvars('-except', initialVars{:}, 'expers_gs', 'forwardDifferenceOnly');

%% Load moments
dbfile = [codepath,'Results/simres.mat'];
load(dbfile,'mapDict');


%% Define moments
if read_from_latex
	texfile = 'intermed_production_20200916auto';
	optionalDefaults = {'mean','u',bench_econ,'2'};
	
	tablePattern = '\\caption{Calibrated Parameters}(.*?)\\label{table:calibration}.+?\\end{table}';
	linePattern = '^\s+\$.*?\\\\';
	pattern = ['\$(?<paramsymbol>.*?)\$\s+&\s+(?<paramdesc>.*?)\s+', ...
		'&\s+\\param{(?<paramname>.*?)}(?<paramidx>\[[0-9]+\])?.*?\s+', ...
		'&\s+(?<momentdesc>.*?)\s+&\s+',...
		'\\simres\{(?<var>.*?)\}(?<optional>(\[.*?\]){0,4})'];
	
	paper = fileread([texfile,'.tex']);
	calibTableText = regexp(paper,tablePattern,'tokens','once');
	
	matches = regexp(calibTableText{1},linePattern,'match','lineanchors');
	
	str = cell(numel(matches),1);
	description = str;
	paramCell = cell(numel(matches), 4);
	
	for ii=1:numel(matches)
		names = regexp(matches{ii},pattern,'names','once');
		[~,key] = replaceVal(names.var,names.optional, 'sim', ...
				matches{ii}, mapDict, optionalDefaults);
		str{ii} = extractAfter( key, ['sim - ',bench_econ,' - ']);
		description{ii} = names.momentdesc;
		paramCell(ii,:) = {names.paramsymbol, names.paramdesc, ...
			names.paramname, names.paramidx(2:end-1)};
	end
	
	paramTable = cell2table(paramCell);
	paramTable.Properties.VariableNames = {'symbol','desc','name','idx'};
	
	% Replace AStarget with fracStarget
	paramTable.name{strcmp(paramTable.name,'AStarget')}='fracStarget';
	
	momentTable = table(str,description);
	momentTable.Properties.RowNames = momentTable.str;
else
	moments = { ...
		'logY','AC';
		'logY','std';
		'Drate','mean';
		'prodIQR','std';
		'logX','std';
		'wagebill_byY','mean';
		'X_byY','mean';
		'KB_byY','mean';
		'LGD','mean';
		'Bbklvg','mean';
		'brupt','mean';
		'Lspr','mean';
		'fracS','mean';
		'fracS','std';
		'eI_rate','mean';
		'rD','mean';
		};
	momentStrings = cell(size(moments,1),1);
	for ii=1:size(moments,1)
		momentStrings{ii} = ['u - ',moments{ii,1},' - ',moments{ii,2}];
	end
end

%% Moment table description shortening
pattern = '( [0-9]+\-[0-9]+|\(.*?\)+|holdings )';
for ii = 1:size(momentTable,1)
	momentTable.description{ii} = ...
		regexprep( momentTable.description{ii}, pattern, '' ); 
end

%% Find parameters for which we actually have results
expertable = cell2table(expers_gs);
expertable.Properties.VariableNames = {'expername', 'deviations', ...
	'name', 'idx', 'suffix', 'direction'};
expertable_U = expertable( expertable.direction=="U", : );
expertable_D = expertable( expertable.direction=="D", : );

paramTable.idx( cellfun('isempty',paramTable.idx) ) = {'1'};
paramTable.idx = str2double( paramTable.idx );

mergedParamTable = innerjoin( paramTable, expertable_U, 'Keys', {'name','idx'} );
expertable_D.Properties.VariableNames{'expername'} = 'expername_D';
expertable_D.Properties.VariableNames{'deviations'} = 'deviations_D';
mergedParamTable = outerjoin( mergedParamTable, expertable_D, 'Keys', {'name','idx'}, ...
	'Type', 'left', 'RightVariables', {'expername_D', 'deviations_D'}, 'MergeKeys', true);
mergedParamTable.Properties.RowNames = cellfun( @(x) ...
	extractBefore( x, length(x) ), mergedParamTable.expername, ...
	'UniformOutput', false);
mergedParamTable.value = zeros( size(mergedParamTable,1), 1);
mergedParamTable.value_D = zeros( size(mergedParamTable,1), 1);

%% Create empty tables
targets_blank = array2table( zeros(size(momentTable,1), size(mergedParamTable,1)) );
targets_blank.Properties.VariableNames = mergedParamTable.Properties.RowNames;
targets_blank.Properties.RowNames = momentTable.str;

targets = cell(2,1);
targets{1} = targets_blank;
targets{2} = targets_blank;

%% Fill tables
for jj=1:size(targets_blank,2)
	econ = targets_blank.Properties.VariableNames{jj};
	param = mergedParamTable.name{jj};
	idx = mergedParamTable.idx(jj);
	fwd_only = ismember( expertable.name{jj}, forwardDifferenceOnly );
	for ii=1:size(targets_blank,1)
		moment = targets_blank.Properties.RowNames{ii};
		targets{2}{ii,jj} = mapDict(['sim - ',econ,'U - ',moment]);
		if fwd_only
			targets{1}{ii,jj} = mapDict(['sim - ',bench_econ,' - ',moment]);
		else
			targets{1}{ii,jj} = mapDict(['sim - ',econ,'D - ',moment]);
		end
	end
	
	% Parameters
	is_square = regexp(param,'sig2') > 0;
	is_beta = regexp(param,'beta') > 0;
	mergedParamTable.value(jj) = mergedParamTable.deviations{jj}{1,2}(idx);
	if is_square
		mergedParamTable.value(jj) = sqrt( mergedParamTable.value(jj) );
	end
	if is_beta
		mergedParamTable.value(jj) = -log( mergedParamTable.value(jj) );
		mergedParamTable.symbol{jj} = ['-\log ', mergedParamTable.symbol{jj}];
	end
	if fwd_only
		mergedParamTable.value_D(jj) = mapDict(sprintf('param - %s - %s',param,idx));
	else
		mergedParamTable.value_D(jj) = mergedParamTable.deviations_D{jj}{1,2}(idx);
		if is_square
			mergedParamTable.value_D(jj) = sqrt( mergedParamTable.value_D(jj) );
		end
		if is_beta
			mergedParamTable.value_D(jj) = -log( mergedParamTable.value_D(jj) );
		end
	end
end

%% Compute elasticities
paramIdx=1:size(targets_blank,2);
numerator_levels = (targets{2}{:,paramIdx} - targets{1}{:,paramIdx});
numerator_pct = numerator_levels * 2 ./ (targets{2}{:,paramIdx} + targets{1}{:,paramIdx});
denominator_levels = mergedParamTable.value( targets_blank.Properties.VariableNames )' - ...
	mergedParamTable.value_D( targets_blank.Properties.VariableNames )';
denominator_pct = denominator_levels * 2 ./ ...
	( mergedParamTable.value( targets_blank.Properties.VariableNames )' + ...
	  mergedParamTable.value_D( targets_blank.Properties.VariableNames )' );
targets_GS_levels = targets{1};
targets_GS_levels{:,paramIdx} = numerator_levels ./ denominator_levels;
targets_GS_pct = targets{1};
targets_GS_pct{:,paramIdx} = numerator_pct ./ denominator_pct;

T = targets_GS_pct;
T.Properties.RowNames = momentTable.description( targets_blank.Properties.RowNames );


%% Plot
close all;
rows = 4;
columns = 4;
total = rows * columns;
pages = ceil( size(T,2) / total );

exclude_targets = [17:23];
include_targets = setdiff(1:size(T,1),exclude_targets);

moments = categorical(T.Properties.RowNames(include_targets));
moments = reordercats(moments,T.Properties.RowNames(include_targets(end:-1:1)));

moments_num = categorical(include_targets);
moments_num = reordercats(moments_num,include_targets(end:-1:1));

for pp=1:pages
	f=figure;
	
	f.PaperOrientation = 'portrait';
	f.PaperPositionMode = 'manual';

	offset = 1.8;
	f.PaperPosition(1) = f.PaperPosition(1) + offset;
	f.PaperPosition(2) = 0;
	f.PaperPosition(1) = offset;
	f.PaperPosition(4) = f.PaperPosition(4)*2.3;
	fig_pos = f.PaperPosition;
	f.PaperSize = [fig_pos(3) + offset fig_pos(4)];
	
	for i=1:total
		idx = total*(pp-1) + i;
		if idx>size(T,2)
			break;
		end
		subplot(rows,columns,i);
		if mod(i,columns)==1
			barh( moments, max(0,[T{include_targets,idx}]) );
			hold on;
			barh( moments, min(0,[T{include_targets,idx}]) );
			hold off;
		else
			barh( moments_num, max(0,[T{include_targets,idx}]) );
			hold on;
			barh( moments_num, min(0,[T{include_targets,idx}]) );
			hold off;
			set(gca,'YTickLabel',[]);
		end

		label = ['$',mergedParamTable.symbol{T.Properties.VariableNames{idx}},'$'];
		title(label,'interpreter','latex');
	end

	print('-dpdf',sprintf('%s_%d.pdf',printfile,pp));
	print('-depsc',sprintf('%s_%d.eps',printfile,pp));
end
