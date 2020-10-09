% Plot Policy Functions, Appendix Figure B.1

global idx transform variableLabel variableName KB BG state;

resfile = 'res_20200910_bench_s130.mat';
printfile = 'Results/polPlot';
state = 5;
KB = 2.13;
BG = 0.71;

idx = cell(2,1);
transform=idx;
variableLabel=idx;
variableName=idx;
idx{1} = 2; % Investment
transform{1} = 0;
variableName{1} = 'X';
variableLabel{1} = 'Investment';
idx{2} = 13; % Capital Requirement Multiplier
transform{2} = 3; % convert to actual multiplier
variableLabel{2} = 'Cap Req Multiplier';
variableName{2} = '$\lambda^I$';
res=load(resfile);

global mobj;
mobj = res.mobj;

f=figure;
for ii=1:2
	subplot(1,2,ii);
	makeOnePlot( f, ii );
end

printFig(f,printfile);

for ii=1:2
	ftmp = figure;
	makeOnePlot( ftmp, ii );
	printFig( ftmp, sprintf('%s%d', printfile, ii) );
end


function makeOnePlot( fig, counter )
	global idx transform variableName variableLabel state KB BG mobj;
	figure(fig);
	mobj.Pfct.plot3D(idx{counter},[3,4],[state,KB,BG],transform{counter});
	xlabel('$L^B$','interpreter','latex');
	ylabel('$W^I$','interpreter','latex');
	zlabel(variableName{counter},'interpreter','latex');
	title(variableLabel{counter});
end

function printFig( fig, printfile )
	figure(fig);
	set(gcf,'PaperUnits','inches');
	set(gcf,'PaperSize',[8 4]);
	set(gcf,'PaperPosition',[0 0 8 4]);
	set(gcf,'PaperPositionMode','manual');
	%set(gcf,'PaperSize',[10*format(1), 3*format(2)]);
	if ~isempty(printfile)
		print('-dpdf',[printfile,'.pdf']);
		print('-depsc',[printfile,'.eps']);
	end
end