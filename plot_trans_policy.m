respath='./';
outpath='./Results/';


%PT_files = dir('PT_res_20170615_*_s140.mat');
%PT_files = PT_files.name;

resfiles = {'res_20200910_xi85_s130'};
labels={'\xi=0.85'};

%resfiles = {'res_20200910_xi95benchgrid_s130'};
%labels={'\xi=0.95'};
		 
grayscale=0;
reportLevels=1;
batchMode=0;

%load([respath,resfile,'.mat']);
%load([respath,'sim_',resfiles{1},'.mat'],'indexmap','varnames');
%load([respath,'GR_',resfile,'.mat']);

% file names for graphs (set to empty for no printing)
printfiles={[outpath,'PT_',resfiles{1}]}; 

close all;

%% Make Welfare Graph

NT_sim = 0;
N_vars = 0;
N_exper = numel(resfiles);

tmpcl = cell(N_exper,1);

for i=1:N_exper
	tmp=load(['PT_',resfiles{i},'.mat']);
	tmp_simseries = tmp.simseries_mean{1};
	
	if size(tmp_simseries,1) ~= NT_sim && NT_sim ~= 0
		error('Number of periods doesn''t match');
	end
	if size(tmp_simseries,2) ~= N_vars && N_vars ~= 0
		error('Number of variables doesn''t match');
	end
	NT_sim=size(tmp_simseries,1);
	N_vars=size(tmp_simseries,2);
	tmpcl{i}=tmp_simseries;
end
indexmap = tmp.indexmap;
simseries_mean = reshape([tmpcl{:}],size(tmpcl{1},1),size(tmpcl{1},2),numel(resfiles));
simseries_mean = permute(simseries_mean,[3,1,2]);

tvec=0:NT_sim-1;

colors={'b-o','k--'};

brsel1 = [indexmap.get('VB'), indexmap.get('VS')];
titles1 = {'Value Fcn, B', 'Value Fcn, S'};
levelind1 = [0,0];
brsel2 = [indexmap.get('C'), indexmap.get('Y')];
titles2 = {'Consumption', 'Output'};
levelind2 = [0,0];

if strcmp(colors{end},'k--')
	simseries_mean(end+1,:,:) = repmat( simseries_mean(1,1,:), 1, NT_sim, 1 );
end

brsel = [brsel1,brsel2];
titles = [titles1,titles2];
levelind_all = [levelind1,levelind2];

brseries_gr = simseries_mean(:,:,brsel);
for v=1:length(brsel)
	thisv=brsel(v);
	if levelind_all(v)
		brseries_gr(:,:,v) = simseries_mean(:,:,thisv);
	else
		brseries_gr(:,:,v) = 100*(simseries_mean(:,:,thisv) ./ simseries_mean(:,1,thisv)-1);
	end
end

makeImpRes(brseries_gr,tvec,titles,colors,[2,2],labels,printfiles{1});
%makeImpRes(simseries_mean(:,:,brsel2),tvec,titles2,colors,[1,2],labels,printfiles{2});