clear;

outpath='./Results/';
N_economy=4; % numbers economies to plot on same graph
respath(1:N_economy)={'./'};
%respath(1:2) = {'../version20170706/'};
bigtitle = 'Bench vs Alternative Pop/Inc Shares: Fin and Non-Fin Recessions';
econ1='bench';
econ2='altshares';
econ3='bench';
econ4='altshares';

%label = {'bench','\xi=\{91%, 93%\}'};
%label = {'no brupt','\xi=93%', '\xi=95%'};
%label = {'bench','no brupt'};
%label = {'bench','\sigma^I=0'};
%label = {'Baseline','\sigma^I=5,\phi_1=0','\sigma^I=5,\phi_1=0,\sigma_\epsilon=0.017,\beta_B=0.9375'};
%label = {'bench','\sigma^I=0','\sigma^I=0, no bank default'};
label = {'bench,fin','altshares,fin','bench,nonfin','altshares,nonfin'};
%label = {'bench','xi=94%','Internal Eq Iss Cost'};
%label = {'bench','no fracS','\phi_1=0'};
%label = {'bench','\phi_1^I=0','\phi_1=0, \sigma^I=0, no bank default'};
%label = {'old','Aug code, July params','Deviation: Macro'};
%label = {'July code, August params','Aug code (94%)'};


if ~exist('resfile','var')
    resfile=['res_20200312_finrecparamcomapre_',econ1,econ2,'_s130'];
    resfile1=['res_20200312_',econ1,'_s130'];
    resfile2=['res_20200312_',econ2,'_s130'];
	resfile3=['res_20200312_',econ3,'_s130'];
	resfile4=['res_20200312_',econ4,'_s130'];
end
outfile=['GR_',resfile];
grayscale=0;
reportLevels=0;
batchMode=0;


load([respath{1},'GR_',resfile1,'.mat']);
simseries_mean_econ{1} = simseries_mean{3}; % financial recession 
indexmap_econ{1} = indexmap;
clear simseries_mean indexmap;

load([respath{2},'GR_',resfile2,'.mat']);
simseries_mean_econ{2} = simseries_mean{3};
indexmap_econ{2} = indexmap;
clear simseries_mean indexmap;

load([respath{3},'GR_',resfile3,'.mat']);
simseries_mean_econ{3} = simseries_mean{2};
indexmap_econ{3} = indexmap;
clear simseries_mean indexmap;

load([respath{4},'GR_',resfile4,'.mat']);
simseries_mean_econ{4} = simseries_mean{2};
indexmap_econ{4} = indexmap;
clear simseries_mean indexmap;

close all;

%% IRFS

outpath2=[outpath,outfile,'_'];
% file names for graphs (set to empty for no printing)
printfiles={[outpath2,'IRF1'],[outpath2,'IRF2'],[outpath2,'IRF3'],[outpath2,'IRFcombined']};        

varsel1 = {'C','X','DWL_byY','rD'};
varsel2 = {'Lrate','brupt','Bbkasset','Bbklbt','Ibklbt','Lspr','BG_g','fracS'};
varsel3 = {'C','X','DWL_byY','Y','gammaG','rD'};

levelind1=[0,0,1,1]; 
levelind2=ones(1,8);
levelind3=ones(1,6);

% brsel1=[indexmap.get('C'),indexmap.get('X'),...
%         indexmap.get('DWL_byY'), indexmap.get('rD')];
%    
% brsel2=[indexmap.get('Lrate'),indexmap.get('brupt'),indexmap.get('Bbkasset'),... 
%         indexmap.get('Bbklbt'),indexmap.get('Ibklbt'),indexmap.get('Lspr')]; 
% 
% brsel3=[indexmap.get('C'),indexmap.get('X'),indexmap.get('DWL_byY'),... 
%         indexmap.get('Bmlbt'),indexmap.get('Ibklbt'),indexmap.get('Lspr')]; 

titles1={'Consumption ','Investment','DWL/GDP','Riskfree rate'}; 
titles2={'Loan loss rate','Bank failure rate','Corporate assets',...
         'Corporate liabilities','Financial liabilities','Loan spread', ...
		 'Government Debt','Saver Share'};
titles3={'Consumption ','Investment','DWL/GDP',...
         'Output','Govt Consumption','Risk-Free Rate'};
     
    
tvec=0:NT_sim-1;    
    
%brsel_all=[brsel1,brsel2,brsel3];
varsel_all=[varsel1,varsel2,varsel3];
levelind_all=[levelind1, levelind2,levelind3];
nvar=length(varsel_all);   
brseries_gr=zeros(N_economy, NT_sim, nvar);
for s=1:N_economy
	indexmap = indexmap_econ{s};
    for v=1:nvar
        thisv=indexmap.get(varsel_all{v});
        if levelind_all(v) 
            brseries_gr(s,:,v) = simseries_mean_econ{s}(:,thisv);
        else
            brseries_gr(s,:,v) = 100*(simseries_mean_econ{s}(:,thisv) ./ repmat(simseries_mean_econ{s}(1,thisv),NT_sim,1)-1);
        end
    end
end

colors={'k-','b-','k:','b:'};

if usejava('desktop')
    makeImpRes(brseries_gr(:,:,1:4),tvec,titles1,colors,[2,2],label,printfiles{1},bigtitle);   
    makeImpRes(brseries_gr(:,:,5:12),tvec,titles2,colors,[2,4],label,printfiles{2},bigtitle);    
    %makeImpRes(brseries_gr(:,:,11:16),tvec,titles3,colors,[2,3],label,printfiles{3},bigtitle);    
	%makeImpRes(brseries_gr(:,:,1:16),tvec,[titles1,titles2,titles3],colors,[4,4],label,printfiles{4},bigtitle);    
end