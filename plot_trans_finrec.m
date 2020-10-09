clear;
outpath='./Results/';

%econ{1}='FLbench';
%econ{2}='FLwithdef';
econ{1}='bench';
%econ{2}='benchriskierrisky';
%econ{3}='bench';
%econ{4}='GMSbankDWLeta';
econ{2}='xi9195';

bigtitle= '';

N_economy=numel(econ); % numbers economies to plot on same graph

respath(1:N_economy)={'./'};
%respath(1:2)={'../version20200310/'};
%respath(2)={'../version20200122/'};
%respath(2)={'../version20191008/'};

experdef(1:N_economy)={'20200910'};
%experdef{1}='20170615';
%experdef{2}='20200312';
%experdef([2,4])={'20200212'};

suffix(1:N_economy)={'s130'};
%suffix{1} = 's140';

label = {'bench','\xi=\{91%, 95%\}'};
%label = {'no brupt','\xi=93%', '\xi=95%'};
%label = {'bench','no brupt'};
%label = {'bench','\sigma^I=0'};
%label = {'Baseline','\sigma^I=5,\phi_1=0','\sigma^I=5,\phi_1=0,\sigma_\epsilon=0.017,\beta_B=0.9375'};
%label = {'bench','\sigma^I=0','\sigma^I=0, no bank default'};
%label = {'\sigma^I=8','\sigma^I=0','\sigma^I=0, no bank default'};
%label = {'bench','xi=94%','Internal Eq Iss Cost'};
%label = {'bench','no fracS','\phi_1=0'};
%label = {'bench','\phi_1^I=0','\phi_1=0, \sigma^I=0, no bank default'};
%label = {'Paper','old constraint','new constraint, Phi=0.85', 'new constraint, Phi=0.65, pibar=0'};
%label = {'Paper','new constraint, Phi=0.85', 'new constraint, Phi=0.65, pibar=0'};
%label = {'old','new, Phi=0.85, pibar=0.007','new, Phi=0.85, pibar=0','new, Phi=0.65, pibar=0.007'};
%label = econ(1:N_economy);
%label = {'Model (1)','Model (2)','Benchmark'};
%label = {'Bench','E[\zeta^I]=0.32','E[\zeta^I]=0.32, \eta^I=0.3','E[\zeta^I]=0.32, \eta^I=0.42'};
%label = {'July code, August params','Aug code (94%)'};
%label = {'Old','Flex div, Ext cost','Flex div','Flex Div, Ext Cost, xi=92'};
%label = {'Bench','\sigma=10'};

resfile=['res_',experdef{1},'_finrec_',strcat(econ{:}),'_',suffix{1}];

simseries_mean_econ = cell(1,N_economy);
indexmap_econ = cell(1,N_economy);
for ii=1:N_economy
    tmp_resfile = sprintf('res_%s_%s_%s',experdef{ii},econ{ii},suffix{ii});
	load([respath{ii},'GR_',tmp_resfile,'.mat']);
	simseries_mean_econ{ii} = simseries_mean{3}; % financial recession 
	indexmap_econ{ii} = indexmap;
	clear simseries_mean indexmap;
end


outfile=['GR_',resfile];
grayscale=0;
reportLevels=0;
batchMode=0;

close all;

%% IRFS

outpath2=[outpath,outfile,'_'];
% file names for graphs (set to empty for no printing)
printfiles={[outpath2,'IRF1'],[outpath2,'IRF2'],[outpath2,'IRF3'],[outpath2,'IRFcombined']};        

%varsel1 = {'C','X','Y','DWL_byY'};
%varsel1 = {'Drate','X','DWL_byY','C'};
varsel1 = {'Drate','X','brupt','DWL_byY'};
% varsel2 = {'Lrate','brupt','Bbkasset','cB','Ibklbt','Lspr'};
%varsel1 = {};
varsel2 = {'BG_g','bailout','eI','bind_lamI','Imasset','Imlvg'};
varsel3 = {'rD','fracS','rB','Lspr','equ_total','DWL_byY'};

levelind1=[1,0,1,1]; 
levelind2=ones(1,6);
levelind3=ones(1,6);

% brsel1=[indexmap.get('C'),indexmap.get('X'),...
%         indexmap.get('DWL_byY'), indexmap.get('rD')];
%    
% brsel2=[indexmap.get('Lrate'),indexmap.get('brupt'),indexmap.get('Bbkasset'),... 
%         indexmap.get('Bbklbt'),indexmap.get('Ibklbt'),indexmap.get('Lspr')]; 
% 
% brsel3=[indexmap.get('C'),indexmap.get('X'),indexmap.get('DWL_byY'),... 
%         indexmap.get('Bmlbt'),indexmap.get('Ibklbt'),indexmap.get('Lspr')]; 

% titles1={'Consumption ','Investment','DWL/GDP','Riskfree rate'}; 
% titles2={'Loan loss rate','Bank failure rate','Corporate assets',...
%          'Borrower Consumption','Financial liabilities','Loan spread'};
titles1={'Corp. default rate','Investment','Bank Bankruptcies','DWL/GDP'}; 
%titles1={'Corp. default rate','Investment','DWL/GDP','Consumption '}; 
titles2={'Govt Debt (End)','Bailouts','eI','bind_lamI','Imasset','Imlvg'};     
% titles3={'Consumption ','Investment','Output',...
%          'Corporate liabilities (=fin. assets)','Frac. Held by Savers','qB'};
titles3={'Risk-Free Rate',...
         'Saver frac','Loan rate','Loan spread','Firm Issuance','DWL'};
     
    
tvec=0:NT_sim-1;    
    
varsel_all=[varsel1,varsel2,varsel3];
%varsel_all=[varsel1];
levelind_all=[levelind1,levelind2,levelind3];
%levelind_all=[levelind2,0,0,0,1,1,1,1,1];
nvar=length(varsel_all);   
brseries_gr=zeros(N_economy, NT_sim, nvar);
for s=1:N_economy
	indexmap = indexmap_econ{s};
    for v=1:nvar
        thisv=indexmap.get(varsel_all{v});
        if ~isempty(thisv)
            if levelind_all(v)
                brseries_gr(s,:,v) = simseries_mean_econ{s}(:,thisv);
            else
                brseries_gr(s,:,v) = 100*(simseries_mean_econ{s}(:,thisv) ./ repmat(simseries_mean_econ{s}(1,thisv),NT_sim,1)-1);
            end
        else
            disp(['Variable "',varsel_all{v},'" not defined in economy ',num2str(s),'! Setting to zero.']); 
        end
    end
end

%colors={'k--','b-','r:','c.-'};
colors={'k-o','b-o','r-o','b-o','c-o'};


if usejava('desktop')
    makeImpRes(brseries_gr(:,:,1:4),tvec,titles1,colors,[2,2],label,printfiles{1},bigtitle);   
    makeImpRes(brseries_gr(:,:,5:10),tvec,titles2,colors,[2,3],label,printfiles{2},bigtitle);    
    makeImpRes(brseries_gr(:,:,11:16),tvec,titles3,colors,[2,3],label,printfiles{3},bigtitle);    
	%makeImpRes(brseries_gr(:,:,1:16),tvec,[titles1,titles2,titles3],colors,[4,4],label,printfiles{4},bigtitle); 
	%makeImpRes(brseries_gr(:,:,7:14),tvec,[titles1,titles2,titles3],colors,[2,4],label,printfiles{4},bigtitle); 
end

if batchMode
    close all;
    clear;
end

return;

%% Make Welfare Graph

N_vars=3; % number welfare variables plotted (subplots)

% if ~grayscale
%     colors = {[0, 0, 1], ...
%             [0, 0.5, 0], ...
%             [1, 0, 0], ...
%             [.7,.7,0]};
% else
%     colors = {[0, 0, 0], ...
%             [0.5, 0.5, 0.5], ...
%             [.3, .3, .3], ...
%             [.2,.2,.2]};    
%     printfile = [printfile,'_bw'];
% end

colors={'k','b','r','c'};
linespec = {'-','-.','--',':'};
linewidth = 2;
%YLim=mat2cell(repmat([-2.5,2],1,N_shock),1,2*ones(1,N_shock));
% YLim=mat2cell(repmat([-10,2],1,N_shock),1,2*ones(1,N_shock));
YLim=mat2cell(repmat([-1.75,0.5],1,N_vars),1,2*ones(1,N_vars));

tvec=0:NT_sim-1;

title_str={'welfare','VB','VS'};

f=figure;
brsel=[indexmap.get('welfare'),indexmap.get('VB'),indexmap.get('VS')];
for s=1:N_vars
    
    subplot(1,N_vars,s);
    hold on;
    h = cell(N_economy,1);
    for i=1:N_economy
%     if reportLevels==1
%         simseries_normalized = simseries_mean{s}(:,brsel);
%     else
        simseries_normalized = 100*(simseries_mean_econ{i}(:,brsel(s)) ./ repmat(simseries_mean_econ{i}(1,brsel(s)),NT_sim,1)-1);
%     end        
        h{i}=plot(tvec,simseries_normalized,linespec{i},'Color',colors{i},'LineWidth',linewidth);
    end
    set(gca,'XLim',[min(tvec),max(tvec)],'YLim',YLim{s},'XTick',min(tvec):round((tvec(end)-tvec(1))/5):max(tvec),'XGrid','on','GridLineStyle','-');
    if s==1
        legend(label{:},'Location','southeast');
        ylabel('Percent Change','Margin',2,'FontSize',14);
    end
    title(title_str{s});
    
    % text box for right panel
    V_max = max(simseries_normalized);
    V_min = min(simseries_normalized);
    if V_min < YLim{s}(1)
       t_bottom=text('Position',[10 -1.25 0], 'String',sprintf('min=%.2f',V_min),'FontSize',10);
    end
    if V_max > YLim{s}(2)
       t_top=text('Position',[10 -1.25 0], 'String',sprintf('max=%.2f',V_max),'FontSize',10);
    end    
end

hold off;
printfile=[outpath,outfile,'_print_welfare'];
%title('Welfare along Median Transition Path');
D = f.PaperPosition; % Returns 1x4 vector [left right width height]
f.PaperSize = [D(3) D(4)];
f.PaperPositionMode = 'auto';
print('-dpdf',[printfile,'.pdf']);
print('-depsc',[printfile,'.eps']);
savefig([printfile,'.fig']);


