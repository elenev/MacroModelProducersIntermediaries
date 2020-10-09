clear;
respath='./';
outpath='./Results/'; %_NBER
if ~exist('resfile','var')
    resfile='res_20200910_bench_s130';
end
outfile=['GR_',resfile];
%outfile=['GR_unc_',resfile];
grayscale=0;
batchMode=0;
plot_MIT_shocks = false;
%plot_shocks = [1,4,3];
plot_shocks = 1:3;

%load([respath,resfile,'.mat']);
%load([respath,'sim_',resfile,'.mat']);
load([respath,'GR_',resfile,'.mat']);


if plot_MIT_shocks
	mit=load([respath,'MIT_',resfile,'.mat']);
	simseries_mean = [simseries_mean; mit.simseries_mean];
	simseries_median = [simseries_median; mit.simseries_median];
	simseries_std = [simseries_std; mit.simseries_std];
	N_shock = size(simseries_mean,1);
end

close all;

%% Cumulative loss
vars = {'Y','C','X'};
brsel = cellfun(@(x)indexmap.get(x),vars);
cumloss = zeros( 2, length(brsel) );
for ii=2:3
	cumloss(ii-1,:) = sum( simseries_mean{ii}(:,brsel) - simseries_mean{1}(:,brsel) );
end
cumloss_pct = cumloss ./ simseries_mean{1}(1,brsel);
fprintf('Cumulative Change As a Percentage Of Steady State\n');
fprintf('%15s %15s %15s\n', 'Variable', 'Non-Fin Rec', 'Fin Rec');
for jj=1:numel(vars)
	fprintf('%15s %15.2f %15.2f\n', vars{jj}, 100*cumloss_pct(:,jj) );
end

%% IRFS

outpath=[outpath,outfile,'_'];
% file names for graphs (set to empty for no printing)
printfiles={[outpath,'IRF1'],[outpath,'IRF2'],[outpath,'IRF3'],[outpath,'IRF4']};        
%printfiles={[outpath,'MIT1'],[outpath,'MIT2'],[outpath,'MIT3'],[outpath,'MIT4']};

brsel1=[indexmap.get('ZA'),indexmap.get('Y'),...
        indexmap.get('C'), indexmap.get('X')];
levelind1=[0,0,0,0];    
brsel2=[indexmap.get('Lrate'),indexmap.get('brupt'),indexmap.get('Bbkasset'),indexmap.get('Bmlbt'),indexmap.get('Lspr'),... 
        indexmap.get('Imasset'),indexmap.get('Ibklbt'),indexmap.get('fracS'),indexmap.get('BG_g'),indexmap.get('rD')]; 
levelind2=ones(1,10);
brsel3=[indexmap.get('rD'),indexmap.get('rB'),indexmap.get('cS'),...
        indexmap.get('dI_rate'),indexmap.get('BG_g'),indexmap.get('DWL_byY')];
levelind3=ones(1,6);
brsel4=[indexmap.get('dI_eff'),indexmap.get('eI_total')];
levelind4=ones(1,2);
brsel5=[indexmap.get('wB'),indexmap.get('wS')];
levelind5=zeros(1,2);
brsel6=[indexmap.get('ZA'),indexmap.get('Y'),indexmap.get('C'),indexmap.get('X'),...
        indexmap.get('Om'),indexmap.get('Drate'),indexmap.get('brupt'),indexmap.get('Lspr')];
levelind6=[zeros(1,4),ones(1,4)];
    
titles1={'TFP ','Output ','Consumption ','Investment '}; 
titles2={'Loan loss rate','Bank failure rate','Corporate assets','Corporate liabilities','Loan spread',...
         'Financial assets','Financial liabilities','Frac Bonds Held by Savers','Gov debt','Safe rate'};
titles3={'Riskfree rate','Loan rate','WI', ...
         'Net payout rate','Gov debt','DWL'};
titles4={'Net Dividends (Div + Eq Repurch. - Eq Iss.)','Equity Issuance'};
titles5={'Wage B','Wage S'};
titles6={'TFP','Output','Consumption','Investment','Risk Shock','Corp. Def. Rate','Bank Def. Rate','Credit Spread'};

    
tvec=0:NT_sim-1;    
    
brsel_all=[brsel1,brsel2,brsel3,brsel4,brsel5,brsel6];
levelind_all=[levelind1, levelind2, levelind3, levelind4,levelind5, levelind6];
nvar=length(brsel_all);   
brseries_gr=zeros(length(plot_shocks), NT_sim, nvar);
for s=1:numel(simseries_mean)
    for v=1:nvar
        thisv=brsel_all(v);
        if levelind_all(v) 
            brseries_gr(s,:,v) = simseries_mean{s}(:,thisv);
        else
            brseries_gr(s,:,v) = 100*(simseries_mean{s}(:,thisv) ./ repmat(simseries_mean{s}(1,thisv),NT_sim,1)-1);
        end
    end
end

colors={'k-o','b-o','r-o','b-o','c-o'};
     
if usejava('desktop')
    makeImpRes(brseries_gr(plot_shocks,:,1:4),tvec,titles1,colors(plot_shocks),[2,2],[],printfiles{1});   
    makeImpRes(brseries_gr(plot_shocks,:,5:14),tvec,titles2,colors(plot_shocks),[2,5],[],printfiles{2});    
    makeImpRes(brseries_gr(plot_shocks,:,15:20),tvec,titles3,colors(plot_shocks),[2,3],[],printfiles{3});  
	makeImpRes(brseries_gr(plot_shocks,:,21:22),tvec,titles4,colors(plot_shocks),[1,2],[],printfiles{4});
	%makeImpRes(brseries_gr(plot_shocks,:,23:24),tvec,titles4,colors(plot_shocks),[1,2],[],printfiles{4});
	%makeImpRes(brseries_gr(plot_shocks,:,25:end),tvec,titles6,colors(plot_shocks),[2,4],[],printfiles{4});
end


%% Make Welfare Graph
% 
% %N_shock = N_shock-1;
% 
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
% linespec = {'-','-.','--',':'};
% linewidth = 2;
% %YLim=mat2cell(repmat([-2.5,2],1,N_shock),1,2*ones(1,N_shock));
% % YLim=mat2cell(repmat([-10,2],1,N_shock),1,2*ones(1,N_shock));
% YLim=mat2cell(repmat([-1.5,0.5],1,N_shock),1,2*ones(1,N_shock));
% 
% tvec=0:NT_sim-1;
% 
% title_str={'Unconditional','Non-financial recession','Financial recession'};
% if N_shock==4
%    title_str = [title_str,'Unexp. WI shock']; 
% end
% 
% f=figure;
% for s=2:N_shock
%     brsel=[indexmap.get('welfare'),indexmap.get('VB'),indexmap.get('VS')];
% %     if reportLevels==1
% %         simseries_normalized = simseries_mean{s}(:,brsel);
% %     else
%         simseries_normalized = 100*(simseries_mean{s}(:,brsel) ./ repmat(simseries_mean{s}(1,brsel),NT_sim,1)-1);
% %     end
%     
%     subplot(1,N_shock-1,s-1);
%     hold on;
%     h = cell(4,1);
%     for i=1:size(simseries_normalized,2)
%         h{i}=plot(tvec,simseries_normalized(:,i),linespec{i},'Color',colors{i},'LineWidth',linewidth);
%     end
%     set(gca,'XLim',[min(tvec),max(tvec)],'YLim',YLim{s},'XTick',min(tvec):round((tvec(end)-tvec(1))/5):max(tvec),'XGrid','on','GridLineStyle','-');
%     if s==2
%         legend('Welfare','V^B','V^S','Location','southeast');
%         ylabel('Percent Change','Margin',2,'FontSize',14);
%     end
%     title(title_str{s});
%     
%     % text box for right panel
%     %VI_max = max(simseries_normalized(:,3));
%     %VI_min = min(simseries_normalized(:,3));
%     %if VI_min < YLim{s}(1)
%     %    t_bottom=text('Position',[0.5473 -2.3 0], 'String',sprintf('V^I_{min}=%.2f',VI_min),'FontSize',10);
%     %end
%     %if VI_max > YLim{s}(2)
%     %    t_top=text('Position',[3.4467 1.8280 0], 'String',sprintf('V^I_{max}=%.2f',VI_max),'FontSize',10);
%     %end
% end
% 
% hold off;
% printfile=[outpath,outfile,'_print_welfare'];
% %title('Welfare along Median Transition Path');
% D = f.PaperPosition; % Returns 1x4 vector [left right width height]
% f.PaperSize = [D(3) D(4)];
% f.PaperPositionMode = 'auto';
% print('-dpdf',[printfile,'.pdf']);
% print('-depsc',[printfile,'.eps']);
% savefig([printfile,'.fig']);
% 
% %N_shock = N_shock+1;
% 
% if batchMode
%     close all;
%     clear;
% end
% 
% 
