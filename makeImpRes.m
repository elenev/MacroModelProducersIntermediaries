function dum=makeImpRes(simseries,tvec,titles,colors,format,lbl,printfile,varargin)

    dum=[];
    sims=size(simseries);
    if length(sims)==2
        npp=1;
        npl=sims(2);
        nper=sims(1);
        simseries=reshape(simseries,[1,size(simseries)]);
    else
        npp=sims(1);
        npl=sims(3);
        nper=sims(2);
    end
    if npl~=prod(format)
        disp('format does not match number of plots');
        return;
    end
    figure;
    for p=1:npl
        h=subplot(format(1),format(2),p);
        hold on;
        for pp=1:npp
            plot(tvec,squeeze(simseries(pp,:,p)),colors{pp},'LineWidth',1.5);
        end
        title(h,titles{p});
        set(h,'XLim',[min(tvec),max(tvec)],'XTick',min(tvec):round((tvec(end)-tvec(1))/5):max(tvec),'XGrid','on','GridLineStyle','-');
        set(h,'FontSize',16);
    end
    if ~isempty(lbl)
           l=legend(lbl{:},'Location','best');
           l.FontSize=14;
	end
	if nargin>7
		sgtitle(varargin{1});
	end
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[4*format(2) 4*format(1)]);
    set(gcf,'PaperPosition',[0 0 4*format(2) 4*format(1)]);
    set(gcf,'PaperPositionMode','manual');
    %set(gcf,'PaperSize',[10*format(1), 3*format(2)]);
    if ~isempty(printfile)
        print('-dpdf',[printfile,'.pdf']);
        print('-depsc',[printfile,'.eps']);
    end

end