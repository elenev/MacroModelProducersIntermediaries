function dum=makeHist3D(expvec,freq,edges,titles,tick)

dum=[];
expbin=cell(2,1);
for d=1:2
    thisexp=expvec{d};
    thisedges=edges{d};
    thisedges(1)=thisedges(1)-1e-4;
    thisedges(end)=thisedges(end)+1e-4;
    [~,thisbin]=histc(thisexp,thisedges);
    expbin{d}=thisbin;    
end
nb=[length(edges{1}),length(edges{2})];
jointfreq=accumarray([expbin{1},expbin{2}],freq,nb)'/sum(freq);
figure;
contour(jointfreq,20,'LineWidth',2);
hold on;
set(gca,'XTick',tick{1},'XTickLabel',cellstr(num2str(edges{1}(tick{1})')),...
    'YTick',tick{2},'YTickLabel',cellstr(num2str(edges{2}(tick{2})')),'FontSize',12);
if ~isempty(titles{1})
    title(titles{1});
end
set(get(gca,'XLabel'),'String',titles{2},'FontSize',12);
set(get(gca,'YLabel'),'String',titles{3},'FontSize',12);


end