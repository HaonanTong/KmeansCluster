close all;
myDir = '.';
myFiles = dir(fullfile(myDir,'cluster*.csv'));

mtrx_mplotData = [];
for i = 1 : length(myFiles)
    [ fig, ngene, expr,plotData, mplotData, agis, agis_new,p{i} ] = ...
        f_plotTable2( myFiles(i).name, [], 'Mean Plot' );
    mtrx_mplotData = [ mtrx_mplotData ; mplotData];
    title('');
end
figure;
plot(mtrx_mplotData','LineWidth',3);
grid on;
    xticks(1:7)
    xticklabels({'0','0.25','0.5','1','4','12','24'})
    title('Plot mean of four clusters derive by kmeans','FontSize',14)
    xlabel('Ethylene treatment(hrs)');
    ylabel('Expression-log2ratio(reference at 0 hrs)');
    set(gca,'fontsize',14);
   
    set(p{2},'Color','red')
    set(p{1},'Color','blue')
    set(p{4},'Color','yellow')
    set(p{3},'Color','green')

figure; 
for i = 1:4
    h(i)=subplot(2,2,i);    %ntitle(''); 
        xlabel('Ethylene treatment(hrs)');
    ylabel({'Expression', 'log2ratio(reference at 0 hrs)'});
end

i=0;
for k = [2 4 1 3]
      i=i+1;
      copyobj(allchild(get(figure(k),'CurrentAxes')),h(i));
      axis(h(k),[0 6 -2 5]);
      set(h(k),'xtick',0:6);
      set(h(k),'XTickLabel',{'0','0.25','0.5','1','4','12','24'});          
end
% suptitle('Clustering Result of KMEANS')
