myDir = '.';
myFiles = dir(fullfile(myDir,'cluster*.csv'));

mtrx_mplotData = [];
for i = 1 : length(myFiles)
    [ fig, ngene, expr,plotData, mplotData, agis, agis_new ] = ...
        f_plotTable2( myFiles(i).name, [], 'Mean Plot' );
    mtrx_mplotData = [ mtrx_mplotData ; mplotData];
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