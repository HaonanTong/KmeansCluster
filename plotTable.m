function [ fig ] = plotTable( csv, savepath )
% Nan
%% Read File
T = readtable(csv,...
     'ReadVariableNames',true);
% summary(T);
Data = table2array(T(:,2:end));
% name = table2array(T(:,1));

%% Plot Data Generation
newData = [];
for i = 1:3:21%7 time points; 3 replicates;
   newData = [newData sum(Data(:,i:i+2),2)];
end
newData = 1/3*newData;
tmp = [];
for i = 2:7
   tmp = [tmp log2( newData(:,i)./(newData(:,1)+.01) )];
end
plotData = tmp;
plotData = [ zeros(size(plotData,1),1) plotData ];

%% Plot
fig = figure;
hold on;axis([0 6 -2 5])
plot(x, plotData,'Color','[.4,.4,.4]');
plot(x,mean( plotData),'Color','r','LineWidth',2);
print(fig,savepath,'-dpng'); 


end

