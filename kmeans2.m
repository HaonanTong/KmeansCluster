%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@Haonan Tong
%PGRP Cluster in kmeans.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Have Differential Expressed genes 1464;
% Data Initalization
T = readtable('Arbdpss_Dffrntl_xprss_gns.csv',...
    'ReadVariableNames',true);
summary(T);
Data = table2array(T(:,2:end));
name = table2array(T(:,1));
% Clustering Data
logData = log(Data);
%normalizedData = log(Data./repmat(Data(:,1),1,size(Data,2)));

% Plot Data Generation
newData = [];

for i = 1:3:21%7 time points; 3 replicates;
   newData = [newData sum(Data(:,i:i+2),2)];
end
newData = 1/3*newData;

tmp = [];
for i = 2:7
   tmp = [tmp log( newData(:,i)./(newData(:,1)+.01) )];
end
plotData = tmp;
plotData = [ zeros(size(plotData,1),1) plotData ];
%% Kmeans
% Kmeans initialization
nclusters = 4; ngenes = zeros(nclusters,1);
str = sprintf('ncluster%d',nclusters);
mkdir('./Result',str);

% Kmeans operation
idx = kmeans(logData,nclusters);


%% Result Analysis
% Print out # of genes for each cluster
for i = 1 : nclusters
    ngenes(i) = sum(idx==i);
end
display(ngenes);

% Print out each clusterto files
table = cell(nclusters,1);
for i = 1 : nclusters
    table{i} = T(idx==i,:);
    writetable(table{i},sprintf('./Result/ncluster%d/cluster%d.csv',nclusters,i))
end


%% Plot
x = 0 : 1 : 6;
for i = 1 : nclusters
    % plot cluster i
    figure;hold on;
    plot(x, plotData(idx==i,:),'Color','r');
end
    



