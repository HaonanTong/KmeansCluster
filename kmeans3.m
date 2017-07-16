%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@Haonan Tong
%PGRP Cluster in kmeans.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Have Differential Expressed genes 1464;
close all;clc;
%% Load Data
T = readtable('Arbdpss_Dffrntl_xprss_gns.csv',...
    'ReadVariableNames',true);
summary(T);
Data = table2array(T(:,2:end));
name = table2array(T(:,1));
% T.Properties.VariableNames

%% log-transformed RPKM values
% add a small constant to those RPKM values at zeros
eps = .001;
[nlen,ncol] = size(Data);
for i = 1 : nlen
    for j = 1 : ncol
        if( Data(i,j) == 0 )
            Data(i,j) = eps;
        end
    end
end
% log transformation
logData = log(Data);

%% Get feature vector
% Take the average transformed data over replicates for each treatment and
% each gene
FeaMtrx = [];

for i = 1:3:21%7 time points; 3 replicates;
   FeaMtrx = [FeaMtrx sum(logData(:,i:i+2),2)];
end
% sample treatment means for each gene
FeaMtrx = 1/3*FeaMtrx;
FeaMtrx = FeaMtrx - repmat(mean(FeaMtrx,2),1,size(FeaMtrx,2));

%% Kmeans
% Kmeans initialization
nclusters = 20; ngenes = zeros(nclusters,1);
str = sprintf('ncluster%d',nclusters);
mkdir('./Result',str);

% Kmeans operation
idx = kmeans(FeaMtrx,nclusters);


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
% Plot Data Generation
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

x = 0 : 1 : 6;
for i = 1 : nclusters
    % plot cluster i
    fig = figure;title(sprintf('%dth Cluster, with %d genes',i,ngenes(i)),'FontSize',16);
    hold on;axis([0 6 -2 5])
    plot(x, plotData(idx==i,:),'Color','[.4,.4,.4]');
    plot(x,mean( plotData(idx==i,:) ),'Color','r','LineWidth',2);
    print(fig,sprintf('./Result/ncluster%d/fig%d',nclusters,i),'-dpng');
end
    

%% Talle Output
%% Output feature vector to a table
TFeaMtrx = array2table(FeaMtrx,'RowNames',name,'VariableNames',{'T0','T1','T2','T3','T4','T5','T6'});
writetable(TFeaMtrx,'./Result/TFeaMtrx.csv','WriteRowNames',true)
 
%% Output plot vector to a table
plotData = array2table(plotData,'RowNames',name,'VariableNames',{'T0','T1','T2','T3','T4','T5','T6'});
writetable(plotData,'./Result/plotData.csv','WriteRowNames',true)










