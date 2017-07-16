%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@Haonan Tong
%PGRP Cluster in kmeans.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Have Differential Expressed genes 1464;
T = readtable('Arbdpss_Dffrntl_xprss_gns.csv',...
    'ReadVariableNames',true);
summary(T);
Data = table2array(T(:,2:end));
name = table2array(T(:,1));

newData = [];
%
for i = 1:3:21%7 time points; 3 replicates;
   newData = [newData sum(Data(:,i:i+2),2)];
end
newData = 1/3*newData;

tmp = [];
for i = 2:7
   tmp = [tmp log( newData(:,i)./newData(:,i-1) )];
end
clear Data;
Data = tmp;

%
nclusters = 10; ngenes = zeros(nclusters,1);
str = sprintf('ncluster%d',nclusters);
mkdir('./Result',str);

idx = kmeans(Data,nclusters);
for i = 1 : nclusters
    ngenes(i) = sum(idx==i);
end
display(ngenes);

table = cell(nclusters,1);
for i = 1 : nclusters
    table{i} = T(idx==i,:);
    writetable(table{i},sprintf('./Result/ncluster%d/cluster%d.csv',nclusters,i))
end


%%



