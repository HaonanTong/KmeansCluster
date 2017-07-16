function [ T_summary, Table_Result ] = f_kmeans( csv, nclusters )
% csv - exel file with 3 replicates for each time point, totally 6 time
% nclusters - # of clusters expected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@Haonan Tong
%PGRP Cluster in kmeans.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('#####################################\n');
fprintf('Kmeans algorithm on profiles from %s\n',csv);
fprintf('By Haonan Tong\n');
fprintf('@PGRP\n');
fprintf('#####################################\n');


%% Load RPKM Data
fprintf('loading file...\n');
T = readtable(csv,...
 'ReadVariableNames',true);

% summary(T);
 
Data = table2array(T(:,2:end));
agis = table2array(T(:,1));
% T.Properties.VariableNames

fprintf('Success!\n');

%% log-transformed RPKM values
fprintf('Transforming RPKM data...\n');
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

fprintf('taking log transformation...\n');
% log transformation
logData = log(Data);

%% Get feature vector
% Take the average transformed data over replicates for each treatment and
% each gene
FeaMtrx = [];

for i = 1:3:21%7 time points; 3 replicates;
   FeaMtrx = [FeaMtrx sum(logData(:,i:i+2),2)];
end

fprintf('Normalizing...\n');

% sample treatment means for each gene
FeaMtrx = 1/3*FeaMtrx;
FeaMtrx = FeaMtrx - repmat(mean(FeaMtrx,2),1,size(FeaMtrx,2));

fprintf('Success!\n');

fprintf('Applying kmeans...\n');
%% Kmeans
% Kmeans initialization
ngenes = zeros(nclusters,1); % initialize # of genes in each cluster;
str = sprintf('ncluster%d',nclusters);
mkdir('./Result',str);

% Kmeans operation
idx = kmeans(FeaMtrx,nclusters);


%% Result Analysis
% Print out # of genes for each cluster
for i = 1 : nclusters
    ngenes(i) = sum(idx==i);
end

T_summary = table((1:nclusters)',ngenes,'VariableNames',...
    {'cluster','ngenes'});

fprintf('Result Analyzing...\n');
fprintf('Summary:\n');

disp(T_summary);

fprintf('Success!\n');


% summary(T_summary);
fprintf('Exporing files...\n');
writetable(T_summary,sprintf('./Result/ncluster%d/T_summary',nclusters))

% Print out each cluster to files
Table_Result = cell(nclusters,1);
for i = 1 : nclusters
    Table_Result{i} = T(idx==i,:);
    writetable(Table_Result{i},sprintf('./Result/ncluster%d/cluster%d.csv',nclusters,i))
end

fprintf('Success!\n');
fprintf('Done!\n');
fprintf('#####################################\n');

%% Talle Output
%% Output feature vector to a table
% TFeaMtrx = array2table(FeaMtrx,'RowNames',name,'VariableNames',{'T0','T1','T2','T3','T4','T5','T6'});
% writetable(TFeaMtrx,'./Result/TFeaMtrx.csv','WriteRowNames',true)
 
%% Output plot vector to a table
% plotData = array2table(plotData,'RowNames',name,'VariableNames',{'T0','T1','T2','T3','T4','T5','T6'});
% writetable(plotData,'./Result/plotData.csv','WriteRowNames',true)

end


