function [ T_summary, Table_Result ] = f_fcmeans( csv, nclusters, threshold )
% csv - exel file with 3 replicates for each time point, totally 6 time
% nclusters - # of clusters expected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@Haonan Tong
%PGRP Hierarchical Clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% threshold for filtering Fuzzy Partition Matrix a nclusters by 1 matrix
% s.t. threshold each clusters.

global myDir


warning('off')

fprintf('#####################################\n');
fprintf('Fuzzy c means algorithm on profiles from %s\n',csv);
fprintf('By Haonan Tong\n');
fprintf('@PGRP\n');
fprintf('#####################################\n');


%% Load RPKM Data
fprintf('loading file...\n');
T = readtable(csv,...
 'ReadVariableNames',true,'ReadRowNames',true);

% summary(T);
 
Data = table2array(T);
[n_gene,~] = size(Data);
fprintf('Totally %d genes\n',n_gene);
% T.Properties.VariableNames

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

fprintf('Applying fuzzy c means...\n');
%% Fuzzy c-means
% initialization
ngenes = zeros(nclusters,1); % initialize # of genes in each cluster;
str = sprintf('ncluster%d',nclusters);
mkdir('./FCMEANS_Result',str);

% fcmeans operation
[~, U] = fcm(FeaMtrx,nclusters);

% find the idx indicator
idx = zeros(n_gene,1);

Ufilter = U .* ( U == repmat(max(U),nclusters,1) );

if isempty(threshold)
    threshold = zeros(nclusters,1);
end

thrmatrix = repmat(threshold, 1, n_gene);
DecisionMatrix = (Ufilter - thrmatrix) > 0; % DecisionMatrix(i,:) contains gene index for cluster i;

for i = 1 : n_gene
    if any(DecisionMatrix(:,i))
        tmp = find(DecisionMatrix(:,i)==1);
        idx(i) = tmp(1);
    end
end

VarNames = cell(nclusters,1);
for i = 1 : nclusters
    VarNames{i} = sprintf('score%d',i);
end
T_cscore = array2table(U','RowNames',T.Properties.RowNames,'VariableNames',VarNames);
writetable(T_cscore,sprintf('%s/T_cscore.csv',myDir),'WriteRowNames',true,'WriteVariableNames',true);


%% Result Analysis
% Print out # of the different spliced versions of a gene
for i = 1 : nclusters
    ngenes(i) = sum(idx==i);
end

T_summary = table((1:nclusters)',ngenes,'VariableNames',...
    {'cluster','ngenes'});

fprintf('Result Analyzing...\n');
fprintf('Summary of different spliced versions of a gene(eg. AT1G17040.1)\n');

disp(T_summary);

% summary(T_summary);
writetable(T_summary,sprintf('%s/T_summary.txt',myDir))

% Print out each cluster to files
ngenes2 = zeros(nclusters,1);

Table_Result = cell(nclusters,1);
for i = 1 : nclusters % derive subtable for each cluster
    Table_Result{i} = T(idx==i,:);
    writetable(Table_Result{i},sprintf('%s/cluster%d.csv',myDir,i),'WriteRowNames',true)
end

fprintf('Unique annotation of gene in TAIR ID...\n');
c_cluster = cell(nclusters,1);
for i = 1 : nclusters % print out anonotation for genes
    fp = fopen(sprintf('%s/Gene-list-cluster%d.txt',myDir,i),'wt');
    fprintf(fp, 'cluster%d\n', i);
    tmp = strtok(Table_Result{i}.Properties.RowNames,'.');
    c_cluster{i} = unique(tmp);
    ngenes2(i) = length(c_cluster{i});
    fprintf(fp, '%s\n', c_cluster{i}{:,:});
    fclose(fp);
end

MainFolder = pwd;
cd(myDir)
unix(sprintf('paste -d ''\t'' Gene-list-cluster*.txt > T_validation%d.txt',nclusters))
cd(MainFolder);

% Print out # of genes for each cluster

T_summary2 = table((1:nclusters)',ngenes2,'VariableNames',...
    {'cluster','ngenes'});

fprintf('Summary of # of gene(eg. AT1G17040)\n');

disp(T_summary2);


% summary(T_summary);
writetable(T_summary2,sprintf('%s/T_summary2.txt',myDir))


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







