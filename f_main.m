function [  ] = f_main( nc_array, pattern )
% @PGRP
% Nan
global myDir

if ~exist('pattern','var') || isempty(pattern)
    fprintf('Not clustering pattern is selected!\n');
    fprintf('Default using kmeans cluster the data into %s clusters\n',num2str(nc_array));
    pattern = 'KMEANS';
end

% nc_array -- a array of interger contains # of clusters you need;
len = length(nc_array);

for i = 1 : len
    nclusters = nc_array(i);
    myDir = sprintf('./%s_Result/ncluster%d',pattern,nclusters); %gets directory
    % fstr = 'Profiles-ANan-DEGs.csv';
    fstr = 'Profiles-ANan-up-regulated.csv';
    rng(1000); % For reproducibility seed
    fprintf('Applying kmeans to %s with number of cluster %d', fstr, nclusters);
    if strcmp(pattern,'KMEANS')
        [ ~, ~ ] = f_kmeans( fstr , nclusters );
    elseif strcmp(pattern,'HC')
            [~,~] = f_hc( fstr, nclusters );
    elseif strcmp(pattern,'FCMEANS')
        [~,~] = f_fcmeans(fstr, nclusters, []);
    end

    %% Plot Tables
    myFiles = dir(fullfile(myDir,'cluster*.csv')); %gets all wav files in struct

    fig = cell(length(myFiles),1);
    for k = 1:length(myFiles)
      baseFileName = myFiles(k).name;
      fullFileName = fullfile(myDir, baseFileName);
      fprintf(1, 'Now reading %s\n Preparing plot...\n', myFiles(k).name);
      [ fig{k}, ~, ~,~, ~, ~ ] = ...
        f_plotTable( sprintf('%s/%s',myDir,myFiles(k).name), [], 'Mean Plot' );
    end
end

end

