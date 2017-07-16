nclusters = 4;rng(1000); % For reproducibility
[ T_summary, Table_Result ] = f_kmeans( 'Profiles-ANan-DEGs.csv', nclusters );

myDir = sprintf('./Result/ncluster%d',nclusters); %gets directory
myFiles = dir(fullfile(myDir,'*.csv')); %gets all wav files in struct

fig = cell(length(myFiles),1);
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', myFiles(k).name);
  [ fig{k}, ~, ~,~, ~, ~ ] = ...
    f_plotTable( sprintf('%s/%s',myDir,myFiles(k).name), [], 'Mean Plot' );
end


