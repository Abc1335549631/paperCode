%% parameters set
dbstop if error
tic
addpath(genpath(pwd))
input.openFile ='D4.mat';
input.clusterSize = 15;
input.alpha = 0.31736;     
input.minS = 59;
input.R1 = 0.48180;        
input.lMax = 46;
input.minSm = 0.03319;     
input.lpH = 15;
input.similarBnd = 0.1;
%% end parameters

hca (input);
toc
rmpath(genpath(pwd));