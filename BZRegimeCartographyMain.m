clear all
close all
clc
% Add subfunction folder to path
addpath("BZRegimeCartographySubFun")
% Prime patternet for consistent result. Ignore as it does nothing for the algorithm
ClassNNEnsDummy = NNEnsTrainFun([1,1,1;2,2,2;3,3,3],[1,0,0;0,1,0;0,0,1],10,1,[]);
rng(0,"twister")

%% Description
% This is the main code for the in-silico BZ reaction regime cartography framework.

%% Regime cartography parameters

% Reactant concentration space parameter
ReactantConcSpaceDim          = 3;
ReactantConcSpaceAxisName     = {'Malonic Acid','Br{O_3}^{-1}','Ce^{3+}'}; % Axis name for plotting

% Input upper and lower concentration bound for reactant. Note that boundary obtain through trial and error
% where the ode equation does not return an error, and is based on the original concentration in the paper.

% ReactantConcUBLB = [Malonic acid lower bound, Malonic acid upper bound;
%                     BrO3- ion lower bound, BrO3- ion upper bound;
%                     Ce3+ ion lower bound, Ce3+ ion upper bound] - concentraton in M (mol dm-3)
ReactantConcUBLB = [1e-2,5e-2;1e-1,5e-1;1e-3,1e-2]; 

% Brion concentration plot parameters
tStep = 4000; % Number of time step for Brion concentration
tEnd  = 4000; % Brion concentration final time point (s)

% Initial sampling loop parameters
lambda        = 3;  % Lambda parameter for LoOP calculation
NNLoOP        = 5;  % Number of neighbours to consider during LoOP calculation
LoOPThreshold = 0.5;% LoOP threshold for the termination of initial sampling loop.
nInitialSampling  = 25; % Number of intial condition to sample from BZ model (before initial sampling loop) ("A" in paper)
nInitialSamplingLoopCond = 5; % Number of conditions to sample from BZ model at every initial sampling loop ("B" in paper) 

% Clustering parameters 
nRegime = 6;   % number of cluster to divide. ("C" in paper)
m       = 1.5; % m value for membership function

% Active learning parameters %
nDesignPerIter = 5;   % total number of design point at each active sampling iteration ("P" in paper)  
nChain         = 25;  % number of MCMC chains which runs in parallel                         
SamplePerCycle = 400; % number of MCMC samples per cycle                              
nCycleMax      = 100; % maximum number of MCMC cycles, set as cut off to prevent cycles to run indefnitely
nCycleMin      = 1;   % minimum number of MCMC cycles
T              = 0.2; % Temperature of gibbs distribution                                   
PropRangeRatio = 0.3; % ratio of range of probing space that will be the range of uniform proposal distribution ("r" in paper)
nEnsemble      = 100; % number of neural network ensembles ("E" in paper)
NNStructure    = [25 25 25]; % structure of neural network ensemble
nActiveLearningIter = 4; % number of active learning iteration

%% Create folders and add relevant folder to path

% Create result folder 
InitialSamplingResultFolderPath = "Result/Initial Sampling/";
if ~exist(InitialSamplingResultFolderPath, 'dir')
   mkdir(InitialSamplingResultFolderPath)
end

SoftClusteringResultFolderPath = "Result/Soft Clustering/";
if ~exist(SoftClusteringResultFolderPath, 'dir')
   mkdir(SoftClusteringResultFolderPath)
end

TestDataResultFolderPath = "Result/Test Data/";
if ~exist(TestDataResultFolderPath, 'dir')
   mkdir(TestDataResultFolderPath)
end

ActiveLearningResultFolderPath = "Result/Active Learning/";
if ~exist(ActiveLearningResultFolderPath, 'dir')
   mkdir(ActiveLearningResultFolderPath)
end

%% Initial Sampling

% Sample initial data using FSSF-fr algorithm

% Note that for the initial concentration range for each reactant are very different, hence we perform FSSF-fr 
% sampling from a normalized space where each axis has a range of 0 to 1 and rescale the sampled points to the 
% appropriate values.

% List of samples from the normalized space
ReactantConcNormList = FSSFfrSamplingFun(nInitialSampling,[0,1;0,1;0,1],[]);
% List of reactant concentration sampled
ReactantConcList = ReactantConcNormList.*(ReactantConcUBLB(:,2) - ReactantConcUBLB(:,1))' + ReactantConcUBLB(:,1)';
% List of Brion concentration sampled
BrionConcList = InSilicoBZReaction(ReactantConcList);
% Convert Brion concentration to feature vector
[BrionConcFeatureVecList,BrionConcFeatureList,Eigenvec,muBrionConc,OscillateNoRescaledCoeff,OscillateFrequencyRescaledCoeff]... 
= BrionConcFeatureVecGeneratorFun(BrionConcList);
CondList = (1:nInitialSampling)'; % List of conditions

EndInitialSamplingLoop = false; % Marker to end initial sampling loop
InitialSamplingLoopCounter = 1;
LoOPRec = [];

% Perform initial sampling loop 
while EndInitialSamplingLoop == false    

    % Perfrom sampling on the input space using FSSFfrSamplingFun
    ReactantConcNormAdd = FSSFfrSamplingFun(nInitialSamplingLoopCond,[0,1;0,1;0,1],ReactantConcNormList);
    ReactantConcAdd = ReactantConcNormAdd.*(ReactantConcUBLB(:,2) - ReactantConcUBLB(:,1))' + ReactantConcUBLB(:,1)';
    BrionConcAdd    = InSilicoBZReaction(ReactantConcAdd); % Run addtional sample through in-silico BZ reaction
    CondAdd         = (CondList(end)+1:CondList(end)+nInitialSamplingLoopCond)';

    % Add sampled points to list
    ReactantConcNormList = [ReactantConcNormList;ReactantConcNormAdd];
    ReactantConcList = ReactantConcNormList.*(ReactantConcUBLB(:,2) - ReactantConcUBLB(:,1))' + ReactantConcUBLB(:,1)';
    BrionConcList = [BrionConcList;BrionConcAdd];
    [BrionConcFeatureVecList,BrionConcFeatureList,Eigenvec,muBrionConc,OscillateNoRescaledCoeff,OscillateFrequencyRescaledCoeff]...
    = BrionConcFeatureVecGeneratorFun(BrionConcList);
    CondList = [CondList;CondAdd];
    
    % Perform LoOP calculation for newly sampled point
    LoOP = LoOPFun(BrionConcFeatureVecList,CondAdd,NNLoOP,lambda);
    LoOPRec = [LoOPRec,LoOP];

    if sum(LoOP>LoOPThreshold) == 0
        EndInitialSamplingLoop = true;
    else
        InitialSamplingLoopCounter = InitialSamplingLoopCounter+1;
    end

end

% Save Flowrate, Brion concentration, Brion concentration feature vector and Brion concentration feature
ReactantConcTable = array2table([CondList,ReactantConcList],'VariableNames',{'Cond','[MA]_0 (M)','BrO3- (M)','Ce3+ (M)'});
writetable(ReactantConcTable,InitialSamplingResultFolderPath + "ReactantConc0.xlsx")
writematrix(BrionConcList,InitialSamplingResultFolderPath + "BrionConc0.csv")
writematrix(BrionConcFeatureVecList,InitialSamplingResultFolderPath + "BrionConcFeatureVec0.csv")
BrionConcFeatureListTable = array2table([CondList,BrionConcFeatureList]);
BrionConcFeatureListTable.Properties.VariableNames = ["Cond","n_osc","f_osc"];
writetable(BrionConcFeatureListTable,InitialSamplingResultFolderPath + "BrionConcFeature0.csv")

% Save LoOP records
save(InitialSamplingResultFolderPath + "LoOP record.mat","LoOPRec","LoOPThreshold","NNLoOP","lambda")

%% Soft clustering

% Perform spectral clustering 
D = pdist(BrionConcFeatureVecList);
DistanceMat = squareform(D);
% Estimate kernal scale using the median of the nearest neighbour of each Abs map. 
NNDist = max(mink(DistanceMat,2,2),[],2);
KernalScale = median(NNDist,"all");   % "\sigma" in paper
S = exp(-(DistanceMat.^2)./(2*KernalScale^2)); % create similarity matrix
[~,~,EigenvalueList] = spectralcluster(S,10,'Distance','precomputed','NumNeighbors',CondList(end));

h1 = figure;
hold on
title("Spectral clustering eigenvalue plot")
plot(1:10,EigenvalueList,".")
xlabel("Eigenvalue no.")
ylabel("Eigenvalue")
grid minor
hold off
saveas(h1,SoftClusteringResultFolderPath +"Spectral clustering eigenvalue plot.fig")

if isempty(nRegime)
    % Take in number of clusters from user is number of clusters not set
    nRegime = input("Number of clusters for C means clustering = ");
end

% Perform membership function allocation

% Calculate regime centroid
RegimeSpectralIdx = spectralcluster(S,nRegime,'Distance','precomputed','NumNeighbors',CondList(end));
RegimeList = unique(RegimeSpectralIdx);
RegimeCentroid = zeros(nRegime,size(BrionConcFeatureVecList,2));

for i = 1:nRegime
    BrionConcFeatureVecRegime = BrionConcFeatureVecList(RegimeSpectralIdx==i,:);
    RegimeCentroid(i,:) = mean(BrionConcFeatureVecRegime,1);
end

% Perform fuzzy based membership function allocation based on regime centroid
DistanceToCenter = pdist2(BrionConcFeatureVecList,RegimeCentroid,"euclidean");
a = DistanceToCenter.^(-2/(m-1));
b = sum((DistanceToCenter.^(-2/(m-1))),2);
RegimeMembershipList = a./b;
[~,RegimeIdxList] = max(RegimeMembershipList,[],2);

% Plot clustering result for each regime
t0Probe = linspace(0,tEnd,tStep);
for i = 1:nRegime
    
    % Identify the number of members within the regime
    RegimeMemberNum = sum(RegimeIdxList==i,"all");
    RegimeMembershipPlot = RegimeMembershipList(RegimeIdxList==i,i);
    CondListRegime = CondList(RegimeIdxList==i,:);

    % Identify Brion time series within regime
    BrionConcListRegime = BrionConcList(RegimeIdxList==i,:);
    
    % Identify number of rows in subplots
    nColumnSubplot = 5; % number of columns in subplot
    nRowsSubplot = ceil(RegimeMemberNum/nColumnSubplot);
    h2 = figure;
    hold on
    % Create subplot
    for j = 1:RegimeMemberNum
        subplot(nRowsSubplot,nColumnSubplot,j)
        plot(t0Probe,BrionConcListRegime(j,:));
        xlim([0,tEnd])
        ylim([-15 -4])
        xlabel("time(s)")
        ylabel("log([Br^{-1}])")
        title(["Condition " + num2str(CondListRegime(j)),"Membership value = "+num2str(round(RegimeMembershipPlot(j),3))])
        grid  minor
    end
    sgtitle("Regime "+num2str(i),fontsize=15)
    set(gcf,'position',[100,100,300*nColumnSubplot,200*nRowsSubplot])
    hold off

    saveas(h2,SoftClusteringResultFolderPath + "Regime"+num2str(i)+".jpg")
    
end

% Extract most representative condition in each regime and plot result
[MostRepresentativeConditionMembership,MostRepresentativeCondition] = max(RegimeMembershipList,[],1);
h3 = figure;
% Create subplot
nColumnSubplot = 5; % number of columns in subplot
nRowsSubplot = ceil(nRegime/nColumnSubplot);
hold on
for j = 1:nRegime
    subplot(nRowsSubplot,nColumnSubplot,j)
    plot(t0Probe,BrionConcList(MostRepresentativeCondition(j),:));
    xlim([0,tEnd])
    ylim([-15 -4])
    xlabel("time(s)")
    ylabel("log([Br^{-1}])")
    title(["Regime "+j,"Condition " + num2str(MostRepresentativeCondition(j))])
    grid  minor
end
sgtitle("Most representative regime ",fontsize=15)
set(gcf,'position',[100,100,300*nColumnSubplot,250*nRowsSubplot])
hold off
saveas(h3,SoftClusteringResultFolderPath + "Most representative condition.jpg")

% Save clustering result
save(SoftClusteringResultFolderPath + "Regime Clustering Record.mat",'RegimeIdxList','RegimeMembershipList','RegimeCentroid')


%% Set up test data
% Generate TestData to identify performance of ensemble neural network (for the termination of active sampling loop)

[ReactantConcTest,RegimeMembershipTest] = BZTestDataCreateFun(@InSilicoBZReaction, ...
                                                              @BrionConcFeatureVecConverterFun, ...
                                                              RegimeCentroid, RegimeList,m, ...
                                                              ReactantConcUBLB,Eigenvec,muBrionConc, ...
                                                              OscillateNoRescaledCoeff, ...
                                                              OscillateFrequencyRescaledCoeff, ...
                                                              ReactantConcSpaceAxisName, ...
                                                              TestDataResultFolderPath);


%% Active learning

ActiveLearningIter = 0; % Counter for active learning loop, start from iteration 0 
EndActiveLearningLoop = false; % Marker to end initial sampling loop    
NNEnsPerformanceRec = []; % set empty matrix for performance record

while EndActiveLearningLoop == false

    % Create subfolder in active learning folder
    ActiveLearningResultSubFolderPath = ActiveLearningResultFolderPath + "Iter"+num2str(ActiveLearningIter + "/");
    if not(exist(ActiveLearningResultSubFolderPath,'dir'))
            mkdir(ActiveLearningResultSubFolderPath)
    end

    if ActiveLearningIter == 0
        
        % Write data extracted from initial sampling into active learning iter 0 subfolder.
        writetable(ReactantConcTable,ActiveLearningResultSubFolderPath + "ReactantConc0.xlsx")
        writematrix(BrionConcList,ActiveLearningResultSubFolderPath + "BrionConc0.csv")
        writematrix(BrionConcFeatureVecList,ActiveLearningResultSubFolderPath + "BrionConcFeatureVec0.csv")
        writetable(BrionConcFeatureListTable,ActiveLearningResultSubFolderPath + "BrionConcFeature0.csv")

    else
        % Perform MCMC sampling
        [CandidatePoints,CandidatePointsUncertainty] = MCMCSamplingFun(ClassNNEns,nChain,SamplePerCycle, ...
                                                                       nCycleMax,nCycleMin,PropRangeRatio, ...
                                                                       T,ReactantConcUBLB, ...
                                                                       ActiveLearningResultSubFolderPath);

        % Suggest feed concentration condition to sample through kmeans stratfied sampling
        ReactantConcIter = ExperimentDesignFun(CandidatePoints,CandidatePointsUncertainty,nDesignPerIter);
        ReactantConcIterNorm = (ReactantConcIter - ReactantConcUBLB(:,1)')./(ReactantConcUBLB(:,2) - ReactantConcUBLB(:,1))';

        % Run suggested condition through in-silico BZ reaction
        CondIter = (CondList(end)+1:CondList(end)+nDesignPerIter)';
        BrionConcIter =  InSilicoBZReaction(ReactantConcIter);
        [BrionConcIterFeatureVec,BrionConcIterFeature] = BrionConcFeatureVecConverterFun(BrionConcIter, ...
                                                                                         Eigenvec,muBrionConc, ...
                                                                                         OscillateNoRescaledCoeff, ...
                                                                                         OscillateFrequencyRescaledCoeff);
       
        % Identify regime membership for newly sampled condition
        [RegimeMembershipIter,RegimeIdxIter] = RegimeAllocFun(BrionConcIterFeatureVec,RegimeCentroid,RegimeList,m);
        
        % Concatenate iter condition with previously sampled conditions
        CondList                = [CondList;CondIter];
        ReactantConcList        = [ReactantConcList;ReactantConcIter];
        ReactantConcNormList    = [ReactantConcNormList;ReactantConcIterNorm];
        BrionConcList           = [BrionConcList;BrionConcIter];
        BrionConcFeatureVecList = [BrionConcFeatureVecList;BrionConcIterFeatureVec];
        BrionConcFeatureList    = [BrionConcFeatureList;BrionConcIterFeature];
        RegimeMembershipList    = [RegimeMembershipList;RegimeMembershipIter];
        RegimeIdxList           = [RegimeIdxList;RegimeIdxIter];

        % Create and save table cotainining reactant concentration, regime index and regime membership 
        % collected during current iteration
        IterRegimeMembershipTable = array2table(RegimeMembershipIter);
        IterRegimeIdxTable = array2table(RegimeIdxIter,"VariableNames","Regime");
        IterReactantTable = array2table([CondIter,ReactantConcIter]);
        IterReactantTable.Properties.VariableNames = ["Cond",string(ReactantConcSpaceAxisName)];
        ActiveLearningIterTable = [IterReactantTable,IterRegimeIdxTable,IterRegimeMembershipTable];
        writetable(IterReactantTable, ActiveLearningResultSubFolderPath + "ReactantConcIter.xlsx")

        % Save Brion conc, Brion conc feature vector, and Brion conc feature for this iter
        writematrix(BrionConcIter,ActiveLearningResultSubFolderPath + "BrionConcIter.csv")
        writematrix(BrionConcIterFeatureVec,ActiveLearningResultSubFolderPath + "BrionConcFeatureVecIter.csv")
        BrionConcIterFeatureTable = array2table([CondIter,BrionConcIterFeature]);
        BrionConcIterFeatureTable.Properties.VariableNames = ["Cond","n_osc","f_osc"];
        writetable(BrionConcIterFeatureTable,ActiveLearningResultSubFolderPath + "BrionConcIterFeature.csv")

    end

    % Train neural network ensemble
    ClassNNEns = NNEnsTrainFun(ReactantConcList,RegimeMembershipList,NNStructure,nEnsemble,ActiveLearningResultSubFolderPath);
    % Plot regime
    RegimeVisualizationFun(ClassNNEns,RegimeList,ReactantConcUBLB,ReactantConcSpaceAxisName,ActiveLearningResultSubFolderPath)
    close all
    % Identify performance of neural network and decide whether to end active learning loop
    NNEnsPerformanceRec = NNEnsPerformanceDiag(ClassNNEns,ReactantConcTest,RegimeMembershipTest,NNEnsPerformanceRec,ActiveLearningResultFolderPath);
    if ActiveLearningIter == nActiveLearningIter
        break
    else
        ActiveLearningIter = ActiveLearningIter + 1;
    end

end

%% Save all data

% Save all condition, membership and concentration
RegimeMembershipListTable = array2table(RegimeMembershipList);
RegimeMembershipListTable.Properties.VariableNames = ["Regime 1 Membership","Regime 2 Membership",...
                                                      "Regime 3 Membership","Regime 4 Membership",...
                                                      "Regime 5 Membership","Regime 6 Membership"];
RegimeIdxTable = array2table(RegimeIdxList,"VariableNames","Regime");
ReactantConcTable = array2table([CondList,ReactantConcList]);
ReactantConcTable.Properties.VariableNames = ["Cond",string(ReactantConcSpaceAxisName)];
FeatureTable = array2table(BrionConcFeatureList);
FeatureTable.Properties.VariableNames = ["f_osc","n_osc"];
AllCondDataTable = [ReactantConcTable,RegimeIdxTable,RegimeMembershipListTable,FeatureTable];
writetable(AllCondDataTable, "Result/AllConditionRec.xlsx")

% Save all condition's Brion concentration and Brion feature vector
writematrix(BrionConcList, "Result/AllConditionBrionConc.csv")
writematrix(BrionConcFeatureVecList, "Result/AllConditionBrionConcFeatureVec.csv")

close all



