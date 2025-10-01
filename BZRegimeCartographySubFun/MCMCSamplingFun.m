%% **Description**
% The MCMC sampling function performs a multichain metropolis MCMC sampling procedure on the multivariate 
% uncertainty distirbution generated from the neural network ensemble. The function takes in an ensemble of 
% neural network along with the limit of the sampling space boundaries and performs the sampling procedure 
% automatically. The convergence diagnostics used is the multivariate potential scale reduction factor (MPSRF) 
% which provides a guidance on when to stop the sampling procedure. 
%%
function [SamplePoints,SamplePointsUncertainty] = MCMCSamplingFun(ClassNNEns,nChain,SamplePerCycle,nCycleMax,nCycleMin,PropRangeRatio,T,xLimit,Folderpath)

%% MCMC parameters 

nDim = size(xLimit,1); % sampling space dimension

% Create starting point for MCMC chain
StartPointsRatio = lhsdesign(nChain,nDim,"iterations",10000);
StartPoints = StartPointsRatio.*repmat((xLimit(:,2) - xLimit(:,1))',nChain,1) + repmat(xLimit(:,1)',nChain,1);

% create file to save MCMC result
ResultFilepath = Folderpath + "MCMC sampling diagnostics";
if exist (ResultFilepath,"file") == 0
    mkdir(ResultFilepath)
end

%% Perform MCMC sampling

tStart = tic;
TargetPdf = @(x) TargetPdfFun(ClassNNEns,x,T); % function handler for target pdf

MCMCSample = MetropolisMCMC(TargetPdf,StartPoints,PropRangeRatio,nChain,nCycleMax,nCycleMin,SamplePerCycle,xLimit,ResultFilepath,tStart);

InfoCell = {"Input dimension no.:" + num2str(nDim);
            "Number of chains: " + num2str(nChain);
            "Maximum number of cycles: " + num2str(nCycleMax);
            "Sample per cycle: " + num2str(SamplePerCycle);
            "Ratio of probing space range which will be the range of uniform distribution: "+ num2str(PropRangeRatio);
            "Scaled gibbs distribution temperature (T): "+ num2str(T)};

writecell(InfoCell,ResultFilepath+"/mcmc.txt") % record MCMC param

% Unwrap and obtain uinique MCMC sampled points and calculate the uncertainty value
SamplePerChain = size(MCMCSample,1);
SamplePoints = zeros(SamplePerChain*nChain,nDim);

for iChain = 1:nChain
    SamplePoints((iChain-1)*(SamplePerChain)+1:iChain*(SamplePerChain),:) = MCMCSample(:,:,iChain);
end

SamplePoints = unique(SamplePoints,'rows');
SamplePointsUncertainty = UncertaintyQuantifyFun(ClassNNEns,SamplePoints);

save(ResultFilepath + "/MCMCSampledPoints.mat","SamplePoints","SamplePointsUncertainty"); % save reult

% Comparison between effectiveness of MCMC and LHS sampling in identifying highly uncertain points. 
LHSPointsRatio = lhsdesign(size(SamplePoints,1),size(xLimit,1),"iterations",100);
LHSPoints = LHSPointsRatio.*repmat((xLimit(:,2) - xLimit(:,1))',size(SamplePoints,1),1) + repmat(xLimit(:,1)',size(SamplePoints,1),1);
LHSPointsUncertainty = UncertaintyQuantifyFun(ClassNNEns,LHSPoints);

% Uncomment to visualize uncertainty of points sampled via LHS as oppose to MCMC sampling
% h2 = figure;
% boxplot([LHSPointsUncertainty,SamplePointsUncertainty],'Notch','on','Labels',{'LHS','MCMC'})
% ylabel("Uncertainty (Mutual Information)")
% ylim([0 1])
% grid minor
% title({'Comparison of uncertainty value','between MCMC and LHS sampling'})
% saveas(h2,ResultFilepath + "/MCMC effectiveness plot.jpg")
% close(h2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Important functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Target probability density function

% **Description**
% The target probability density function (target pdf) is the target distribution at which the MCMC sampler 
% will sample from. The target pdf = exp(uncertainty distribution/T) (A gibbs distribution). The main purpose 
% of converting the uncertainty distribution to a gibbs distribution is to prevent location with zero
% probabilities hence the MCMC sampler can jump from one mode to the other. Note that it is not necessary to
% scale the function into an actual pdf, hence strictly speaking the variable pdf is just a function 
% proportional to the actual pdf.


function pdf = TargetPdfFun(ClassNNEns,x,T)

pdf1 = UncertaintyQuantifyFun(ClassNNEns,x); % calculate uncertainty (mutual information) of x
pdf = exp(pdf1./T); % gibbs distribution

end

%% Metropolis MCMC function

% **Description**
% This version of the Metropolis MCMC function is design for the sampling of mcmc points for a predefined
% number of cycles. It runs multiple chains and performs multivariate potential scale reduction factor
% (MPRSF). Once the MPSRF reaches below a desired value (typically set at < 1.2), the sampling process stops
% and the function returns the sampled points.


function MCMCSample = MetropolisMCMC(TargetPdf,StartPoints,PropRangeRatio,nChain,nCycleMax,nCycleMin,SamplePerCycle,xLimit,ResultFilepath,tStart)

xProp = @(x) PropDistUni(x,PropRangeRatio,xLimit); % proposal distribution
Pdf = @(x) TargetPdf(x); % target distribution
MPSRFLimit = 1.2; % MPSRF threshold to stop sampling 
MCMCSample = [];
MPSRFList = [];
CumulativeTime = [];


for i = 1:nCycleMax

    MCMCSamplePerCycle = MetropolisAlgo(StartPoints,SamplePerCycle,Pdf,xProp,nChain); % run metropolis algorithm

    MCMCSample = cat(1,MCMCSample,MCMCSamplePerCycle);
    StartPoints = permute(MCMCSamplePerCycle(end,:,:),[3,2,1]); % final point of current cycle set as starting point for next cycle

    % MPSRF diagnostics
    MPSRF = MPSRFDiagnFun(MCMCSample,TargetPdf);
    MPSRFList = cat(1,MPSRFList,MPSRF);

    % Record of cumulative time taken for each cycle (in minutes)
    CumulativeTime = cat(1,CumulativeTime,toc(tStart)/60);
    
    % Save MCMC results in result folder
    save(ResultFilepath+"/MPSRF.mat","MPSRFList")
    save(ResultFilepath+"/MCMC_sample.mat","MCMCSample")
    save(ResultFilepath+"/Cumulative_time.mat","CumulativeTime")
   
    % Plot MPSRF vs iteration 
    h1 = figure;
    hold on
    plot(1:i,MPSRFList)
    title(["Multivariate potential scale reduction factor vs number of cycle" "Time taken = " + num2str(toc(tStart)/60) + "minutes"])
    ylabel("MPSRF")
    xlabel("Cycle no.")
    grid minor
    ylim([1,min(MPSRFList(1),MPSRFList(end).*2)])
    xlim([0,Inf])
    hold off
    saveas(h1,ResultFilepath + "/MPSR vs cycle no.jpg")
    close(h1)

    % Stop cycle when MPSRF less than MPSRFLimit AND when it runs for at least a certain number of cycles as
    % defined by nCycleMin
    if i > nCycleMin && MPSRF < MPSRFLimit
       break
    end
end


%% Proposal distribution function

% **Description**
% The proposal distribution function samples from the proposal distribution. the proposal distirbution is a
% uniform distribution where the minimum and maximum range of the dstribution will never exceed the boudaries
% of the sampling space boundaries. 

function xProp = PropDistUni(xLast,PropRangeRatio,xLimit)

HalfRange = (PropRangeRatio.*(xLimit(:,2) - xLimit(:,1)))/2;

UniUpperBound = min(xLast+HalfRange',xLimit(:,2)');
UniLowerBound = max(xLast-HalfRange',xLimit(:,1)');

xProp = (rand(size(xLast)).*(UniUpperBound-UniLowerBound)) + UniLowerBound;

end

end

%% Multivariate potential scale reduction factor diagnostic function

% **Description**
% The multivariate potential scale reduction factor (MPSRF) diagnostic function calculates the MPSRF from the
% sampled MCMC points. The MPSR formula follows the paper "General Methods for Monitoring Convergence of 
% Iterative Simulations" by Brooks and Gelman- section 3 and 4 - a general method for monitoring convergence 
% on an unknown distribution. The method uses the volume of a convex hull to identify the within chain and
% between chain dispersion (MPSRF = convex hull volume across all points/average convex hull volume of a 
% single chain). Only points within the 70% interval is selected for the volume calculation. When MPSRF
% falls  below 1.2, the distribution is assumed to converge.

function MPSRF = MPSRFDiagnFun(MCMCSample,TargetPdf)

% Remove final iteration for odd number samples
if mod(size(MCMCSample,1),2)~=0
    MCMCSample(end,:,:) = [];
end

[N,D,M]=size(MCMCSample); % N - number of points per chain, D - dimension, M - number of chains

% Discard the first half sequence of all chains to avoid burn-in period
MCMCSample = MCMCSample(N/2+1:end,:,:);

% Identify convex hull volume for each individual chain
VolumeChain = zeros(M,1);
for i  = 1:M

    MCMCSampleChain = MCMCSample(:,:,i);
    
    % Extract points which lies within the 70% interval
    y = TargetPdf(MCMCSampleChain);
    Qupper = quantile(reshape(y,[],1),0.85);
    Qlower = quantile(reshape(y,[],1),0.15);
    MCMCSampleChain = MCMCSampleChain(y<Qupper & y>Qlower,:);
    
    % display("pause now")
    % while i==1
    % end

    % Calculate convex hull volume
    [~,VolumeChain(i)] = convhulln(unique(MCMCSampleChain,'rows'));
end

VolumeChainMean = mean(VolumeChain,"all");

% Identify convex hull volume across all data points
MCMCSampleAll = zeros(M*(N/2),D);
for i  = 1:M
    MCMCSampleAll(((i-1)*(N/2)+1):i*(N/2),:) = MCMCSample(:,:,i);
end

% Extract points which lies within the 70% interval
y = TargetPdf(MCMCSampleAll);
Qupper = quantile(reshape(y,[],1),0.85);
Qlower = quantile(reshape(y,[],1),0.15);
MCMCSampleAll = MCMCSampleAll(y<Qupper & y>Qlower,:);

% Calculate convex hull volume for all pints
[~,VolumeAll] = convhulln(unique(MCMCSampleAll,'rows'));

MPSRF = VolumeAll/VolumeChainMean;

end

%% Metropolis algorithm

% **Description**
% The metropolis algorithm for MCMC sampling. 

function SamplePoints = MetropolisAlgo(StartPoints,nSamples,PdfFun,PropFun,nChain)

% check if proposal disibution working
PropPoint = PropFun(StartPoints);
if ~isfinite(PropPoint)
    error("Proposal distribution returns NaN value")
end

SamplePoints = zeros(nSamples,size(StartPoints,2),nChain);
 
CurrPoint = StartPoints;  %Current Point for each chain
RandList = rand(nChain,nSamples); % list of random number between 0 and 1

for i = 1:nSamples
    
    PropPoint = PropFun(CurrPoint); % Identify next set of proposed points based on current points
    AcceptRatio = PdfFun(PropPoint)./PdfFun(CurrPoint);
    AcceptRatio(isnan(AcceptRatio)) = 1;
    
    % Determine if to accept or reject the proposal point
    AccIdx = RandList(:,i) <= AcceptRatio;
    CurrPoint(AccIdx,:) = PropPoint(AccIdx,:);

    % Record new points
    SamplePoints(i,:,:) = permute(CurrPoint,[3,2,1]);
 
end

end

end
