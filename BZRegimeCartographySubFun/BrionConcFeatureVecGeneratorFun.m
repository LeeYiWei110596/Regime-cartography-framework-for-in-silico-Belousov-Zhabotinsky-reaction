

%% Description
% A function that converts Brion concentration time series data into a feature vector used for the initial 
% sampling stopping point and clustering. It performs feature extraction through concatenating the first n
% principal compoent score that explained > 99.5% of the spread of the data, the oscillation frequency and the 
% number of oscillations in the Brion concentration times series data into a single feature vector. The 
% frequency and number of oscilations are scaled so that it mataches the range of the first PCA score.

function [BrionConcFeatureVec,BrionConcFeature,Eigenvec,muBrionConc,OscillateNoRescaleCoeff,OscillateFrequencyRescaleCoeff] = BrionConcFeatureVecGeneratorFun(BrionConcList)

%% Perform PCA on Abs map data

% Perform pca
[coeff,score,~,~,explained,muBrionConc] = pca(BrionConcList);

% Identify the first n principal component where explained is > ExplainThres
Explain = 0;
nExplain = 0; % number of principal component where the explained value is just > ExplainThres
ExplainThres = 99.5;
for i = 1:size(explained,1)
    Explain = Explain + explained(i);
    if Explain > ExplainThres
        nExplain = i;
        break
    end
end

BrionConcPCAScore = score(:,1:nExplain);
Eigenvec = coeff(:,1:nExplain);

%% Identify the number of peaks and oscilating distance in each input

OscillateNo = zeros(size(BrionConcList,1),1);
OscillateTimeFrame = zeros(size(BrionConcList,1),1);
OscillateFrequency = zeros(size(BrionConcList,1),1);

for i = 1:size(BrionConcList,1)

    [pks,locs] = findpeaks(BrionConcList(i,:),"MinPeakWidth",10,"MinPeakHeight",-10);
    OscillateNo(i) = size(pks,2)-1; % not inlcuding the final, smaller peak

    if numel(locs)==1
        OscillateTimeFrame(i) = 0; 
        OscillateFrequency(i) = 0;
    else
        OscillateTimeFrame(i) = locs(end)-locs(1);
        OscillateFrequency(i) = (OscillateNo(i))./OscillateTimeFrame(i);
    end
    
end

max1stPC = max(BrionConcPCAScore(:,1),[],"all");
min1stPC = min(BrionConcPCAScore(:,1),[],"all");


% Rescale number of oscillations
minOscillateNo = min(OscillateNo,[],"all");
maxOscillateNo = max(OscillateNo,[],"all");
OscillateNoRescaled = ((OscillateNo - minOscillateNo)./(maxOscillateNo - minOscillateNo)).*(max1stPC - min1stPC) + min1stPC;
OscillateNoRescaleCoeff = [(max1stPC - min1stPC)/(maxOscillateNo - minOscillateNo), min1stPC-minOscillateNo*(max1stPC - min1stPC)/(maxOscillateNo - minOscillateNo)];

% Rescale frequence of oscillation
minOscillateFrequency = min(OscillateFrequency,[],"all");
maxOscillateFrequency = max(OscillateFrequency,[],"all");
OscillateFrequencyRescaled = ((OscillateFrequency - minOscillateFrequency)./(maxOscillateFrequency - minOscillateFrequency)).*(max1stPC - min1stPC) + min1stPC;
OscillateFrequencyRescaleCoeff = [(max1stPC - min1stPC)/(maxOscillateFrequency - minOscillateFrequency), min1stPC-minOscillateFrequency*(max1stPC - min1stPC)/(maxOscillateFrequency - minOscillateFrequency)];


BrionConcFeatureVec = [BrionConcPCAScore,2.*OscillateNoRescaled,2.*OscillateFrequencyRescaled]; 
BrionConcFeature = [OscillateFrequency,OscillateNo];

end
