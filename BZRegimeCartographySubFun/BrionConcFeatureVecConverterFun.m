%% Description
% A function that converts Brion conentration time series data into a feature vector based on the PCA
% parameters and rescale coefficient created using "BrionConcFeatureVecGeneratorFun"

function [BrionConcFeatureVec,BrionConcFeature] = BrionConcFeatureVecConverterFun(BrionConcList,Eigenvec,muBrionConc, ...
                                                                 OscillateNoRescaleCoeff, ...
                                                                 OscillateFrequencyRescaleCoeff)

% Convert BrionConcList into PCA score
BrionConcPCAScore = (BrionConcList - muBrionConc)*Eigenvec;

% Identify number of peaks in Brion concentration time serie
OscillateNo = zeros(1,size(BrionConcList,1))';
OscillateTimeFrame = zeros(1,size(BrionConcList,1))';
OscillateFrequency = zeros(1,size(BrionConcList,1))';

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


OscillateNoRescaled = OscillateNo.*(OscillateNoRescaleCoeff(1)) + OscillateNoRescaleCoeff(2);
OscillateFrequencyRescaled = OscillateFrequency.*(OscillateFrequencyRescaleCoeff(1)) + OscillateFrequencyRescaleCoeff(2);

% Crete final feature vector
BrionConcFeatureVec = [BrionConcPCAScore,2.*OscillateNoRescaled,2.*OscillateFrequencyRescaled];
BrionConcFeature = [OscillateFrequency,OscillateNo];


end
