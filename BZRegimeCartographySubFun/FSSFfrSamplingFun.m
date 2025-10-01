
%% Description %%
% This function samples points from the input space based on the FSSF-fr paper by Shang, B., & Apley, D. W. 
% (2020) (doi: 10.1080/00224065.2019.1705207). The code can suggest sampling points with or without taking 
% into consideration the points that has been sampled so far.

%%
function DesignedPoints = FSSFfrSamplingFun(nSample,xLimit,SampledPoints)

%% Main codes

nDim = size(xLimit,1); % Calculate number of dimenisons
SpaceMidPoint = sum(xLimit,2)./2; % Midpoint of the sampling space

% Use latin hypercube to generate candidate points spreaded throughout the sampling space
% Note: for the number of candidate points proposed by LHS, Shang, B., & Apley, D. W. (2020) proposed the
% formula 1000*nDim+2*n where n is the number of points to sample. For our case, we set this number to 10 000
% which way exceeds the proposed number, for simplicity. 
rng(0,"twister") % Set random number generator to ensure similar LHS output
CandidatePointsLHS = lhsdesign(10000,nDim,"iterations",1000); % Use the formula suggested by Shang, B., & Apley, D. W. (2020) to identify number of candidate points
CandidatePoints = CandidatePointsLHS.*(xLimit(:,2) - xLimit(:,1))' + xLimit(:,1)';

% Calculate distance of candidate points to imaginary reflected points at boundaries. 
DistReflectLBLHS = CandidatePoints-(xLimit(:,1))'; % Distance between LHS candidate points and the lower boundaries.
DistReflectUBLHS =  (xLimit(:,2))' - CandidatePoints; % Distance between LHS candidate points and the upper boundaries.
DistReflectLBLHS(DistReflectLBLHS>(SpaceMidPoint'-(xLimit(:,1))')) = 0; 
DistReflectUBLHS(DistReflectUBLHS>(SpaceMidPoint'-(xLimit(:,1))')) = 0;
DistReflectLHS = 2.*((2*nDim).^0.5).*(DistReflectLBLHS + DistReflectUBLHS); % Keep the shortest distance between the two boundaries

% Identify design points sequentially. After a point is selected as the design point, that point is included 
% within the lists of experimentally sampled points while calculating the minimum distance during the 
% selection of the next design point.This is to ensure that the next design point does not lie close to the 
% current design point.

DesignedPoints = []; % Points that have been selected as the design points for exploration

% Select exploration points sequentially
for i = 1:nSample

    % Calculate minimum distance
    if i ==1
        if ~isempty(SampledPoints)
            % min Dist for sampling with sampled points
            minDistSampled = min(pdist2(CandidatePoints,SampledPoints),[],2); % Minimum distance between candidate points and sampled points
            minDist = min([minDistSampled,DistReflectLHS],[],2); % Minimum distance between each candidate points to sampled points + reflected imaginary points
        else
            % min Dist for blank sampling
            minDist = min(DistReflectLHS,[],2); % Minimum distance between each candidate points to reflected imaginary points
        end
    else
        DistToDesign = pdist2(CandidatePoints,DesignedPoints(end,:)); % Distance between candidate points and newly identified design point
        minDist = min(minDist,DistToDesign); % Update minimum distance record based on distance to newly choosen design point
    end

    % Identify design point
    [~,idxCandidate] = max(minDist,[],1);
    DesignedPoints  = cat(1,DesignedPoints,CandidatePoints(idxCandidate,:));

end

end

