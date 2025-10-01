%% Description %%
% This function selects the next set of points to select from experiment platform from the candidate points 
% based on candidate points uncertainty. The algorithm below is similar to a K-means stratified sampling, but 
% instead of sampling from each strata, we choose the candidate point that is the closest to the centroid of
% the K means cluster. 

%%
function DesignPoints = ExperimentDesignFun(CandidatePoints,CandidatePointsUncertainty,nDesignTotal)

% Identify the 75% quantile of the uncertainty values as cutoff
UncertaintyCutOff = quantile(CandidatePointsUncertainty,0.75);
CandidatePoints = CandidatePoints(CandidatePointsUncertainty>UncertaintyCutOff,:);

% Perform K means clustering on the candidate exploit points to identify strata. Number of clusters depends on
% the number of design points to sample next in the experimental platform
[~,Centroid] = kmeans(CandidatePoints,nDesignTotal,"Replicates",10,"MaxIter",10000);

% The design point is the closest candiate points to each centroid
[~,I]= pdist2(CandidatePoints,Centroid,"euclidean","Smallest",1);
DesignPoints = CandidatePoints(I,:);

end