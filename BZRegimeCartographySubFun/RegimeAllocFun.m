%% Description %%
% This function allocates the regime membership base on the membership allocation function. 

function [RegimeMembership,RegimeIdx] = RegimeAllocFun(yCluster,RegimeCentre,RegimeList,m)

% Distance to center
DistanceToCenter = pdist2(yCluster,RegimeCentre,"euclidean");

a = DistanceToCenter.^(2/(m-1));
b = sum(1./(DistanceToCenter.^(2/(m-1))),2);

RegimeMembership = 1./(a.*b);
RegimeIdx = onehotdecode(RegimeMembership,RegimeList,2,"double");

end