%% Description %%
% This function performs the Local Outlier Probability (LoOP) using the formula proposed by Kriegel et.al 
% (2009) (doi: https://doi.org/10.1145/1645953.1646195). The function takes in a set of data and calculated
% the LoOP for specfic points within the dataset. 

%%
function LoOP = LoOPFun(y,yIdxProbe,NN,lambda)

LoOP = zeros(size(yIdxProbe,1),1);
yBase = y(1:yIdxProbe(1)-1,:);
yProbe = y(yIdxProbe,:);

for i = 1:size(yProbe,1)

    yi = [yBase;yProbe(i,:)];

    % Calculate k nearest neighbour distance for all abs map samples
    [NNIdx,distance] = knnsearch(yi,yi,"K",NN+1);
    distance = distance(:,2:end);
    NNIdx = NNIdx(:,2:end);

    % Calculate pDist for each point
    pDist = lambda*((sum(distance.^2,2)./(NN)).^0.5);

    % Calculate PLOF for each sample point
    PLOF = (pDist./mean(pDist(NNIdx),2))-1;
    nPlOF = lambda*((mean(PLOF.^2,"all"))^0.5);

    % Calculate outlier probability of additional sample
    LoOP(i) = max(0,erf(PLOF(end)/(nPlOF*(2^0.5))));

end


end
