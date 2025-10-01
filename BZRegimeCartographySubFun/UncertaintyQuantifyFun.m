%% Description
% The uncertainty quantification function utilizes the mutual information of the classification neural network 
% ensemble (ClassNNEns) to quantify the uncertainty for a particular input. Refer to Y. Gal, "Uncertainty in 
% deep learning" phd thesis. 

%%
function Uncertainty = UncertaintyQuantifyFun(ClassNNEns,x)

% Obtain probability score for each points within x
Score = NNEnsPredict(ClassNNEns,x);
ScoreAvg = mean(Score,3,"omitnan");

% General mutual information method, comment out if not used
Score = Score+eps; % add eps so 0*log2(0) returns a very small value close to zero instead of NaN
PredEntropy = -sum(ScoreAvg.*log2(ScoreAvg),2);
ExpEntropy  = -mean(sum(Score.*log2(Score),2),3,"omitnan");

% Calculate the final uncertainty value using mutual information formula
Uncertainty = PredEntropy - ExpEntropy;
Uncertainty(Uncertainty<= 0) = 0; % to ensure no negative values within the uncertainty vector

end