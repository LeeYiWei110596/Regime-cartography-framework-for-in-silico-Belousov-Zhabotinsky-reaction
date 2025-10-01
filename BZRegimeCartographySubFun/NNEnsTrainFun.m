%% Description
% The neural netwotrk ensemble training function trains an ensemble of classification neural network based on 
% the user defined neural network strucutre and number of members within an ensemble.

function ClassNNEns = NNEnsTrainFun(xTrain,yTrainScore,nNodes,nEnsemble,FolderPath)


%% Create neural network ensemble
ClassNNEns = cell(nEnsemble,1);
TrainRecNNEns = cell(nEnsemble,1);

for i = 1:nEnsemble
    [Mdl,TrainRec] = NNTrain(nNodes,xTrain,yTrainScore);
    ClassNNEns(i,1) = {Mdl};
    TrainRecNNEns(i,1) = {TrainRec};
end

if ~isempty(FolderPath)
    save(FolderPath + "ClassNNEns.mat","ClassNNEns") % save neural network ensemble
    save(FolderPath + "ClassNNEnsTrainingRecord.mat","TrainRecNNEns") % save the training record of the ensemble
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Important function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Neural network training function

% **Description**
% The neural netwotrk training function trains a single classification neural network based on the user defined 
% neural network strucutre. Neural network initialization has been changed to Glorot initialization ,as it is
% one of the most common ways people initialize neural network ensembles (as oppose to the default Nguyen-Widrow 
% initialization used in Matlab's original code.


function [Mdl,TrainRec] = NNTrain(nNodes,xTrain,yTrainScore)

% Neural network training setting
net = patternnet(nNodes,'trainrp'); % set training algorithm - use trainrp

%% comment out for normal random initial sampling
% Configure network
net = configure(net,xTrain',yTrainScore'); 

% Set neural network initial weight through glorot initialization
nnStruct = [size(xTrain,2),nNodes,size(yTrainScore,2)];
nHiddenLayer = size(nNodes,2)+1;

% Set bias to zero
for layerNum = 1:nHiddenLayer 
    net.b{layerNum,1} = zeros(nnStruct(layerNum+1),1);
end

% Set weight of input layer through glorot initialization
net.IW{1,1} = Glorot([nnStruct(2),nnStruct(1)],nnStruct(1),nnStruct(2));

% Set weight of hidden layer through glorot initialization
for layerNum = 1:nHiddenLayer-1
    net.LW{layerNum+1,layerNum} = Glorot([nnStruct(layerNum+2),nnStruct(layerNum+1)],nnStruct(layerNum+1),nnStruct(layerNum+2));
end
%

net.performFcn = 'crossentropy'; 
net.performParam.regularization = 0;
net.trainParam.epochs = 10000; % maximum training epoch
net.trainParam.showWindow = false;  % surpress GUI
net.trainParam.showCommandLine = false;
net.trainParam.show = 1000;
net.trainParam.goal = 0;
net.trainParam.min_grad = 1e-4;
net.divideParam.trainRatio = 1; % train, validation, test ratio division. Always set to 1:0:0
net.divideParam.valRatio   = 0;
net.divideParam.testRatio  = 0;
net.trainParam.lr = 0.01;
[Mdl,TrainRec] = train(net,xTrain', yTrainScore',"useGPU","no"); 

function weights = Glorot(sz,numOut,numIn)

Z = 2*rand(sz) - 1; % Change range of uniform distribution to [-1 1]
bound = sqrt(6 / (numIn + numOut));

weights = bound.*Z;

end

end
end
