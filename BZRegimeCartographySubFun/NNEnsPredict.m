%% Description
% This function generates the prediciton of the neural network ensemble for input x.

function Score = NNEnsPredict(ClassNNEns,x)

% dummy round to identify score matrix size
NumEns = size(ClassNNEns,1); % number of ensemble members
MdlDummy = ClassNNEns{1,1};
ScoreDummy = MdlDummy(x(1,:)');
Score = zeros(size(x,1),size(ScoreDummy,1),NumEns);


for i  = 1:NumEns
    Mdl = ClassNNEns{i,1};
    Score1 = Mdl(x');
    Score(:,:,i) = Score1';
end

end