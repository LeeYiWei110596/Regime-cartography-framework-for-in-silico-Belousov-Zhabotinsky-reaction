%% Description
% This function creates test data for the Batch BZ models

%%
function [ReactantConcTest,RegimeMembershipTest] = BZTestDataCreateFun(BatchBZSimulation, ...
                                                                       BrionConcFeatureVecConverterFun, ...
                                                                       RegimeCentroid, RegimeList, m, ...
                                                                       ReactantConcUBLB, Eigenvec, ...
                                                                       muBrionConc, ...
                                                                       OscillateNoRescaledCoeff, ...
                                                                       OscillatePeriodRescaledCoeff, ...
                                                                       ReactantConcSpaceAxisName, ...
                                                                       TestDataFolderpath)
% Create test data through latin hypercube sampling
nReactantConcTest = 1000; 
CondTest = (601:(600+nReactantConcTest))'; % condition 600 and above are test conditions
nDim = size(ReactantConcUBLB,1);
ReactantConcNormTest = lhsdesign(nReactantConcTest,nDim,"iterations",100);

% Input feed concentration through Batch BZ simulation handle
ReactantConcTest = ReactantConcNormTest.*(ReactantConcUBLB(:,2) - ReactantConcUBLB(:,1))' + ReactantConcUBLB(:,1)';
BrionConcTest = BatchBZSimulation(ReactantConcTest);

% Create feature vector for Brion concentration 
[BrionConcTestFeatureVec,BrionConcTestFeature] = BrionConcFeatureVecConverterFun(BrionConcTest,Eigenvec,muBrionConc,OscillateNoRescaledCoeff,OscillatePeriodRescaledCoeff);

% Idetify regime memberhsip for each point
[RegimeMembershipTest,RegimeIdxTest] = RegimeAllocFun(BrionConcTestFeatureVec,RegimeCentroid,RegimeList,m);

% Save test data reactant concentration into excel sheet
ReactantConcTestTable = array2table([CondTest,ReactantConcTest]);
ReactantConcTestTable.Properties.VariableNames = ["Cond",string(ReactantConcSpaceAxisName)];
writetable(ReactantConcTestTable,TestDataFolderpath + "/ReactantConcTest.xlsx")

% Save test data Brion conc, Brion conc feature vector and Brion conc feature into excel sheet
writematrix(BrionConcTest,TestDataFolderpath + "/BrionConcTest.csv")
writematrix(BrionConcTestFeatureVec,TestDataFolderpath + "/BrionConcTestFeatureVec.csv")
BrionConcTestFeatureTable = array2table([CondTest,BrionConcTestFeature]);
BrionConcTestFeatureTable.Properties.VariableNames = ["Cond","f_osc","n_osc"];
writetable(BrionConcTestFeatureTable,TestDataFolderpath + "/BrionConcTestFeature.csv")

% Save all data into a test data rec file
RegimeIdxTestTable = array2table(RegimeIdxTest);
RegimeIdxTestTable.Properties.VariableNames = "Regime Idx";
RegimeMembershipTestTable = array2table(RegimeMembershipTest);
RegimeMembershipTestTable.Properties.VariableNames = ["Regime 1 Membership","Regime 2 Membership",...
                                                      "Regime 3 Membership","Regime 4 Membership",...
                                                      "Regime 5 Membership","Regime 6 Membership"];
TestDataRecTable = [ReactantConcTestTable,RegimeIdxTestTable,RegimeMembershipTestTable,BrionConcTestFeatureTable(:,2:end)];
writetable(TestDataRecTable,TestDataFolderpath + "/TestDataRec.xlsx")

end