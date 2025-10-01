%% Description %%
% This function calculates the performance of the current neural network by comparing with test
% dataset. The performance metrics (accuracy, precision, recall and macro F1-scores) are calculated using the
% formula proposed by Binaghi et.al. (https://doi.org/10.1016/S0167-8655(99)00061-6) and Harju et.al.
% (https://doi.org/10.48550/arXiv.2309.13938) so it better reflects the performance of the soft classification 
% nature of the neural network ensembles.

%%
function NNEnsPerformanceRec = NNEnsPerformanceDiag(ClassNNEns,InputTest,RegimeMembershipTest,NNEnsPerformanceRec,ActiveLearningFolderPath)

%% Calculate prediction from curent iter ClassNNENs

RegimeMembershipPred = NNEnsPredict(ClassNNEns,InputTest);
RegimeMembershipPredAvg = mean(RegimeMembershipPred,3,"omitnan"); % reference probability

%% Calculate difference in prediction by ClassNNEns and ground truth and record accuracy, precision, recall, and F1-score

% Accuracy calculation follows formula of overall accuracy by Binaghi et.al. (https://doi.org/10.1016/S0167-8655(99)00061-6)
% Precision, recall and F1 score calculated based on method proposed by Harju et.al. (https://doi.org/10.48550/arXiv.2309.13938)

FuzzySetIntersect = min(RegimeMembershipPredAvg,RegimeMembershipTest);

% Calculate overall accuracy
Accuracy = (sum(FuzzySetIntersect,"all")/size(RegimeMembershipPredAvg,1))*100;

% Calculate Macro-averaged precision
Precision = sum(FuzzySetIntersect,1)./sum(RegimeMembershipPredAvg,1);
MacroAvgPrecision = mean(Precision,"all");
PrecisionRec = [Precision,MacroAvgPrecision];

% Calculate Macro-averaged recall
Recall = sum(FuzzySetIntersect,1)./sum(RegimeMembershipTest,1);
MacroAvgRecall = mean(Recall,"all");
RecallRec = [Recall,MacroAvgRecall];

% Calculate Macro-averaged F1 score
F1score = 2.*sum(FuzzySetIntersect,1)./sum(RegimeMembershipTest+RegimeMembershipPredAvg,1);
MacroAvgF1score = mean(F1score,"all");
F1scoreRec = [F1score,MacroAvgF1score];

% Note: within the recall, precision and F1-score matrix, the first 8 columns are the individual class metric
% and the final column is the macro averaged of all class metrics.

if isempty(NNEnsPerformanceRec)
    NNEnsPerformanceRec.Accuracy  = Accuracy;
    NNEnsPerformanceRec.Precision = PrecisionRec;
    NNEnsPerformanceRec.Recall    = RecallRec;
    NNEnsPerformanceRec.F1score   = F1scoreRec;
else
    NNEnsPerformanceRec.Accuracy  = [NNEnsPerformanceRec.Accuracy;Accuracy];
    NNEnsPerformanceRec.Precision = [NNEnsPerformanceRec.Precision;PrecisionRec];
    NNEnsPerformanceRec.Recall    = [NNEnsPerformanceRec.Recall;RecallRec];
    NNEnsPerformanceRec.F1score   = [NNEnsPerformanceRec.F1score;F1scoreRec];
end

%% Plot the evolution of NNEnsAccuracyRec

AccuracyFig = figure;
hold on
plot([0:(size(NNEnsPerformanceRec.Accuracy,1)-1)]',NNEnsPerformanceRec.Accuracy,".-","LineWidth",1.5,"MarkerSize",15)
grid minor
title("NN ensemble accuracy vs Iteration")
xlabel("Iteration")
ylabel("Accuracy")
ytickformat('percentage')
ylim([-inf 100])
hold off

PrecisionFig = figure;
hold on
plot([0:(size(NNEnsPerformanceRec.Precision,1)-1)]',NNEnsPerformanceRec.Precision(:,1:end-1),"--.","LineWidth",1,"MarkerSize",12)
plot([0:(size(NNEnsPerformanceRec.Precision,1)-1)]',NNEnsPerformanceRec.Precision(:,end),".-","LineWidth",1.5,"MarkerSize",15)
grid minor
title("NN ensemble precision vs Iteration")
xlabel("Iteration")
ylabel("Precision")
ylim([-inf 1])
legend(["Regime 1","Regime 2","Regime 3","Regime 4","Regime 5","Regime 6","Macro-averaged"],'Location','eastoutside')
hold off

RecallFig = figure;
hold on
plot([0:(size(NNEnsPerformanceRec.Recall,1)-1)]',NNEnsPerformanceRec.Recall(:,1:end-1),"--.","LineWidth",1,"MarkerSize",12)
plot([0:(size(NNEnsPerformanceRec.Recall,1)-1)]',NNEnsPerformanceRec.Recall(:,end),".-","LineWidth",1.5,"MarkerSize",15)
grid minor
title("NN ensemble recall vs Iteration")
xlabel("Iteration")
ylabel("Recall")
ylim([-inf 1])
legend(["Regime 1","Regime 2","Regime 3","Regime 4","Regime 5","Regime 6","Macro-averaged"],'Location','eastoutside')
hold off

F1scoreFig = figure;
hold on
plot([0:(size(NNEnsPerformanceRec.F1score,1)-1)]',NNEnsPerformanceRec.F1score(:,1:end-1),"--.","LineWidth",1,"MarkerSize",12)
plot([0:(size(NNEnsPerformanceRec.F1score,1)-1)]',NNEnsPerformanceRec.F1score(:,end),".-","LineWidth",1.5,"MarkerSize",15)
grid minor
title("NN ensemble F1-score vs Iteration")
xlabel("Iteration")
ylabel("F1-score")
ylim([-inf 1])
legend(["Regime 1","Regime 2","Regime 3","Regime 4","Regime 5","Regime 6","Macro-averaged"],'Location','eastoutside')
hold off

%% Create folder and save results


saveas(AccuracyFig,ActiveLearningFolderPath+"Accuracy vs iteration.jpg")
saveas(PrecisionFig,ActiveLearningFolderPath+"Precision vs iteration.jpg")
saveas(RecallFig,ActiveLearningFolderPath+"Recall vs iteration.jpg")
saveas(F1scoreFig,ActiveLearningFolderPath+"F1-score vs iteration.jpg")
close("all")
save(ActiveLearningFolderPath+"NNEnsPerformanceRec.mat","NNEnsPerformanceRec")


end