%% Description
% Create 1 x 5 regime plot for ClassNNEns, sliced across different axis

function RegimeVisualizationFun(ClassNNEns,RegimeList,ReactantConcUBLB,ReactantConcSpaceAxisName,Filepath)

RegimeVisualizePointsPerDim = 100;
nCollage = 5;
TransitionThreshold = 0.8;

RegimePlotFilepath = Filepath + "/Regime Plots";
if exist (RegimePlotFilepath,"file") == 0
    mkdir(RegimePlotFilepath)
end

%% Generate regime plot without cut off

RegimePlotFun(ClassNNEns,RegimeList,ReactantConcUBLB,RegimeVisualizePointsPerDim,ReactantConcSpaceAxisName,0,RegimePlotFilepath)

%% Generate regime plot with cut off

RegimePlotFun(ClassNNEns,RegimeList,ReactantConcUBLB,RegimeVisualizePointsPerDim,ReactantConcSpaceAxisName,TransitionThreshold,RegimePlotFilepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Important functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Regime plot function

% Description
% Creates regime plots, which slices across the [MA]_0, [BrO3-]_0 and [Ce3+]_0 axis

function RegimePlotFun(ClassNNEns,RegimeList,ReactantConcSpaceBound,RegimeVisualizePointsPerDim,ReactantConcSpaceAxisName,TransitionThreshold,RegimePlotFilepath)

%% Create points to predict

% Create feed concentration plotting points
Af1DFull = linspace(ReactantConcSpaceBound(1,1),ReactantConcSpaceBound(1,2),RegimeVisualizePointsPerDim);
Bf1DFull = linspace(ReactantConcSpaceBound(2,1),ReactantConcSpaceBound(2,2),RegimeVisualizePointsPerDim);
Cf1DFull = linspace(ReactantConcSpaceBound(3,1),ReactantConcSpaceBound(3,2),RegimeVisualizePointsPerDim);


for iCollage = nCollage

    % Create feed concentration spacing in the z direction
    AfZdir = linspace(ReactantConcSpaceBound(1,1),ReactantConcSpaceBound(1,2),iCollage);
    BfZdir = linspace(ReactantConcSpaceBound(2,1),ReactantConcSpaceBound(2,2),iCollage);
    CfZdir = linspace(ReactantConcSpaceBound(3,1),ReactantConcSpaceBound(3,2),iCollage);

    for zAxisIdx = 1:3 % 1,2,3 being Af,Bf and Cf as the z direction respectively
    
        Af1D = Af1DFull;
        Bf1D = Bf1DFull;
        Cf1D = Cf1DFull;
    
        if zAxisIdx == 1
            Af1D = AfZdir;
            zList = AfZdir;
            xAxisIdx = 2;
            yAxisIdx = 3;
        elseif zAxisIdx == 2
            Bf1D = BfZdir;
            zList = BfZdir;
            xAxisIdx = 1;
            yAxisIdx = 3;
        elseif zAxisIdx == 3
            Cf1D = CfZdir;
            zList = CfZdir;
            xAxisIdx = 1;
            yAxisIdx = 2;
        end
    
        % Create plotting points
        [AF,BF,CF] = ndgrid(Af1D,Bf1D,Cf1D);
        Af = reshape(AF,[],1);
        Bf = reshape(BF,[],1);
        Cf = reshape(CF,[],1);
        clear AF BF CF

    
        MembershipPred = NNEnsPredict(ClassNNEns,[Af,Bf,Cf]); % Identify membership from NN ens
        MembershipPredAvg = mean(MembershipPred,3,"omitnan");
        MaxMembershipPred = max(MembershipPredAvg,[],2);
        RegimeIdx = double(string(onehotdecode(MembershipPredAvg ,RegimeList,2)));
    
        h = figure;
        h.Position = [0 0 180*iCollage 180];
        counter = 1;
     
        for z = zList
    
            if zAxisIdx == 1
                RegimeIdxPlot = RegimeIdx(Af==z,:);
                MaxMembershipPlot = MaxMembershipPred(Af==z,:);
            elseif zAxisIdx == 2
                RegimeIdxPlot = RegimeIdx(Bf==z,:);
                MaxMembershipPlot = MaxMembershipPred(Bf==z,:);
            elseif zAxisIdx == 3
                RegimeIdxPlot = RegimeIdx(Cf==z,:);
                MaxMembershipPlot = MaxMembershipPred(Cf==z,:);
            end
            
            [R,G,B] = RegimeColor(RegimeIdxPlot,MaxMembershipPlot,TransitionThreshold); % Assign color based on cluster index
            
            % Reshape and resize RGB vectors to generate image matrix of suitable size.
            R =flipud(rot90(reshape(R,RegimeVisualizePointsPerDim,[]),1));           
            R = imresize(R,[size(R,1)*2,size(R,1)*3]);
            G =flipud(rot90(reshape(G,RegimeVisualizePointsPerDim,[]),1));           
            G = imresize(G,[size(G,1)*2,size(G,1)*3]);
            B =flipud(rot90(reshape(B,RegimeVisualizePointsPerDim,[]),1));           
            B = imresize(B,[size(B,1)*2,size(B,1)*3]);
            RegimeRGBPlot = cat(3,R,G,B);
    
            % Create regime collage
            subplot(1,iCollage,counter)
            hold on
            grid minor
            imshow(RegimeRGBPlot,"Interpolation","bilinear")
            axis on
            ax = gca;
            ax.YDir = 'normal';
            set(ax,'XTick',linspace(0,size(RegimeRGBPlot,2),5))
            set(ax,'XTickLabel',round(linspace(ReactantConcSpaceBound(xAxisIdx,1),ReactantConcSpaceBound(xAxisIdx,2),5),2,"significant"))
            set(ax,'YTick',linspace(0,size(RegimeRGBPlot,1),5))
            set(ax,'YTickLabel',round(linspace(ReactantConcSpaceBound(yAxisIdx,1),ReactantConcSpaceBound(yAxisIdx,2),5),2,"significant"))
            xlim([0.5,size(RegimeRGBPlot,2)])
            ylim([0.5,size(RegimeRGBPlot,1)])
            xlabel(ReactantConcSpaceAxisName{xAxisIdx} + "(mol/dm^{3})","FontSize",7);
            ylabel(ReactantConcSpaceAxisName{yAxisIdx} + "(mol/dm^{3})","FontSize",7);    
            title([ReactantConcSpaceAxisName{zAxisIdx} , " = " + num2str(z) + "mol/dm^{3}"],"FontSize",8) 
            hold off
    
            counter = counter + 1;
        end

        ImageName = "z axis = "+ReactantConcSpaceAxisName{zAxisIdx}+", 1x" + num2str(iCollage) +"CutOff = "...
                     + num2str(TransitionThreshold)+ " Regime plot.jpeg";
        saveas(h,RegimePlotFilepath+"/"+ImageName)
    end
end


%% Cluster color function

% Description
% Set color based on regime.

function [R,G,B] = RegimeColor(ClusterIdx,MaxMembershipPlot,TransitionThreshold)

ClusterColor1  = [  0 225   0];  % Lime
ClusterColor2  = [  0   0 255];  % Blue
ClusterColor3  = [255   0 255];  % Fuchsia
ClusterColor4  = [255 255   0];  % Yellow
ClusterColor5  = [  0 191 255];  % Deep sky blue
ClusterColor6  = [255   0   0];  % Red


ClusterColorMatrix = [ClusterColor1;ClusterColor2;ClusterColor3;ClusterColor4;ClusterColor5;ClusterColor6];

R = ClusterColorMatrix(ClusterIdx,1)./255;
G = ClusterColorMatrix(ClusterIdx,2)./255;
B = ClusterColorMatrix(ClusterIdx,3)./255;

% Set RGB for Max membership values < cut off limit to 0
R(MaxMembershipPlot<TransitionThreshold) = 0;
G(MaxMembershipPlot<TransitionThreshold) = 0;
B(MaxMembershipPlot<TransitionThreshold) = 0;

end


end

end