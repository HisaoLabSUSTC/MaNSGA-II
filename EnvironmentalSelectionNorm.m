function [Population,FrontNo,CrowdDis] = EnvironmentalSelectionNorm(Population,N,IdealPoint,alpha)
% The environmental selection of MaNSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    PopObj = Population.objs; 
    [~,M] = size(PopObj);
    NDSol = PopObj(NDSort(PopObj,Population.cons,1)==1,:);
    NadirPoint = max(NDSol,[],1);
    range = NadirPoint - IdealPoint; 
    if range == 0
        range = 10e-6;
    end
    PopObj = (PopObj - IdealPoint)./range;
    PopObj = (1-alpha)*PopObj+(alpha/M)*sum(PopObj,2);
    IdealPoint = zeros(1,M);
    [FrontNo,MaxFNo] = NDSort(PopObj,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    if MaxFNo > 1
        %% Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(PopObj,FrontNo);
    
        %% Select the solutions in the last front based on their crowding distances
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next(Last(Rank(1:N-sum(Next)))) = true;
    
        %% Population for next generation
        Population = Population(Next);
        FrontNo    = FrontNo(Next);
        CrowdDis   = CrowdDis(Next);
    else
        %% Get the solutions in the first front
        Population = Population(FrontNo==1);
        PopObj = PopObj(FrontNo==1,:);
        FrontNo = ones(1,size(PopObj,1));
           
        %% Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(PopObj,FrontNo);
        
        %% Select the solutions in the first front using distance-based subset selection 
        [~,selind] = DisSel(PopObj,IdealPoint,N);
        
        %% Population for next generation 
        Population = Population(selind);
        FrontNo    = FrontNo(selind);
        CrowdDis   = CrowdDis(selind);
    end
end