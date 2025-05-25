classdef MaNSGAII_Norm < ALGORITHM
% <2025> <multi> <real/integer/label/binary/permutation> <constrained/none>
% Many-objective Nondominated sorting genetic algorithm II
% a --- 0.01 --- tradeoff 

%------------------------------- Reference --------------------------------
% L.M.Pang, H.Ishibuchi, K.Deb, K.Shang, MaNSGA-II: Many-objective NSGA-II, 
% IEEE Transactions on Emerging Topics in Computational Intelligence, 2025.
% (This is a version with normalization mechanism)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %%Parameter Setting 
            a = Algorithm.ParameterSet(0.01);
            alpha = (a * Problem.M)/(1 - a + a * Problem.M);
            
            %% Generate random population
            Population = Problem.Initialization();
            IdealPoint = min(Population.objs);
            [~,FrontNo,CrowdDis] = EnvironmentalSelectionNorm(Population,Problem.N,IdealPoint,alpha);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                IdealPoint = min([IdealPoint;Offspring.objs]);
                [Population,FrontNo,CrowdDis] = EnvironmentalSelectionNorm([Population,Offspring],Problem.N,IdealPoint,alpha);
            end
        end
    end
end