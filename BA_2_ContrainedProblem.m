clc;
clear;
close all;

%% Problem Definition
ObjFunction = @(x) objective(x);    % Objective Function

Dims = 2;                            % Number of Decision Variables
VarMin = [-100 -100];                        % Decision Variables Lower Bound
VarMax = [100 100];                         % Decision Variables Upper Bound

%% Bees Algorithm Parameters
MaxIt = 2000;                        % Maximum Number of Iterations
n = 20;                         % Number of Scout Bees
nep = 20;                         % Number of Recruited Bees for Elite Sites
Shrink = 1;

%% Initialization
Bees.Position = [];
Bees.Cost = [];
Bees = repmat(Bees, n, 1);

% Generate Initial Solutions
for i = 1:n
    Bees(i).Position = unifrnd(VarMin, VarMax, [1, Dims]);
    Bees(i).Cost = ObjFunction(Bees(i).Position);
end
recruitment = ceil(linspace(nep,1,n));
size= linspace(0,1,n);

% Array to Hold Best Cost Values
OptCost = zeros(MaxIt, 1);

%% Bees Algorithm Local and Global Search
for it = 1:MaxIt
    % All Sites (Exploitation and Exploration)
    for i = 1:n
        bestNewBee.Cost = inf;
        assigntment = linspace(0,1,recruitment(i));
        %assigntment = D_Tri_real_array(0,size(i),1,1,recruitment(i));
        for j = 1:recruitment(i)
            newBee.Position = foraging(Bees(i).Position, assigntment(j), VarMin, VarMax, Shrink);            newBee.Cost = ObjFunction(newBee.Position);
            if newBee.Cost < bestNewBee.Cost
                bestNewBee = newBee;
            end
        end
        
        if bestNewBee.Cost < Bees(i).Cost
            Bees(i) = bestNewBee;
        end
    end
    
    % Update Best Solution Ever Found
    [~, RankOrder] = sort([Bees.Cost]);
    Bees = Bees(RankOrder);
    OptSol = Bees(1);
    OptCost(it) = OptSol.Cost;
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(OptCost(it)) ': Best Pos = ' num2str(OptSol.Position)]);
    Shrink = Shrink * 0.99;
end

function z = objective(x)
    z = (x(1) - 3)^2 + (x(2) - 3)^2 + penalty(x);
end

function p = penalty(x)
    constraint1 = @(x) x(1) + x(2) - 4;
    constraint2 = @(x) x(1) - 3 * x(2) - 1;
    c1 = constraint1(x);
    c2 = constraint2(x);
    
    penaltyFactor = 1e4;
    p = penaltyFactor * max(0, c1)^2 + penaltyFactor * min(0, c2)^2;
end

function y = foraging(x,assign, VMn, VMx, Shrink)
    nVar = numel(x);
    ngh = triangular(0,assign,1);
    r = ngh * (VMx-VMn) * Shrink;

    y = x + unifrnd(-r, r, size(x));
    
    for i=1:nVar
        y(i) = max(y(i), VMn(i));
    end
end

function random_number = triangular(min, mode, max)
    r = rand();
    prop_mode = (mode - min) / (max - min);
    if r <= prop_mode
        random_number = min + sqrt(r * (max - min) * (mode - min));
    else
        random_number = max - sqrt((1 - r) * (max - min) * (max - mode));
    end
end