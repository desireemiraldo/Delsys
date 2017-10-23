Var = {'ACCF','ACCPitchF','ACCRollF','ACCYawF',...
    'GYRF','GYRPitchF','GYRRollF','GYRYawF',...
    'MAGF','MAGPitchF','MAGRollF','MAGYawF','Cte'};%,'AccF.^2'}; %Used in linear combination

Combinations = repmat(struct('Features',{0}),8191,1);

n = 0;
for k = 1: length(Var)
    % -- Combinatorial analysis to find the best set of variables
    % to be used in linear combination
    CombinatoricsInd = nchoosek(1:1:length(Var),k);
    
    for kk = 1 : size (CombinatoricsInd,1)
        n = n+1;
        Features = CombinatoricsInd(kk,:);
        
        Combinations(n).Features = Var(Features);
    end
end

save('FeaturesComb1.mat','Combinations')