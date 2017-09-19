function [] = SavingCombinatorics(Name,numTrials, ResultsStruct)

for i = 1:(length(ResultsStruct)/numTrials)
    TPos = 0; FPos = 0; TNeg = 0; FNeg = 0;
    
    for j = (i-1)*numTrials+1: i*numTrials
        TPos = TPos + ResultsStruct(j).TP;
        FPos = FPos + ResultsStruct(j).FP;
        TNeg = TNeg + ResultsStruct(j).TN;
        FNeg = FNeg + ResultsStruct(j).FN;
    end
    
    ResultsCombinatorics(i) = struct('Trials',{ResultsStruct(i*numTrials).Trial},'k',...
        {ResultsStruct(i*numTrials).k},'Features',{ResultsStruct(i*numTrials).Features},'TP',...
        {TPos},'FP',{FPos},'TN',{TNeg},'FN',{FNeg},'beta',{ResultsStruct(i).beta});
end

save(['ResultCombinatorics_',Name,'.mat'],'ResultsCombinatorics')
end