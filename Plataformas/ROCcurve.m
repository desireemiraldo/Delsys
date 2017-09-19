function [Sensib,Specif] = ROCcurve(ResultsStruct)

Sensib = NaN (length(ResultsStruct),1);
Specif = NaN (length(ResultsStruct),1);

for i = 1: length(ResultsStruct)
    
    Sensib(i) = ResultsStruct(i).TP/...
        (ResultsStruct(i).TP + ResultsStruct(i).FN);
    
    Specif(i) = ResultsStruct(i).TN/...
        (ResultsStruct(i).FP + ResultsStruct(i).TN);
end

figure;
plot(1-Specif,Sensib,'.','MarkerSize',12)
xlim([0 1]); ylim = ([0 1]);
line([0 1], [0 1],'LineStyle','--','color','k')
xlabel('1-Specif');
ylabel('Sensib');

end