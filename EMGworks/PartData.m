function [indTr,indTs] = PartData(instant,Ntest)

N = length(instant);


[Trainning, Test] = crossvalind('LeaveMOut', N, Ntest);

indTr = find(Trainning==1);

indTs = find(Test==1);

end


