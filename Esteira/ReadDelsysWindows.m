function [Trials,deltaT, Label] = ReadDelsysWindows(FilePath, winPts)

A =  importdata(FilePath);
Label = [A.textdata(1,:),num2cell(A.data(1,:))];
A =  importdata(FilePath,',',1);
A.textdata(1,:) = [];
Data = [A.textdata, num2cell(A.data)];

linesNumber = 78*winPts;
BlocksNumber = int8((length(Data)-1)/linesNumber);

%Trials = zeros(linesNumber,size(Data,2),BlocksNumber);

index = 1; 

while index <= BlocksNumber 
    Trials(:,:,index) = Data(linesNumber*(index-1)+1 : linesNumber*index,:);
    Trials(:,4,index) = num2cell(str2double(Trials(:,4,index)));
    index = index +1;
end

delta(:,1,1) = Trials(1,4,:);
delta(:,2,1) = Trials(end,4,:);

deltaT = cell2mat(delta);

end