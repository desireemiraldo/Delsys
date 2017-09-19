function [ToeOff, HeelStrike] = ReshapeInstants1(deltaT, instant,Name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

A = strcmp(instant.textdata(:,1),Name);
ind = find(A==1)-1; %first row is label

ToeOff = NaN;
HeelStrike = NaN;

for k = 1 : length(deltaT)
    
    first = deltaT(k,1);
    last = deltaT(k,2);
    
    i = find(instant.data(ind(1),:)>=first & instant.data(ind(1),:)<=last); %edit when using both legs
    
    NewHS = instant.data(ind(1),i(find(rem(i,2)==1))); %- deltaT(k,1);
    NewTO = instant.data(ind(1),i(find(rem(i,2)==0))); %- deltaT(k,1);
       
    ColNumber = size(ToeOff,2);
    
    if length(NewTO)> ColNumber
        ToeOff(:,ColNumber + 1 :length(NewTO)) = NaN;
    end
    if length(NewTO)< ColNumber
        ToeOff(k,:) = NaN;
    end
    ToeOff(k,1:length(NewTO))  = NewTO;
    
    
    ColNumber = size(HeelStrike,2);
    
    if length(NewHS)> ColNumber
        HeelStrike(:,ColNumber + 1 :length(NewHS)) = NaN;
    end
    if length(NewHS)< ColNumber
        HeelStrike(k,:) = NaN;
    end
    HeelStrike(k,1:length(NewHS))  = NewHS;

end



end

