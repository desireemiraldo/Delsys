function Variable = SelectVar(Data, ChannelType, Signal)

[Sample,Channel,Trials] = size(Data);

for k = 1:Trials
    index = 0;
    for i = 1: Sample
        if  strcmp(Data(i,1),ChannelType) && strcmp(Data(i,2),Signal)
            index = index +1;
            Variable(index,:,k) = cell2mat(Data(i,4:end,k));
            
        end
    end
end

% for k = 1:Trials
%     Variable(:,1,k) = Variable(:,1,k) - Variable(1,1,k);
% end

end