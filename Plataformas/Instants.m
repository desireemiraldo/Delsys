function [HeelContact, ToeOff] = Instants(instant,File)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
HeelContact = [];
k=1;

if strcmp(instant.textdata(1,1),'﻿File')
    instant.textdata = instant.textdata(2:end,:);
end
for i=1:2:size(instant.textdata,1)
    if strcmp(instant.textdata(i,1),File)
        for j = 1 :2: size(instant.data,2)
            HeelContact(1,k) = instant.data(i,j);
            HeelContact(2,k) = instant.data(i+1,j);
            
            ToeOff(1,k) = instant.data(i,j+1);
            ToeOff(2,k) = instant.data(i+1,j+1);
            
            k=k+1;
        end
    end
end
for i = 1: size(HeelContact,1)
    for j = 1: size(HeelContact,2)
        if HeelContact(i,j) == -1
            HeelContact(i,j) = 0;
        end
        if ToeOff(i,j) == -1
            ToeOff(i,j) = 0;
        end
        
    end
end
end
               