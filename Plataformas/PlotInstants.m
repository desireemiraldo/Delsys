function PlotInstants( instant,File )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


for i=1:length(instant.textdata)
   if strcmp(instant.textdata(i,1),File)
       if strcmp(instant.textdata(i,2),'r')
               for j = 1:2: size(instant.data,2)
                   if instant.data(i,j)>0
                       HeelContact = line([instant.data(i,j) instant.data(i,j)], [-2 1000],...
                           'Color',[0 0 1]);
                       ToeOff = line([instant.data(i,j+1) instant.data(i,j+1)], [-2 1000],...
                           'Linewidth',1,'Linestyle','--','Color',[0 0 1]);
                       set(HeelContact,'Clipping','off')
                       set(ToeOff,'Clipping','off')
                   end
               end
       end
        if strcmp(instant.textdata(i,2),'l')
               for j = 1:2: size(instant.data,2)
                   if instant.data(i,j)>0
                       HeelContact = line([instant.data(i,j) instant.data(i,j)], [-5 1000],...
                           'Color',[1 0 0]);
                       ToeOff =line([instant.data(i,j+1) instant.data(i,j+1)], [-5 1000],...
                           'Linewidth',1,'Linestyle','--','Color',[1 0 0]);
                       set(HeelContact,'Clipping','off')
                       set(ToeOff,'Clipping','off')
                   end
               end
       end
   end
end

end

