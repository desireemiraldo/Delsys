function [ResultsCombinatorics] = combinatorics(Var,Sensors,p,y,pct,ToeOff,Fy,delay,time,indTr,indTs)

LinearCombination = NaN(size(time));

TotalComb = 0;
for k = 1: length(Var)
    b = nchoosek(length(Var),k);
    TotalComb = TotalComb + b;
end

ResultsCombinatorics = struct('TrialsTr',{0},'Features',{0},'Threshold',{0},'TP',{0},...
    'FP',{0},'TN',{0},'FN',{0},'beta',{0});
ResultsCombinatorics = repmat(ResultsCombinatorics,TotalComb,1);

%% --

n = 0;
for k = 1: length(Var)
    % -- Combinatorial analysis to find the best set of variables
    % to be used in linear combination
    CombinatoricsInd = nchoosek(1:1:length(Var),k);
    
    for kk = 1 : size (CombinatoricsInd,1)
        TPos = 0; TNeg = 0; FPos = 0; FNeg = 0;
        
        Features = CombinatoricsInd(kk,:);
        
            
            pTr = p(:,Features,indTr);
            yTr = y(:,indTr);
            
            pTs = p(:,Features,indTs);
            
            beta = ols(pTr,yTr); % Coefs for linear combination
            betaM = mean(beta,2);
            
            % --- Applying Linear Combination
            
            for j = 1:length(indTs)
                for i = 1 : length(Sensors)
                    
                    index = indTs(j)+(i-1)*length(indTs);
                    
                    TP = 0; FP = 0; TN = 0; FN = 0;
                    n = n+1; %disp(n)
                    
                    LinearCombination(:,index) = pTs(:,:,(j)+(i-1)*length(indTs))*betaM;
                    
                    %% Plot Fy
                    
%                     figure(index)
%                     subplot(2,1,1); plot(Fy(:,1,index)+delay, Fy(:,2,index),'linewidth',1.5,'Color',[0.3 0.3 0.3]);
%                     hold on;
%                     plot(time(:,index),20*y(:,index),'r','linewidth',1.5)
%                     xlim([time(1,index),time(end,index)]);
%                     xticklabels({})
%                     ylabel('For�a [kg/N]','fontsize',12)
%                     %title([' Trial ', num2str(index)]);
%                     
%                     subplot(2,1,2); plot(time(:,index),LinearCombination(:,index),'linewidth',1.5);
%                     ylim([min(LinearCombination(:,index))*1.1 max(LinearCombination(:,index))*1.1]);
%                     xlim([time(1,index),time(end,index)]);
%                     ylabel('Combina��o Linear','fontsize',12)
%                     xlabel('Tempo [s]','fontsize',12)
%                     
                    %  --- Checking the combination's quality
                    threshold = (max(max(LinearCombination(:,:))))*pct;
                    [pks,locs] = findpeaks(LinearCombination(:,index),time(:,index),'MinPeakHeight',threshold);
                    
%                     for z = 1:length(locs)
%                         Line = line([locs(z) locs(z)], [-0.15 3],'Linewidth',1.5,'Linestyle','--','Color',[0 0 0]);
%                         set(Line,'Clipping','off')
%                     end
%                     
                    ind = [Fy(1,1,index), ToeOff(index,:),Fy(end,1,index)+0.1];
                    ind(isnan(ind))=[];
                    
                    for w = 2: length(ind)
                        center = ((ind(w)-0.1) + ind(w-1))/2;
                        radius = abs(ind(w)-0.1 - center);
                        
                        if w~= length(ind)
                            A = rangesearch(locs,ind(w)-0.05,0.05);
                            if isempty(A{1})
                                FN = FN+1;
                            else
                                TP = TP+1;
                            end
                        end
                        
                        if center >= ind(1)
                            A = rangesearch(locs,center,radius);
                            if isempty(A{1})
                                TN = TN+1;
                            else
                                FP = FP+1;
                            end
                        end
                    end
                    
                    %---
                    TPos = TPos + TP;
                    TNeg = TNeg + TN;
                    FPos = FPos + FP;
                    FNeg = FNeg + FN;
                    
                    % --- save
                    %                 ResultsindTrs(n) = struct('Trial',{indTs(j)},...
                    %                     'Features',{Var(Features)},'Locs',{locs},'Threshold',...
                    %                     {pct},'TP',{TP},'FP',{FP},'TN',{TN},'FN',{FN},'beta',{beta});
                    %                 % structSave = struct('teste',{File,k,Features,TP,FN,TN,FP,beta});
                    
                    % keyboard %breakpoint
                end
            end
      
        %disp(n/length(indTs))
        ResultsCombinatorics(n/length(indTs)) = struct('TrialsTr',{indTr},...
            'Features',{Var(Features)},'Threshold',{pct},'TP',{TPos},...
            'FP',{FPos},'TN',{TNeg},'FN',{FNeg},'beta',{beta});
    end
end


end

