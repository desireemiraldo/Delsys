%% Gait Analysis (w/ accelerometer)
clear all; clc; close all

%% -- About the files

Path = '.\Piloto\';

Folder = {'RNW\'};

Files = {'RNW_Calçado_Confortavel_Rep_1.1',...
    'RNW_Calçado_Confortavel_Rep_2.2',...
    'RNW_Calçado_Confortavel_Rep_3.3',...
    'RNW_Calçado_Confortavel_Rep_4.4',...
    'RNW_Calçado_Confortavel_Rep_5.5',...
    'RNW_Descalco_Confortavel_Rep_1.17',...
    'RNW_Descalco_Confortavel_Rep_2.18',...
    'RNW_Descalco_Confortavel_Rep_4.20',...
    'RNW_Descalco_Confortavel_Rep_5.21',...
    'RNW_Descalco_Confortavel_Rep_6.22'};

Sensor = {[1,2]};

Win = {'Gauss'}; %,'Rect'};

%% -- Initializing some variables

% --
Signal = {'ACC X', 'ACC Y', 'ACC Z',...
    'Gyro X', 'Gyro Y', 'Gyro Z',...
    };%'Mag X', 'Mag Y', 'Mag Z',};

Var = {'ACCF','ACCXF','ACCYF','ACCZF',...
    'GyroF','GyroXF','GyroYF','GyroZF',...
    };%'MagF','MagXF','MagYF','MagZF','Cte'};%,'AccF.^2'}; %Used in linear combination

TotalComb = 0;
for k = 1: length(Var)
    b = nchoosek(length(Var),k);
    TotalComb = TotalComb + b;
end

% -- Filter
Fs = 148.14;
[t,Wn] = buttord(10/(Fs/2),20/(Fs/2),1,60);
[b,a] = butter(t, Wn);

% -- Standard deviation to be used in Linear Combination
sd = 50e-3;

% --
deltaT = 2.75; %seconds/window

% --

RC = struct('TrialsTr',{0},'Features',{0},'Threshold',{0},'TP',{0},...
    'FP',{0},'TN',{0},'FN',{0},'beta',{0});
RC = repmat(RC,TotalComb*101,1);

% % p = NaN(2*winPts,length(Var),2*22);
% % y = NaN(2*winPts,2*22);
% % ForceY = NaN(2*60*FsFP,8);
% % timeACC = NaN(size(y));
ToeOff = [];
HeelStrike = [];

for Sub = 1: length(Folder)
    for w = 1: length(Win)
        
        instant = importdata([Path,Folder{Sub},'Instantes.txt'],';');
        
        %% APAGAR DEPOIS
        for x = 2:11
            instant.textdata(x,1) = strrep(instant.textdata(x,1),'Calcado','Calçado');
        end
        %% 
        
        % Data = BuildingTrials([Path,Folder{Sub}],Files, 3, deltaT, instant);
        % save('RNW_data.mat','Data');
        load RNW_data.mat
        
        delta(:,1,1) = Data.data(1,1,:);
        delta(:,2,1) = Data.data(end,1,:);
        
        
        for i = 1:length(Signal)
            
            VarName = strrep(Signal{i},' ','');
            eval([VarName, ' = SelectVar2(Data,Sensor{Sub}(1),Signal(i));']);
                        
            temp = eval(VarName);
            %Filtering
            eval([VarName 'F = [temp(:,1,:),filtfilt(b,a,temp(:,2:end,:))];']);
        end
        
        %EMG = = ReadEMGworksWindows(FilePath,1,'EMG', deltaT);
        Cte = ones(size(eval(VarName),1),size(eval(VarName),2), size(eval(VarName),3));
        
        % FSR data
        Ftoe =  SelectVar2(Data,11,'ACC Y');
        
        
        tic
        
        %% -- Resultants
        
        % Delsys
        ACC = [ACCX(:,1,:), sqrt(ACCX(:,2:end,:).^2 + ACCY(:,2:end,:).^2 + ACCZ(:,2:end,:).^2)];
        ACCF = [ACCX(:,1,:), sqrt(ACCXF(:,2:end,:).^2 + ACCYF(:,2:end,:).^2 + ACCZF(:,2:end,:).^2)];
        
        Gyro = [GyroX(:,1,:), sqrt(GyroX(:,2:end,:).^2 + GyroY(:,2:end,:).^2 + GyroZ(:,2:end,:).^2)];
        GyroF = [GyroX(:,1,:), sqrt(GyroXF(:,2:end,:).^2 + GyroYF(:,2:end,:).^2 + GyroZF(:,2:end,:).^2)];
        
        %         Mag = [MagX(:,1,:), sqrt(MagX(:,2:end,:).^2 + MagY(:,2:end,:).^2 + MagZ(:,2:end,:).^2)];
        %         MagF = [MagX(:,1,:), sqrt(MagXF(:,2:end,:).^2 + MagYF(:,2:end,:).^2 + MagZF(:,2:end,:).^2)];
        
        %% ARRUMAR RESHAPE
       
        
        %[NewInstant] = ReshapeInstants(delta, instant,[Folder,Name]);
        [TO, HS] = ReshapeInstants2(delta, instant,Files);
        
        TO = TO + 36e-3; %Delay Delsys ACC
        HS = HS  + 36e-3;
        
        %% Linear combination of different variables apllied in one trial
        
        % --- Building inputs for Orthogonal Least Squares Algorithm
        % --- (ols.m) implemented by Renato Naville Watanabe
        
        % -- Initializing variables for linear combination
        
        pp = NaN(size(Data.data,1),length(Var),size(Data.data,3));
        yy = NaN(size(Data.data,1),size(Data.data,3));
        
        first = NaN(1,size(Data.data,3));
        last = NaN(1,size(Data.data,3));
        
        for j = 1 : size(Data.data,3)
            % for i = 1 : length(Sensors)
                first(j) = min([HS(j,:),TO(j,:)],[],2);
                last(j) = max([HS(j,:),TO(j,:)],[],2);
                
                tempTO = (TO(j,:));
                tempTO(isnan(tempTO))=[];
                if strcmp(Win(w),'Gauss')
                    stimulWin = sum(exp(-0.5*((Data.data(:,1,j) - (tempTO - sd))/(sd/3)).^2),2);
                end
                
                if strcmp(Win(w),'Rect')
                    stimulWin = zeros(length(Data.data),1);
                    for kk = 1: size(TO,2)
                        stimulWin(floor((TO(j,kk)-2*sd)*Fs):floor(TO(j,kk)*Fs)) = 1;
                    end
                end
                yy(:,j) = stimulWin;
                
                for jj = 1 : length (Var)
                    pp(:,jj,:) = eval(Var{jj});
                    
                end
            % end
        end
    end
    p(:,:,(Sub-1)*(length(Sensors)*size(Data.data,3))+1: Sub*(length(Sensors)*size(Data.data,3))) = pp;
    y(:,(Sub-1)*(length(Sensors)*size(Data.data,3))+1: Sub*(length(Sensors)*size(Data.data,3))) = yy;
    ForceY((Sub-1)*60*FsFP +1:Sub*60*FsFP,:) = Fy;
    
    ToeOff = [ToeOff;TO];
    HeelStrike = [HeelStrike;HS];
    timeACC(:,(Sub-1)*(length(Sensors)*size(Data.data,3))+1: Sub*(length(Sensors)*size(Data.data,3)),1) = ACC(:,1,:);
    
end
%% --- TESTANDO AS JANELAS DE ESTÍMULO
%         figure;
%         plot(Fy(:,1),Fy(:,3),'k')
%         hold on
%         plot(timeACC,1100*y)

%% --  -- Select trials for trainning and test
[indTr,indTs] = PartData(ToeOff,length(ToeOff)*8/11);
% % % % indTr = [1;3;4;10;15;20;21;22;33;34;35;36];
% % % % indTs = [2;5;6;7;8;9;11;12;13;14;16;17;18;19;23;24;25;26;27;28;29;30;31;32;37;38;39;40;41;42;43;44];

%% --- Combinatorial Analysis

n=0; %t=0;
for pct = 0: 0.01 : 1
    n = n+1;
    disp(pct)
    [ResultsCombinatorics] = combinatorics1(Var,Sensors,p,y,pct,ForceY,ToeOff,timeACC,indTr,indTs);
    
    % RT(t+1:t+length(ResultsTrials)) = ResultsTrials;
    RC((n-1)*TotalComb + 1 : n*TotalComb) = ResultsCombinatorics;
    
    % t = ;
    
    
    
end
% save(['RTdelay2_',Name,'.mat'],'RT')
save(['RCdelayMag_',Name,'.mat'],'RC')

% RT = struct('Trial',{0},'Features',{0},'Locs',{0},'Threshold',{0},...
%     'TP',{0},'FP',{0},'TN',{0},'FN',{0},'beta',{0});
% RT = repmat(RT,TotalComb*101*32,1);

% RC = struct('TrialsTr',{0},'Features',{0},'Threshold',{0},'TP',{0},...
%     'FP',{0},'TN',{0},'FN',{0},'beta',{0});
% RC = repmat(RC,TotalComb*101,1);




toc

%% -- Plots

% % figure;
% % subplot(3,1,1); plot(Fy(:,1), Fy(:,2:end));
% % legend({'1','2','3','4','5','6','7'})
% %
% % %Delsys
% % subplot(3,1,2); plot(ACC(:,1), ACC(:,ShankL),ACCF(:,1), ACCF(:,ShankL)); ylabel('Shank L')
% % title('Resultant'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACC(:,ShankR),ACCF(:,1), ACCF(:,ShankR)); ylabel('Shank R')
% % % Kistler
% % %subplot(3,1,2); plot(ACC1(:,1),ACC1(:,ShankL)); ylabel('Shank L')
% % %title('Resultant')
% % %subplot(3,1,3); plot(ACC2(:,1),ACC2(:,ShankL)); ylabel('Shank R')
% % ylim([-2 2]);
% % PlotInstants( instant, File )
% %
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% %
% % subplot(3,1,2); plot(ACC(:,1), ACCX(:,ShankL),ACC(:,1), ACCXF(:,ShankL)); ylabel('Shank L')
% % title('X'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCX(:,ShankR),ACCF(:,1), ACCXF(:,ShankR)); ylabel('Shank R')
% % ylim([min(ACCX(:,ShankL))-1, max(ACCX(:,ShankL))+1]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(ACC(:,1), ACCY(:,ShankL), ACC(:,1), ACCYF(:,ShankL));  ylabel('Shank L')
% % title('Y'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCY(:,ShankR), ACCF(:,1), ACCYF(:,ShankR)); ylabel('Shank R')
% % ylim([min(ACCY(:,ShankL))-1, max(ACCY(:,ShankL))+1]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(ACC(:,1), ACCZ(:,ShankL), ACC(:,1), ACCZF(:,ShankL)); ylabel('Shank L');
% % title('Z'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCZ(:,ShankR), ACCF(:,1), ACCZF(:,ShankR));ylabel('Shank R');
% % ylim([min(ACCZ(:,ShankL))-1, max(ACCZ(:,ShankL))+1]);
% % PlotInstants( instant, File )
