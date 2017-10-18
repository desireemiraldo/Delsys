%% Gait Analysis (w/ accelerometer)
clear all; clc; close all

%% -- About the files

Path = '.\Piloto\';

Folder = {'RNW\'};

Files = {'RNW_Cal�ado_Confortavel_Rep_1.1',...
    'RNW_Cal�ado_Confortavel_Rep_2.2',...
    'RNW_Cal�ado_Confortavel_Rep_3.3',...
    'RNW_Cal�ado_Confortavel_Rep_4.4',...
    'RNW_Cal�ado_Confortavel_Rep_5.5',...
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
    'GyroF','GyroXF','GyroYF','GyroZF', 'Cte'};...
    %'MagF','MagXF','MagYF','MagZF'};%,'AccF.^2'}; %Used in linear combination

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

delay = 0;

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
        %% Loading data
        
        instant = importdata([Path,Folder{Sub},'Instantes.txt'],';');
        
        %% APAGAR DEPOIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for x = 2:11
            instant.textdata(x,1) = strrep(instant.textdata(x,1),'Calcado','Cal�ado');
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %[Data,numWin] = BuildingTrials2([Path,Folder{Sub}],Files, 3, deltaT);
        % save('RNW_data.mat','Data');
        load RNW_data.mat
        numWin = [0;7;14;21;28;34;41;48;54;61;68];
        
        time = Data.data(:,1,:);
        delta(:,1,1) = time(1,1,:);
        delta(:,2,1) = time(end,1,:);
        
        [TO, HS] = ReshapeInstants2(delta, instant,Files,numWin);
        
        TO = TO + delay; % equipaments delay
        HS = HS  + delay;
        
        
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
        
        % Mag = [MagX(:,1,:), sqrt(MagX(:,2:end,:).^2 + MagY(:,2:end,:).^2 + MagZ(:,2:end,:).^2)];
        % MagF = [MagX(:,1,:), sqrt(MagXF(:,2:end,:).^2 + MagYF(:,2:end,:).^2 + MagZF(:,2:end,:).^2)];
        
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
        end
        for jj = 1 : length (Var)
            pp(:,jj,:) = eval([Var{jj},'(:,2,:)']);
            
        end
        % end
    end
end

%% --- TESTANDO AS JANELAS DE EST�MULO
%         figure;
%         plot(Ftoe(:,1),Ftoe(:,2),'k')
%         hold on
%         plot(time,1100*y)

%% --  -- Select trials for trainning and test
[indTr,indTs] = PartData(TO,floor(length(TO)*0.35));

%% --- Combinatorial Analysis

n=0; %t=0;
for pct = 0: 0.01 : 1
    n = n+1;
    disp(pct)
    [ResultsCombinatorics] = combinatorics2(Var,1,pp,yy,pct,TO,Ftoe,delay,time,indTr,indTs);
    
    % RT(t+1:t+length(ResultsTrials)) = ResultsTrials;
    RC((n-1)*TotalComb + 1 : n*TotalComb) = ResultsCombinatorics;
   
end
Name = strrep(Folder{Sub},'\','sensor1');
% save([Name,'.mat'],'RT')
save([Name,'.mat'],'RC')

% RT = struct('Trial',{0},'Features',{0},'Locs',{0},'Threshold',{0},...
%     'TP',{0},'FP',{0},'TN',{0},'FN',{0},'beta',{0});
% RT = repmat(RT,TotalComb*101*32,1);

% RC = struct('TrialsTr',{0},'Features',{0},'Threshold',{0},'TP',{0},...
%     'FP',{0},'TN',{0},'FN',{0},'beta',{0});
% RC = repmat(RC,TotalComb*101,1);




toc

%% -- Plots

% % figure;
% % subplot(3,1,1); plot(Ftoe(:,1), Ftoe(:,2:end));
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
% % plot(Ftoe(:,1), Ftoe(:,2:end)); legend({'1','2','3','4','5','6','7'})
% %
% % subplot(3,1,2); plot(ACC(:,1), ACCX(:,ShankL),ACC(:,1), ACCXF(:,ShankL)); ylabel('Shank L')
% % title('X'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCX(:,ShankR),ACCF(:,1), ACCXF(:,ShankR)); ylabel('Shank R')
% % ylim([min(ACCX(:,ShankL))-1, max(ACCX(:,ShankL))+1]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Ftoe(:,1), Ftoe(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(ACC(:,1), ACCY(:,ShankL), ACC(:,1), ACCYF(:,ShankL));  ylabel('Shank L')
% % title('Y'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCY(:,ShankR), ACCF(:,1), ACCYF(:,ShankR)); ylabel('Shank R')
% % ylim([min(ACCY(:,ShankL))-1, max(ACCY(:,ShankL))+1]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Ftoe(:,1), Ftoe(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(ACC(:,1), ACCZ(:,ShankL), ACC(:,1), ACCZF(:,ShankL)); ylabel('Shank L');
% % title('Z'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCZ(:,ShankR), ACCF(:,1), ACCZF(:,ShankR));ylabel('Shank R');
% % ylim([min(ACCZ(:,ShankL))-1, max(ACCZ(:,ShankL))+1]);
% % PlotInstants( instant, File )
