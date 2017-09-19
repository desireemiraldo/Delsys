%% Gait Analysis (w/ accelerometer)
clear all; clc; close all

%% -- About the files
% % Path = '.\Dropfoot\';
% % % -- Settings
% % Name = 'DCM\Descalco25_Dual_170831_1';
% % csv = '-Delsys 1.csv';

paths = {'.\Dropfoot\','.\Dropfoot\'};
folders = {'DCM\','DCM\'};
names = {'Descalco25_Dual_170831_1','Tenis25_Dual_170831_1'};
ext = {'-Delsys 1.csv','-Delsys 1.csv'};
shanksL = [2,2];
shanksR = [1,1];
sensor = {{'ShankR'},{'ShankR'}};

Win = {'Gauss'}; %,'Rect'};

%% -- Initializing some variables

% -- Sample Frequencies
Fs = 148.39; % Delsys sensor
FsFP = 300;  % force plates

% --
ChannelType = 'AUX';
Signal = {'IM ACC Pitch', 'IM ACC Roll', 'IM ACC Yaw',...
    'IM GYR Pitch', 'IM GYR Roll', 'IM GYR Yaw',...
    'IM MAG Pitch', 'IM MAG Roll', 'IM MAG Yaw',};

Var = {'ACCF','ACCPitchF','ACCRollF','ACCYawF',...
    'GYRF','GYRPitchF','GYRRollF','GYRYawF',...
    'MAGF','MAGPitchF','MAGRollF','MAGYawF','Cte'};%,'AccF.^2'}; %Used in linear combination

TotalComb = 0;
for k = 1: length(Var)
    b = nchoosek(length(Var),k);
    TotalComb = TotalComb + b;
end

% -- Filter
[t,Wn] = buttord(10/(Fs/2),20/(Fs/2),1,60);
[b,a] = butter(t, Wn);

% -- Standard deviation to be used in Linear Combination
sd = 50e-3;

% --
winPts = 202;

% --

% RT = struct('Trial',{0},'Features',{0},'Locs',{0},'Threshold',{0},...
%     'TP',{0},'FP',{0},'TN',{0},'FN',{0},'beta',{0});
% RT = repmat(RT,TotalComb*101*32,1);

RC = struct('TrialsTr',{0},'Features',{0},'Threshold',{0},'TP',{0},...
    'FP',{0},'TN',{0},'FN',{0},'beta',{0});
RC = repmat(RC,TotalComb*101,1);

p = NaN(2*winPts,length(Var),2*22);
y = NaN(2*winPts,2*22);
ForceY = NaN(2*60*FsFP,8);
timeACC = NaN(size(y));
ToeOff = [];
HeelStrike = [];

for Sub = 1: length(names)
    for w = 1: length(Win)
        
        Path = paths{Sub};
        Folder = folders{Sub};
        Name = names {Sub};
        csv = ext{Sub};
        ShankR = shanksR(Sub); %ShankL = shanksL(Sub);
        Sensors = sensor{Sub};
        
        
        FilePath = [Path,Folder,Name];
        
        %% Loading Data
        
        [Trials,deltaT,TrLabel] = ReadDelsysWindows([FilePath,csv], winPts);
        
        instant = importdata('Instantes_gait1.txt',';');
        
        %[NewInstant] = ReshapeInstants(deltaT, instant,[Folder,Name]);
        [TO, HS] = ReshapeInstants1(deltaT, instant,[Folder,Name]);
        
        TO = TO + 70e-3; %Delay Delsys ACC
        HS = HS  + 70e-3;
        
        tic
        
        for i = 1:length(Signal)
            VarName = strrep(strrep(Signal{i},'IM ',''),' ','');
            eval([VarName '= SelectVar(Trials, ChannelType, Signal(i));']);
            temp = eval(VarName);
            %Filtering
            eval([VarName 'F = [temp(:,1,:),filtfilt(b,a,temp(:,2:end,:))];']);
        end
        
        %EMG = SelectVar(Trials, 'IMEMG', 'EMG');
        
        Cte = ones(size(eval(VarName),1),size(eval(VarName),2), size(eval(VarName),3));
        
        % Cortex data
        Forces = importdata([FilePath,'.forces']);
        Fy = (Forces.data(:,1) -1)/FsFP;
        for i =1 : length(Forces.colheaders)
            if strcmp(Forces.colheaders{i}(1:2), 'FY')
                Fy = [Fy, Forces.data(:,i)];
            end
        end
        
        
        Fy(:,6:8) = Fy(:,6:8)/10;
        
        %% -- Resultants
        
        % Delsys
        ACC = [ACCPitch(:,1,:), sqrt(ACCPitch(:,2:end,:).^2 + ACCRoll(:,2:end,:).^2 + ACCYaw(:,2:end,:).^2)];
        ACCF = [ACCPitch(:,1,:), sqrt(ACCPitchF(:,2:end,:).^2 + ACCRollF(:,2:end,:).^2 + ACCYawF(:,2:end,:).^2)];
        
        GYR = [GYRPitch(:,1,:), sqrt(GYRPitch(:,2:end,:).^2 + GYRRoll(:,2:end,:).^2 + GYRYaw(:,2:end,:).^2)];
        GYRF = [GYRPitch(:,1,:), sqrt(GYRPitchF(:,2:end,:).^2 + GYRRollF(:,2:end,:).^2 + GYRYawF(:,2:end,:).^2)];
        
        MAG = [MAGPitch(:,1,:), sqrt(MAGPitch(:,2:end,:).^2 + MAGRoll(:,2:end,:).^2 + MAGYaw(:,2:end,:).^2)];
        MAGF = [MAGPitch(:,1,:), sqrt(MAGPitchF(:,2:end,:).^2 + MAGRollF(:,2:end,:).^2 + MAGYawF(:,2:end,:).^2)];
        
        
        %% Linear combination of different variables apllied in one trial
        
        % --- Building inputs for Orthogonal Least Squares Algorithm
        % --- (ols.m) implemented by Renato Naville Watanabe
        
        % -- Initializing variables for linear combination
        
        pp = NaN(length(ACC),length(Var),(length(Sensors)*size(Trials,3)));
        yy = NaN(length(ACC),(length(Sensors)*size(Trials,3)));
        
        first = NaN(1,(length(Sensors)*size(Trials,3)));
        last = NaN(1,(length(Sensors)*size(Trials,3)));
        
        for j = 1 : size(Trials,3)
            for i = 1 : length(Sensors)
                first(j) = min([HS(j,:),TO(j,:)],[],2);
                last(j) = max([HS(j,:),TO(j,:)],[],2);
                
                tempTO = (TO(j,:));
                tempTO(isnan(tempTO))=[];
                if strcmp(Win(w),'Gauss')
                    stimulWin = sum(exp(-0.5*((ACC(:,1,j) - (tempTO - sd))/(sd/3)).^2),2);
                end
                
                if strcmp(Win(w),'Rect')
                    stimulWin = zeros(length(ACC),1);
                    for kk = 1: size(TO,2)
                        stimulWin(floor((TO(j,kk)-2*sd)*Fs):floor(TO(j,kk)*Fs)) = 1;
                    end
                end
                yy(:,j) = stimulWin;
                
                for jj = 1 : length (Var)
                    pp(:,jj,:) = eval([Var{jj},'(:,eval(Sensors{i})+1,:)']);
                    
                end
            end
        end
    end
    p(:,:,(Sub-1)*(length(Sensors)*size(Trials,3))+1: Sub*(length(Sensors)*size(Trials,3))) = pp;
    y(:,(Sub-1)*(length(Sensors)*size(Trials,3))+1: Sub*(length(Sensors)*size(Trials,3))) = yy;
    ForceY((Sub-1)*60*FsFP +1:Sub*60*FsFP,:) = Fy;
    
    ToeOff = [ToeOff;TO];
    HeelStrike = [HeelStrike;HS];
    timeACC(:,(Sub-1)*(length(Sensors)*size(Trials,3))+1: Sub*(length(Sensors)*size(Trials,3)),1) = ACC(:,1,:);

end
%% --- TESTANDO AS JANELAS DE ESTÍMULO
%         figure;
%         plot(Fy(:,1),Fy(:,3),'k')
%         hold on
%         plot(timeACC,1100*y)

%% --  -- Select trials for trainning and test
[indTr,indTs] = PartData(ToeOff,length(ToeOff)*8/11);


indTr = [1;3;4;10;15;20;21;22;33;34;35;36];
indTs = [2;5;6;7;8;9;11;12;13;14;16;17;18;19;23;24;25;26;27;28;29;30;31;32;37;38;39;40;41;42;43;44];
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
save(['RCdelayMAG_',Name,'.mat'],'RC')

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
% % subplot(3,1,2); plot(ACC(:,1), ACCPitch(:,ShankL),ACC(:,1), ACCPitchF(:,ShankL)); ylabel('Shank L')
% % title('Pitch'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCPitch(:,ShankR),ACCF(:,1), ACCPitchF(:,ShankR)); ylabel('Shank R')
% % ylim([min(ACCPitch(:,ShankL))-1, max(ACCPitch(:,ShankL))+1]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(ACC(:,1), ACCRoll(:,ShankL), ACC(:,1), ACCRollF(:,ShankL));  ylabel('Shank L')
% % title('Roll'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCRoll(:,ShankR), ACCF(:,1), ACCRollF(:,ShankR)); ylabel('Shank R')
% % ylim([min(ACCRoll(:,ShankL))-1, max(ACCRoll(:,ShankL))+1]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(ACC(:,1), ACCYaw(:,ShankL), ACC(:,1), ACCYawF(:,ShankL)); ylabel('Shank L');
% % title('Yaw'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCYaw(:,ShankR), ACCF(:,1), ACCYawF(:,ShankR));ylabel('Shank R');
% % ylim([min(ACCYaw(:,ShankL))-1, max(ACCYaw(:,ShankL))+1]);
% % PlotInstants( instant, File )
