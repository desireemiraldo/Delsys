%% Gait Analysis (w/ accelerometer)
clear all; clc; close all

%% -- About the files

Path = '.\Data\';

Folder = {'S001\','S001\','S002\','S003\','S004\','S005\','S006\','S007\','S008\','S009\','S010\'};

right = {[1,2,11],[1,2,11],[4,3,11],[3,4,11],[3,4,11],[3,4,11],[3,4,11],[3,4,11],[3,4,11],[3,4,11],[3,4,11]};

left = {[1,2,11],[1,2,11],[6,5,12],[5,6,12],[5,6,12],[5,6,12],[5,6,12],[5,6,12],[5,6,12],[5,6,12],[5,6,12]};

Win = {'Gauss'};%,'Rect'};

Name = char(inputdlg('Enter the name for result files:',...
    'RESULT FILES - NAME', [1 50]));
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
Fs = 148.1481;
[t,Wn] = buttord(10/(Fs/2),20/(Fs/2),1,60);
[b,a] = butter(t, Wn);

% -- Standard deviation to be used in Linear Combination
sd = 50e-3;

% --
deltaT = 2.75; %seconds/window

delay = 0;

% --
Trials = 0; TempTO = 0;
for w = 1: length(Win)
    for Sub = 1: 2%length(Folder)
        Files = dir([Path,Folder{Sub},'*Rep*.csv']);
        Files = {(Files(:).name)}';
        shank = {'right', 'left'};
        side = randperm(2,1);
        
        Sensor = eval(shank{side});
        
        for numSensor = 1: 1 % length(Sensor{Sub})
            
            %% Loading data
            
            instant = importdata([Path,Folder{Sub},'Instantes_gait.txt'],'\t');
            
            
            % % % % % % % % % %             load RNW_data.mat
            % % % % % % % % % %             numWin = [0;7;14;21;28;34;41;48;54;61;68];
            load('.\Piloto\RNW\.matData\RNWvelocities_data.mat');
            load('.\Piloto\RNW\.matData\RNWvelocities_numWin.mat');
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %             [Data,numWin] = BuildingTrials([Path,Folder{Sub}],Files, 3, deltaT);
            %             save([Name,'_data.mat'],'Data');
            %             save([Name,'_numWin.mat'],'numWin');
            %             movefile ([Name,'_data.mat'],[Path,Folder{Sub},'.matData'])
            
            time = Data.data(:,1,:);
            delta = nan(size(time,3),2);
            delta(:,1) = time(1,1,:);
            delta(:,2) = time(end,1,:);
            
            [newTO, newHS] = ReshapeInstants(delta, instant,Files,numWin);
            
            newTO = newTO + delay; % equipaments delay
            newHS = newHS  + delay;
            
            
            for i = 1:length(Signal)
                
                VarName = strrep(Signal{i},' ','');
                eval([VarName, ' = SelectVar(Data,Sensor{Sub}(numSensor),Signal(i));']);
                
                temp = eval(VarName);
                %Filtering
                eval([VarName 'F = [temp(:,1,:),filtfilt(b,a,temp(:,2:end,:))];']);
            end
            
            %EMG = = ReadEMGworksWindows(FilePath,1,'EMG', deltaT);
            Cte = ones(size(eval(VarName),1),size(eval(VarName),2), size(eval(VarName),3));
            
            % FSR data
            newFtoe =  SelectVar(Data,11,'ACC Y');
            
            
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
            
            newP = NaN(size(Data.data,1),length(Var),size(Data.data,3));
            newY = NaN(size(Data.data,1),size(Data.data,3));
            
            first = NaN(1,size(Data.data,3));
            last = NaN(1,size(Data.data,3));
            
            for j = 1 : size(Data.data,3)
                % for i = 1 : length(Sensors)
                first(j) = min([newHS(j,:),newTO(j,:)],[],2);
                last(j) = max([newHS(j,:),newTO(j,:)],[],2);
                
                tempTO = (newTO(j,:));
                tempTO(isnan(tempTO))=[];
                if strcmp(Win(w),'Gauss')
                    stimulWin = sum(exp(-0.5*((Data.data(:,1,j) - (tempTO - sd))/(sd/3)).^2),2);
                end
                
                if strcmp(Win(w),'Rect')
                    stimulWin = zeros(length(Data.data),1);
                    for kk = 1: size(newTO,2)
                        if isfinite(newTO(j,kk))
                            ind1 = floor(((newTO(j,kk)-delta(j,1)) -2*sd)*Fs);
                            ind2 =  floor((newTO(j,kk)-delta(j,1))*Fs);
                            
                            if ind1<=0
                                ind1 = 1;
                            end
                            
                            stimulWin(ind1:ind2,1) = 1;
                        end
                    end
                end
                newY(:,j) = stimulWin;
            end
            for jj = 1 : length (Var)
                newP(:,jj,:) = eval([Var{jj},'(:,2,:)']);
                
            end
            % end
        end
        
        p = NaN(size(Data.data,1),length(Var),size(Data.data,3) + Trials);
        y = NaN(size(Data.data,1),size(Data.data,3) + Trials);
        TO = NaN (size(Data.data,3) + Trials, max([size(newTO,2),size(TempTO,2)]));
        Ftoe = NaN(size(Data.data,1),2,size(Data.data,3) + Trials);
        
        %% Ftoe
        if Trials ~=0
            p(:,:,1:Trials) = TempP;
            p(:,:,Trials+1:end) = newP;
            y(:,1:Trials) = TempY;
            y(:,Trials+1:end) = newY;
            TO(1:Trials,1:size(TempTO,2)) = TempTO;
            TO(Trials+1:end,1:size(newTO,2)) = newTO;
            Ftoe(:,:,1:Trials) = TempFtoe;
            Ftoe(:,:,Trials+1:end) = newFtoe;
        else
            p = newP;
            y = newY;
            TO = newTO;
            Ftoe = newFtoe;
        end
        TempP = p;
        TempY = y;
        TempTO = TO;
        TempFtoe = Ftoe;
        Trials = Trials + size(Data.data,3);
        clear newP; clear newY; clear newTO; clear newFtoe;
        
    end
end
time = Ftoe(:,1,:);

%% --- TESTANDO AS JANELAS DE ESTÍMULO
% figure;
% plot(Ftoe(:,1),Ftoe(:,2),'k')
% hold on
% plot(time,1100*y)

%% --  -- Select trials for trainning and test
if numSensor == 1
    [indTr,indTs] = PartData(TO,floor(length(TO)*0.35));
end
%             indTr = [1;2;3;4;5;6;7;10;11;13;16;17;18;19;22;24;25;26;28;31;32;33;34;35;37;38;39;41;42;44;45;46;47;48;49;51;53;55;56;58;62;64;65;66;67];
%             indTs = [8;9;12;14;15;20;21;23;27;29;30;36;40;43;50;52;54;57;59;60;61;63;68];

%% --- Combinatorial Analysis

n=0; %t=0;
RC = struct('TrialsTr',{0},'Features',{0},'Threshold',{0},'TP',{0},...
    'FP',{0},'TN',{0},'FN',{0},'beta',{0});
RC = repmat(RC,TotalComb*101,1);

% RT = struct('Trial',{0},'Features',{0},'Locs',{0},'Threshold',{0},...
%     'TP',{0},'FP',{0},'TN',{0},'FN',{0},'beta',{0});
% RT = repmat(RT,TotalComb*101*length(indTr),1);

for pct = -0.5: 0.01 : 1
    n = n+1;
    disp(pct)
    [ResultsCombinatorics] = combinatorics(Var,1,p,y,pct,TO,Ftoe,delay,time,indTr,indTs);
    
    % RT(t+1:t+length(ResultsTrials)) = ResultsTrials;
    RC((n-1)*TotalComb + 1 : n*TotalComb) = ResultsCombinatorics;
    
end
nameRes = [Name,'sensor',num2str(Sensor{Sub}(numSensor)),'.mat'];
% save([Name,'.mat'],'RT')
save(nameRes,'RC')

movefile (nameRes,'Resultados')

clear RC; %clear RT

toc
%         end
%     end
% end
