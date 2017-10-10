paths = {'.\Piloto\','.\Piloto\','.\Piloto\','.\Piloto\','.\Piloto\','.\Piloto\'};
folders = {'RNW\','RNW\','RNW\','RNW\','RNW\','RNW\'};
names = {'Piloto_profRenato__Plot_and_Store_Rep_1.5',...
    'Piloto_profRenato__Plot_and_Store_Rep_2.6',...
    'Piloto_profRenato__Plot_and_Store_Rep_2.1',...
    'Piloto_profRenato_Tenis_Plot_and_Store_Rep_1.2',...
    'Piloto_profRenato_Tenis_Plot_and_Store_Rep_1.3',...
    'Piloto_profRenato_Tenis_Plot_and_Store_Rep_2.4'};
ext = {'.xls','.xls','.xls','.xls','.xls','.xls'};



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

% RT = struct('Trial',{0},'Features',{0},'Locs',{0},'Threshold',{0},...
%     'TP',{0},'FP',{0},'TN',{0},'FN',{0},'beta',{0});
% RT = repmat(RT,TotalComb*101*32,1);


ToeOff = [];
HeelStrike = [];

for Sub = 1: length(names)
        
        Path = paths{Sub};
        Folder = folders{Sub};
        Name = names {Sub};
        xls = ext{Sub};       
        
        FilePath = [Path,Folder,Name,xls];
        
        %% Loading Data
              
        for i = 1:length(Signal)
            VarName = strrep(Signal{i},' ','');
            eval([VarName '=  ReadEMGworksWindows(FilePath,1,Signal(i), deltaT);']);
            %% Sensor = numero do sensor na torre delsys
            
            temp = eval(VarName);
            %Filtering
            eval([VarName 'F = [temp(:,1,:),filtfilt(b,a,temp(:,2:end,:))];']);
        end
        
        EMG =  ReadEMGworksWindows(FilePath,1,'EMG', deltaT);
       
        % FSR data
        Fheel = ReadEMGworksWindows(FilePath,11,'ACC X', deltaT);
        Ftoe =ReadEMGworksWindows(FilePath,11,'ACC Y', deltaT);
        


        %% -- Resultants
        
        % Delsys
        ACC = [ACCX(:,1,:), sqrt(ACCX(:,2:end,:).^2 + ACCY(:,2:end,:).^2 + ACCYaw(:,2:end,:).^2)];
        ACCF = [ACCX(:,1,:), sqrt(ACCXF(:,2:end,:).^2 + ACCYF(:,2:end,:).^2 + ACCZF(:,2:end,:).^2)];
        
        Gyro = [GyroX(:,1,:), sqrt(GyroX(:,2:end,:).^2 + GyroY(:,2:end,:).^2 + GyroZ(:,2:end,:).^2)];
        GyroF = [GyroX(:,1,:), sqrt(GyroXF(:,2:end,:).^2 + GyroYF(:,2:end,:).^2 + GyroZF(:,2:end,:).^2)];
        
%         Mag = [MagX(:,1,:), sqrt(MagX(:,2:end,:).^2 + MagY(:,2:end,:).^2 + MagZ(:,2:end,:).^2)];
%         MagF = [MagX(:,1,:), sqrt(MagXF(:,2:end,:).^2 + MagYF(:,2:end,:).^2 + MagZF(:,2:end,:).^2)];

figure;
subplot(4,1,1)
plot(ACCF(:,1), ACCF(:,2),Fheel(:,1),Fheel(:,2),Ftoe(:,1),Ftoe(:,2))
xlim([0 20]); %ylim([0 10]);
legend('ACC','heel',   'toe')
title('ACC')

subplot(4,1,2)
plot(ACCXF(:,1), ACCXF(:,2))

subplot(4,1,3)
plot(ACCYF(:,1), ACCYF(:,2))

subplot(4,1,4)
plot(ACCZF(:,1), ACCZF(:,2))

figure;
subplot(4,1,1)
plot(GyroF(:,1), GyroF(:,2),EMG(:,1),EMG(:,2))
xlim([0 20]); %ylim([0 10]);
legend('Gyro','EMG')
title('Gyro')

subplot(4,1,2)
plot(GyroXF(:,1), GyroXF(:,2))

subplot(4,1,3)
plot(GyroYF(:,1), GyroYF(:,2))

subplot(4,1,4)
plot(GyroZF(:,1), GyroZF(:,2))

end