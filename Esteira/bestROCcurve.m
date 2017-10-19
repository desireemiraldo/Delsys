% clc; clear all; close all
Path1 = 'C:\Users\BMClab\Downloads\Desiree\Drop Foot Project\Identificacao de evento\Processamento\Resultados\';
load('RNWsensor1.mat')

indTr = [1;2;3;4;5;6;7;10;11;13;16;17;18;19;22;24;25;26;28;31;32;33;34;35;37;38;39;41;42;44;45;46;47;48;49;51;53;55;56;58;62;64;65;66;67];
indTs = [8;9;12;14;15;20;21;23;27;29;30;36;40;43;50;52;54;57;59;60;61;63;68];

%load([Path1,'MelhoresCombinacoes_9var.mat'])
len = 511;

% load('RCdelay_13var_Desiree.mat');
% load('MelhoresCombinacoes_13var.mat');
% len = 8191;

i = Order(1).Index;

[Sensib,Specif] = ROCcurve(RC);

figure(1); cmap = colormap(jet(len));
xlim([0 1]); ylim = ([0 1]);
xlabel('Taxa de falso positivo');
ylabel('Taxa de verdadeiro positivo');
title ('Curva ROC da melhor combinação');
hold on

ind = i:len:len*101;
plot(1-Specif(ind),Sensib(ind),'color',cmap(i,:))%,'.','MarkerSize',12);

line([0 1], [0 1],'LineStyle','--','color','k','linewidth',2);


%% melhor threshold

[mini,pct] = min(sqrt((1-Specif(ind)).^2 + (1-Sensib(ind)).^2));

pct = pct/100;

% % indTr = [1;3;4;10;15;20;21;22;33;34;35;36];
% % indTs = [2;5;6;7;8;9;11;12;13;14;16;17;18;19;23;24;25;26;27;28;29;30;31;32;37;38;39;40;41;42;43;44];

%%

LinearCombination = NaN(size(time));


%% -- Dividing forces
trialsNumber = length(ToeOff);
ind = 0;
while rem(size(ForceY(1:end-ind,1)),trialsNumber)~=0
    ind = ind+1;
end

[row,col] = size(ForceY(1:end-ind,:));
len = row/trialsNumber;

Fy = zeros(len,col,trialsNumber);

for i = 1: trialsNumber
    Fy(:,:,i) = ForceY((i-1)*len +1 : i*len,:);
end


%% -- Select trials for trainning and test
% [indTr,indTs] = PartData(ToeOff,16);

%% --

n = 0;

Features = [2 3 4 5 6 7 8 9 10];


pTr = p(:,Features,indTr);
yTr = y(:,indTr);

pTs = p(:,Features,indTs);

beta = ols(pTr,yTr); % Coefs for linear combination
betaM = mean(beta,2);

% --- Applying Linear Combination

for j = 1:length(indTs)
    for i = 1 : length(Sensors)
        
        index = indTs(j)+(i-1)*length(indTs);
        
        n = n+1; %disp(n)
        
        LinearCombination(:,index) = pTs(:,:,(j)+(i-1)*length(indTs))*betaM;
        
        %% Plot Fy
        kg = 54;
        
        figure(index)
        subplot(2,1,1); plot(Fy(:,1,index)+70e-3, Fy(:,3,index)/kg,'linewidth',1.5,'Color',[0.3 0.3 0.3]);
        hold on;
        plot(timeACC(:,index),20*y(:,index),'r','linewidth',1.5)
        xlim([timeACC(1,index),timeACC(end,index)]);
        xticklabels({})
        ylabel('Força [kg/N]','fontsize',12)
        %title([' Trial ', num2str(index)]);
        
        subplot(2,1,2); plot(timeACC(:,index),LinearCombination(:,index),'linewidth',1.5);
        ylim([min(LinearCombination(:,index))*1.1 max(LinearCombination(:,index))*1.1]);
        xlim([timeACC(1,index),timeACC(end,index)]);
        ylabel('Combinação Linear','fontsize',12)
        xlabel('Tempo [s]','fontsize',12)
        
        %  --- Checking the combination's quality
        threshold = (max(LinearCombination(:,index)))*pct;
        [pks,locs] = findpeaks(LinearCombination(:,index),timeACC(:,index),'MinPeakHeight',threshold);
        
        for z = 1:length(locs)
            Line = line([locs(z) locs(z)], [-0.15 1.5],'Linewidth',1.5,'Linestyle','--','Color',[0 0 0]);
            set(Line,'Clipping','off')
        end
        
    end
end



