clc; clear all; close all;

try
    parpool('local', 1);
catch
    disp('Parpool failed to launch, continuing on single core')
end



rundate = char(datetime('today'));

NT = 8;
SNRsS = logspace(0,2,NT);
SNRsT = logspace(1,6,NT);

NsnrS = size(SNRsS,2);
NsnrT = size(SNRsT,2);

Nnoise = 96;

% N = 4; %Number of Realizations
Ns = ceil(linspace(2400,2400,NT));

%% Initialize and Run

corr_corr_1DLimits_S = zeros(1,NsnrS);
corr_uncorr_1DLimits_S = zeros(1,NsnrS);

corr_corr_2DLimits_S = zeros(1,NsnrS);
corr_uncorr_2DLimits_S = zeros(1,NsnrS);

corr_corr_1DLimits_T = zeros(1,NsnrT);
corr_uncorr_1DLimits_T = zeros(1,NsnrT);

corr_corr_2DLimits_T = zeros(1,NsnrT);
corr_uncorr_2DLimits_T = zeros(1,NsnrT);
%% COMMENT OUT THIS SECTION IF YOU DO NOT HAVE THE AARHUSINV FORWARD
g = waitbar(0,['TEM Analysis',num2str(i),' out of ',num2str(NsnrT)]);
a = tic;

for i = 1:NsnrT
a2 = tic;

    SNR = SNRsT(i);

    if i == 1
        redo_priors = 1;
    else
        redo_priors = 0;
    end


    waitbar((4*(i-1)+1)/(4*NsnrT),g,['TEM: Noise: Corr -- Assumption: Corr -- SNR: ',num2str(SNR),' (',num2str(i),'/',num2str(NsnrT),')']);
    [corr_corr_T,corr_corr_O,ccTs,ccOs,AccT,AccO] = Primary(SNR,Ns(i),'IsCorr_AssumedCorr',1,rundate,redo_priors,1,'TEM',Nnoise);

    waitbar((4*(i-1)+2)/(4*NsnrT),g,['TEM: Noise: Corr -- Assumption: UnCorr -- SNR: ',num2str(SNR),' (',num2str(i),'/',num2str(NsnrT),')'])
    [corr_uncorr_T,corr_uncorr_O,cuTs,cuOs,AcuT,AcuO] = Primary(SNR,Ns(i),'IsCorr_AssumedUnCorr',1,rundate,redo_priors,0,'TEM',Nnoise);

    Acu1DsTEM{i} = AcuO;
    Acc1DsTEM{i} = AccO;
    Acu2DsTEM{i} = AcuT;
    Acc2DsTEM{i} = AccT;

    corr_corr_1DLimits_T(i) = corr_corr_O;
    corr_uncorr_1DLimits_T(i) = corr_uncorr_O;

    corr_corr_2DLimits_T(i) = corr_corr_T;
    corr_uncorr_2DLimits_T(i) = corr_uncorr_T;

    cc1DLsTEM{i} = ccOs;
    cu1DLsTEM{i} = cuOs;

    cc2DLsTEM{i} = ccTs;
    cu2DLsTEM{i} = cuTs;
end

delete(g)
%%
g = waitbar(0,['SEIS Analysis',num2str(i),' out of ',num2str(NsnrS)]);

for i = 1:NsnrS
    SNR = SNRsS(i);
    if i == 1
        redo_priors = 1;
    else
        redo_priors = 0;
    end
    waitbar((4*(i-1)+3)/(4*NsnrS),g,['Seis: Noise: Corr -- Assumption: Corr -- SNR: ',num2str(SNR),' (',num2str(i),'/',num2str(NsnrS),')']);
    [corr_corr_T,corr_corr_O,ccTs,ccOs,AccT,AccO] = Primary(SNR,Ns(i),'IsCorr_AssumedCorr',1,rundate,redo_priors,1,'seis',Nnoise);

    waitbar((4*(i-1)+4)/(4*NsnrS),g,['Seis: Noise: Corr -- Assumption: UnCorr -- SNR: ',num2str(SNR),' (',num2str(i),'/',num2str(NsnrS),')'])
    [corr_uncorr_T,corr_uncorr_O,cuTs,cuOs,AcuT,AcuO] = Primary(SNR,Ns(i),'IsCorr_AssumedUnCorr',1,rundate,redo_priors,0,'seis',Nnoise);
    
    Acu1Ds{i} = AcuO;
    Acc1Ds{i} = AccO;
    Acu2Ds{i} = AcuT;
    Acc2Ds{i} = AccT;

    corr_corr_1DLimits_S(i) = corr_corr_O;
    corr_uncorr_1DLimits_S(i) = corr_uncorr_O;

    corr_corr_2DLimits_S(i) = corr_corr_T;
    corr_uncorr_2DLimits_S(i) = corr_uncorr_T;

    cc1DLs{i} = ccOs;
    cu1DLs{i} = cuOs;

    cc2DLs{i} = ccTs;
    cu2DLs{i} = cuTs;

 SNRtimes(i) = toc(a2);

end

timer1 = toc(a);
delete(g)
disp(['Time Consumed:',num2str(timer1/3600),' hours'])

keyboard
%% PLOT RESULTS

plottype = 'TEM';
%plottype = 'seis'

switch plottype
    case 'TEM'
        SNRs = SNRsT;
        corr_corr_1DLimits = corr_corr_1DLimits_T;
        corr_corr_2DLimits = corr_corr_2DLimits_T;
        corr_uncorr_1DLimits = corr_uncorr_1DLimits_T;
        corr_uncorr_2DLimits = corr_uncorr_2DLimits_T;
        ymin = 1;
        ymax = 100;
        xl = 'Signal-To-Noise Ratio (SNR_{log})';
    case 'seis'
        SNRs = SNRsS;
        corr_corr_1DLimits = corr_corr_1DLimits_S;
        corr_corr_2DLimits = corr_corr_2DLimits_S;
        corr_uncorr_1DLimits = corr_uncorr_1DLimits_S;
        corr_uncorr_2DLimits = corr_uncorr_2DLimits_S;
        ymin = 0.5;
        ymax = 50;
        xl = 'Signal-To-Noise Ratio (SNR)';
end
fact = log10(SNRs(end)./SNRs(1));
xmin = SNRs(1)/fact;
xmax = SNRs(end)*fact;

figure()
hold on
set(gca,'yscale','log','xscale','log')
box on
plot(SNRs,corr_corr_1DLimits,'.b','MarkerSize',20,'displayname','Noise is correlated and assumed correlated, 1D Prior')
plot(SNRs,corr_corr_2DLimits,'.r','MarkerSize',20,'displayname','Noise is correlated and assumed correlated, 2D Prior')
plot(SNRs,corr_uncorr_1DLimits,'.g','MarkerSize',20,'displayname','Noise is correlated and assumed uncorrelated, 1D Prior')
plot(SNRs,corr_uncorr_2DLimits,'.k','MarkerSize',20,'displayname','Noise is correlated and assumed uncorrelated, 2D Prior')

grid on
xlim([xmin xmax])
ylim([ymin ymax])
hold off
legend
xlabel(xl,'FontSize',8)
ylabel('Resolution Limit [m]','FontSize',18)