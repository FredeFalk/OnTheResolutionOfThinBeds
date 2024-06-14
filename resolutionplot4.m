function K2=resolutionplot4(Ms,NUM,snr,mode,ref,rundate,it,type)

M = Ms;

%% Ms contains realizations of the posterior with 3 numbers (3 rows)
% The first value is an integer identifying which of the reference models
% the sample is accepted for.
% The second value contains the sampled depth to wedge
% The third value contains the sampled wedge thickness

k = size(M,1);

CF = figure('Position',[100 100 700 500]);

NREAL = size(M,1)/100;
Ms2 = zeros(3,100,NREAL);

for i = 1:100

    indcs = find(M(:,1)==i);
    Nmatch = numel(indcs);

    Ms2(1,i,1:Nmatch) = M(indcs,1);
    Ms2(2,i,1:Nmatch) = M(indcs,2);
    Ms2(3,i,1:Nmatch) = M(indcs,3);

    %% Determine if reference wedge depth is within the 2.5th and 97.5th percentiles of the sampled depths
    ReferenceDepth = ref;
    ReferenceThickness = 0.1*i;
    RefThicks(i) = ReferenceThickness;
    
    NumUnique = numel(unique(squeeze(Ms2(2,i,1:Nmatch))));
    Uniques(i) = NumUnique;

    DepthLowBound = prctile(Ms2(2,i,1:Nmatch),2.5);
    DepthUpBound = prctile(Ms2(2,i,1:Nmatch),97.5);

    DepthCriteria = and(DepthLowBound < ReferenceDepth,DepthUpBound > ReferenceDepth);

    ThickLowBound = prctile(Ms2(3,i,1:Nmatch),2.5);
    ThickUpBound = prctile(Ms2(3,i,1:Nmatch),97.5);
    
    DUB(i) = DepthUpBound;
    DLB(i) = DepthLowBound;
    TUB(i) = ThickUpBound;
    TLB(i) = ThickLowBound;

    ThicknessCriteria = and(ThickLowBound < ReferenceThickness,ThickUpBound > ReferenceThickness);

    SamplingCriteria(i) = NumUnique >= 100;
    ResolutionCriteria(i) = and(DepthCriteria,ThicknessCriteria);
    %%

end

%%%% Make Matrix %%%%%
%Number of rows per thickness of wedge
MsCur = M;

NR = 2000;
MaxD = 100;
INC= MaxD/NR;
Mat = zeros(NR,100);
SumV = zeros(1,100);

k = size(M,1);

for i = 1:k

    colno = MsCur(i,1); %Read column number from the sample
    SumV(colno) = SumV(colno)+1; %Add count for normalization later

    for j = 1:NR
        D = j*INC;
        if D > MsCur(i,2)
            if D < MsCur(i,2)+MsCur(i,3)
                Mat(j,colno) = Mat(j,colno)+1;
            end
        end
    end
end

for i = 1:max(MsCur(:,1))
    Mat(:,i) = Mat(:,i)./SumV(i);
end

MaMat = max(Mat);
xs = linspace(0,10,NUM);
xs = xs(1:length(MaMat));

CombinedCriteria = and(MaMat(1:(length(MaMat)))>0.95,ResolutionCriteria);
FinalCriteria = and(CombinedCriteria,SamplingCriteria);

K = find(CombinedCriteria);
Ks = find(FinalCriteria);

if size(Ks,2)==0
    K2 = NaN;
else
    K2 = xs(Ks(1));
end

%% Graphics for export
CF = figure('WindowState','maximized');
ax1 = subplot(10,1,[1,2,3,4]);
imagesc(0.1:0.1:10,0:-0.1:(-100),Mat)
hold on
plot([0 10 10 0],[-10 -10 -20 -10],'--r','LineWidth',2)
set(gca,'ydir','normal')
title(['SNR: ',num2str(snr),' __ ','Lim:',num2str(K2),' Last Sounding ',num2str(Uniques(end)),' Unique Models'],'Interpreter','none')
ylim([-40,0])
ylabel('Depth [m]')
yyaxis right
plot(0.1:0.1:10,MaMat,'-k','LineWidth',2)
ylabel('Maximum Pixelwise Probability')
ylim([-0.1,1.1])
set(ax1,'ycolor','black')
set(gca,'FontSize',6)

ax2 = subplot(10,1,[5,6]);
plot(0.1:0.1:10,ones(1,100),'-k','HandleVisibility','off')
plot(0.1:0.1:10,zeros(1,100),'-k','HandleVisibility','off')
hold on
plot(0.1:0.1:10,FinalCriteria,'-r','LineWidth',2,'DisplayName','Wedge Resolved & Sampled')
plot(0.1:0.1:10,CombinedCriteria,'--r','DisplayName','Accuracy Criterion')
plot(0.1:0.1:10,ResolutionCriteria,'--r','DisplayName','Wedge Resolved')
plot(0.1:0.1:10,MaMat>0.95,'--b','DisplayName','Precision Criterion')
grid on
ylim([-0.1,1.1])
legend('Location','SouthOutside','Orientation','horizontal')
ylabel('True / False')
yyaxis right
plot(0.1:0.1:10,Uniques,'-c','DisplayName','Number of Unique Models')

ax3 = subplot(10,1,[7,8]);
plot(0.1:0.1:10,DLB,'-r','HandleVisibility','off','LineWidth',2)
hold on
plot(0.1:0.1:10,DLB*0+ref,'-b','DisplayName','Reference Depth','LineWidth',2)
plot(0.1:0.1:10,DUB,'-r','DisplayName','Posterior 95% Confidence Interval','LineWidth',2)
legend('Location','SouthOutside','Orientation','horizontal')
grid on
ylim([-10,30])
ylabel('Depth [m]')
set(gca,'FontSize',6)

ax4 = subplot(10,1,[9,10]);
plot(0.1:0.1:10,TLB,'-r','HandleVisibility','off','LineWidth',2)
hold on
plot(0.1:0.1:10,RefThicks,'-b','DisplayName','Reference Thickness','LineWidth',2)
plot(0.1:0.1:10,TUB,'-r','DisplayName','Posterior 95% Confidence Interval','LineWidth',2)
grid on
ylim([-10 50])
ylabel('Thickness [m]')
legend('Location','SouthOutside','Orientation','horizontal')
set(gca,'FontSize',6)
%%

localdir = pwd;
folder = [localdir,'\',type,'\Results\',rundate,'\',mode,'\1D_Prior\'];

if it > 0
    name = [folder,'1D_SNR_',num2str(snr),'_Noise_',num2str(it),'_lim_',num2str(K2),'.png'];
else
    name = [folder,'1D_SNR_',num2str(snr),'_lim_',num2str(K2),'.png'];
end

if ~isfolder(folder)
    mkdir(folder)
end

exportgraphics(CF,name,'Resolution',300)
close(CF)

disp('The resolution limit in the 1D case is:')
disp(num2str(K2))

end