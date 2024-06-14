function [K2, marginals]=resolutionplot5(m,snr,mode,true_pars,rundate,it,type)
prem = m;

indsies=find(m(:,1)>0);
m=m(indsies,:);

R = 2000; %Grid resolution
%create grid

xs = linspace(-100,300,R);
ys = linspace(0,100,R);

%% Determine if reference wedge is within the 2.5th and 97.5th percentiles of the sampled wedges
ReferenceDepth = true_pars(1);
ReferenceTheta = true_pars(2);
ReferenceX0 = true_pars(3);
ReferenceLen = true_pars(4);
ReferenceInc = true_pars(5);

SampledDepths = m(:,7);
SampledThetas = m(:,8);
SampledX0s = m(:,10);
SampledLens = m(:,12);
SampledIncs = m(:,13);

DepthLowBound = prctile(SampledDepths,2.5);
DepthUpBound = prctile(SampledDepths,97.5);

SampledModels = numel(unique(SampledDepths));
SamplingCriteria = SampledModels >= 100;

DepthCriteria = and(DepthLowBound < ReferenceDepth,DepthUpBound > ReferenceDepth);

for i = 1:R
    curx = xs(i);
    curthicks = max(0,tand(SampledThetas).*(curx-SampledX0s));

    ThickLowBound = prctile(curthicks,2.5);
    ThickUpBound = prctile(curthicks,97.5);

    TLB(i) = ThickLowBound;
    TUB(i) = ThickUpBound;
    DUB(i) = DepthUpBound;
    DLB(i) = DepthLowBound;

    refthick = tand(ReferenceTheta).*(curx-ReferenceX0);
    RefThicks(i) = refthick;

    ThicknessCriteria(i) = and(ThickLowBound<=refthick,ThickUpBound>=refthick);

    ResolutionCriteria(i) = and(DepthCriteria,ThicknessCriteria(i));
end



tic
% CF=figure('Position',[100 100 700 500])
% for s = 1:4
%
P = zeros(R,R);
%     t = s+1;
%     if s < 4
%         m = prem(t:t,:);
%     else
%         m = prem;
%     end

for k = 1:size(m,1)
    for i = 1:R
        x = xs(i);
        A = m(k,13);

        th = tand(m(k,8))*(m(k,13)*(x-m(k,10)));
        th = (th+abs(th))/2;

        if A*x < A*(m(k,10)+A*m(k,12))

            for j = 1:R
                y = ys(j);
                if y>m(k,7)
                    if y<m(k,7)+th
                        P(i,j) = P(i,j)+1;
                    end
                end
            end
        end
    end
end

toc

P = P/size(m,1);

XS = (xs+100)*0.1;
MaMat = max(P');

K = find(and(MaMat>0.95,ResolutionCriteria));

if SamplingCriteria
    if isempty(K)
        K2 = NaN;
    else
        K2 = XS(K(1));
    end
else
    K2 = NaN;
end


%% Graphics for export
CF = figure('WindowState','maximized');
ax1 = subplot(10,1,[1,2,3,4]);
imagesc(XS,-ys,P')
hold on
plot([0 10 10 0],[-10 -10 -20 -10],'--r','LineWidth',2)
set(gca,'ydir','normal')
title(['SNR: ',num2str(snr),' __ ','Lim:',num2str(K2),' Unique Models: ',num2str(SampledModels)],'Interpreter','none')
ylim([-40,0])
xlim([0 10])
ylabel('Depth [m]')
yyaxis right
plot(XS,MaMat,'-k','LineWidth',2)
ylabel('Maximum Pixelwise Probability')
ylim([-0.1,1.1])
set(ax1,'ycolor','black')

ax2 = subplot(10,1,[5,6]);
plot(XS,ones(1,R),'-k','HandleVisibility','off')
plot(XS,zeros(1,R),'-k','HandleVisibility','off')
hold on
plot(XS,ResolutionCriteria.*MaMat>0.95,'-r','LineWidth',2,'DisplayName','Wedge Resolved')
plot(XS,ResolutionCriteria,'--r','DisplayName','Accuracy Criterion')
plot(XS,MaMat>0.95,'--b','DisplayName','Precision Criterion')
grid on
ylim([-0.1,1.1])
xlim([0 10])
legend('Location','SouthOutside','Orientation','horizontal')
ylabel('True / False')
set(gca,'FontSize',6)

ax3 = subplot(10,1,[7,8]);
plot(XS,DLB,'-r','HandleVisibility','off','LineWidth',2)
hold on
plot(XS,DLB*0+ReferenceDepth,'-b','DisplayName','Reference Depth','LineWidth',2)
plot(XS,DUB,'-r','DisplayName','Posterior 95% Confidence Interval','LineWidth',2)
legend('Location','SouthOutside','Orientation','horizontal')
grid on
ylim([-10,30])
xlim([0 10])
ylabel('Depth [m]')
set(gca,'FontSize',6)

ax4 = subplot(10,1,[9,10]);
plot(XS,TLB,'-r','HandleVisibility','off','LineWidth',2)
hold on
plot(XS,RefThicks,'-b','DisplayName','Reference Thickness','LineWidth',2)
plot(XS,TUB,'-r','DisplayName','Posterior 95% Confidence Interval','LineWidth',2)
grid on
ylim([-10 50])
xlim([0 10])
ylabel('Thickness [m]')
legend('Location','SouthOutside','Orientation','horizontal')
set(gca,'FontSize',6)
%%
localdir = pwd;
folder = [localdir,'\',type,'\Results\',rundate,'\',mode,'\2D_Prior\'];

if it > 0
    name = [folder,'2D_SNR_',num2str(snr),'_Noise_',num2str(it),'_lim_',num2str(K2),'.png'];
else
    name = [folder,'2D_SNR_',num2str(snr),'_lim_',num2str(K2),'.png'];
end

if ~isfolder(folder)
    mkdir(folder)
end

exportgraphics(CF,name,'Resolution',300)
close(CF)

display('The resolution limit in the 2D case is:')
display(num2str(K2))

end