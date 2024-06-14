function Prior1D = translate2Dto1D(Prior2D,Prior1D,type)
%% THIS CODE TRANSLATES 2D Prior to 1D Prior, while also producing a publication-figuere
Num2D = size(Prior2D,1);
Num1D = size(Prior1D,1);
NVERT = 1000;
MarginDepths = linspace(0,100,NVERT);
MarginalDistribution = MarginDepths*0;


%Number of 1D prior models for each 2D prior model
N1D2D = Num1D/Num2D;

P2D = zeros(NVERT,N1D2D);

Prior1D = reshape(Prior1D,[8,N1D2D,Num2D]);

%Determine the extent of the 2D prior
Minx = min(Prior2D(:,3)+Prior2D(:,5).*Prior2D(:,4));
Maxx = max(Prior2D(:,3)+Prior2D(:,5).*Prior2D(:,4));
xWidth = Maxx-Minx;

%Produce vector for xvalues for P2D and P1D
XS = linspace(Minx,Maxx,N1D2D);

%START TRANSLATING (Also calculate P2D, the pixelwise probability of the 2D Prior)
for i = 1:Num2D
    %Take one 2D realization
    Cur2D = Prior2D(i,:);
    
    %Where's the current 2D wedge located along x?
    CurWedgeExtent = [Cur2D(3), Cur2D(3)+Cur2D(5)*Cur2D(4)];
    CurMin = min(CurWedgeExtent);
    CurMax = max(CurWedgeExtent);

    %Take "N1D2D" random x-coordinate points
    points = XS;

    %Find the current wedge depth at all points
    
    PointsAboveWedge = and(points>=CurMin,points<=CurMax);
    
    Thicks = PointsAboveWedge.*(points-Cur2D(3))*Cur2D(5)*tand(Cur2D(2));
    Depths = Cur2D(1).*ones(1,N1D2D);

    Prior1D(1,:,i) = Depths;
    Prior1D(2,:,i) = Thicks;
    Prior1D(3:5,:,i) = repmat(Cur2D(6:8)',[1,N1D2D]);
    Prior1D(6:8,:,i) = repmat(Cur2D(9:11)',[1,N1D2D]);

    %Fill in P2D
    for ii = 1:N1D2D
        LogicalBelowTop = Depths(ii) < MarginDepths;
        LogicalAboveBottom = (Depths(ii)+Thicks(ii)) > MarginDepths;

        indices2D = and(LogicalBelowTop,LogicalAboveBottom);
        P2D(indices2D,ii) = P2D(indices2D,ii)+1;
    end
end

Prior1D=reshape(Prior1D,[8, Num1D]);
P2D = P2D/Num2D;

save([type,'Distribution1D.mat'],"Prior1D",'-v7.3')
%% PRODUCTION FIGURE
Fig3a = figure('Units','centimeters','Position',[0,0,18.3,17]);

% HISTROGRAMS FOR 2D PRIOR PARAMETERS
w = 7;
h = 9;

subplot(h,w,10)
hist(Prior2D(:,5),10)
xlim([-2,2])
ylim([0,Num2D])
title('Slope Direction')
set(gca,'YTick',[])

subplot(h,w,8)
hist(Prior2D(:,4),30)
xlim([min(Prior2D(:,4)),max(Prior2D(:,4))])
title('Length [m]')
set(gca,'YTick',[])

subplot(h,w,9)
hist(Prior2D(:,3),30)
xlim([min(Prior2D(:,3)),max(Prior2D(:,3))])
title('Tip X [m]')
set(gca,'YTick',[])

subplot(h,w,15)
hist(Prior2D(:,1),30)
xlim([min(Prior2D(:,1)),max(Prior2D(:,1))])
title('Depth [m]')
set(gca,'YTick',[])

subplot(h,w,16)
hist(Prior2D(:,2),30)
xlim([min(Prior2D(:,2)),max(Prior2D(:,2))])
title('Inclination [^o]')
set(gca,'YTick',[])

% HISTROGRAMS FOR 1D PRIOR PARAMETERS
subplot(16,5,[4,9,14]+10)
hist(Prior1D(1,:),30)
xlim([min(Prior1D(1,:)),max(Prior1D(1,:))])
title('Depth [m]')
set(gca,'YTick',[])

subplot(16,5,[5,10,15]+10)
hist(Prior1D(2,:),30)
ylim([0 8E5*10])
xlim([0,max(Prior1D(2,:))])
title('Thickness [m]')
set(gca,'YTick',[])
set(gca,'XTick',[0 20 40 60])
%Calculate Marginal Distribution for P(Wedge|depth) and 2D Prior Distribution
for i = 1:Num1D/100
    Cur1D = Prior1D(:,i);

    LogicalBelowTop = Cur1D(1) < MarginDepths;
    LogicalAboveBottom = (Cur1D(1)+Cur1D(2)) > MarginDepths;

    indices = and(LogicalBelowTop,LogicalAboveBottom);

    MarginalDistribution(indices) = MarginalDistribution(indices)+1;
end

MarginalDistribution = MarginalDistribution/(Num1D/100);

P1D = zeros(NVERT,NVERT);

for i = 1:NVERT
    P1D(:,i) = MarginalDistribution;
end

MaxP = round(100*max(P2D(:)))/100;

subplot(15,16,[97:101,113:117,129:133,145:149])
imagesc(XS,MarginDepths,P2D)
xlabel('X [m]')
ylabel('Depth [m]')
clim([0,MaxP])

%
subplot(15,16,[107:111,123:127,139:143,155:159])
imagesc(XS,MarginDepths,P1D)
xlabel('X [m]')
ylabel('Depth [m]')
clim([0,MaxP])
%
MarginalDistribution2D = sum(P2D,2)/N1D2D;
%
subplot(15,16,[102,118,134,150])
fill(MarginalDistribution2D,MarginDepths,[0,0,0.6])
set(gca,'ydir','reverse')
set(gca,'YTick',[],'XTick',[0,0.1])
yyaxis right
ylabel('Marginal PDF - F(Depth)')
set(gca,'YColor','k')
set(gca,'YTick',[])
%
subplot(15,16,[112,128,144,160])
fill(MarginalDistribution,MarginDepths,[0,0,0.6])
set(gca,'ydir','reverse')
set(gca,'YTick',[],'XTick',[0,0.1])
yyaxis right
ylabel('Marginal PDF - F(Depth)')
set(gca,'YColor','k')
set(gca,'YTick',[])

subplot(30,16,[177:181])
imagesc(0:256,0:1,(1:256))
set(gca,'XAxisLocation','Top')
set(gca,'xtick',[0,256],'ytick',[],'xticklabel',{'0',num2str(MaxP)})
title('Pixelwise P(Wedge)','VerticalAlignment','top')

subplot(30,16,[187:191])
imagesc(0:256,0:1,1:256)
set(gca,'XAxisLocation','Top')
set(gca,'xtick',[0,256],'ytick',[],'xticklabel',{'0',num2str(MaxP)})
title('Pixelwise P(Wedge)','VerticalAlignment','top')

%Posterior plots
subplot(15,16,[97:101,113:117,129:133,145:149]+80)
imagesc(XS,MarginDepths,P2D*0)
xlabel('X [m]')
ylabel('Depth [m]')
clim([0,MaxP])

subplot(15,16,[97:101,113:117,129:133,145:149]+90)
imagesc(XS,MarginDepths,P2D*0)
xlabel('X [m]')
ylabel('Depth [m]')
clim([0,MaxP])

subplot(15,16,[112,128,144,160]+70)
set(gca,'ydir','reverse')
set(gca,'YTick',[],'XTick',[0,0.1])
yyaxis right
ylabel('Marginal PDF - F(Depth)')
set(gca,'YColor','k')
set(gca,'YTick',[])
box on

subplot(15,16,[112,128,144,160]+80)
set(gca,'ydir','reverse')
set(gca,'YTick',[],'XTick',[0,0.1])
yyaxis right
ylabel('Marginal PDF - F(Depth)')
set(gca,'YColor','k')
set(gca,'YTick',[])
box on

%Draw Bounding Boxes
annotation(Fig3a,"rectangle",[0.52 0.025 0.46 0.875],'LineWidth',2)
annotation(Fig3a,"rectangle",[0.02 0.025 0.46 0.875],'LineWidth',2)

annotation(Fig3a,'textbox',[0.02 0.9 0.46 0.1],'string',"\rho_{2D}(\bf{m}\rm_{2D})",'FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none')
annotation(Fig3a,'textbox',[0.52 0.9 0.46 0.1],'string',"\rho_{1D}(\bf{m}\rm_{1D})",'FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none')

exportgraphics(Fig3a,'Fig3WithPosteriorVectorformat.eps','Resolution',300)

%% PRODUCTION FIGURE
Fig3b = figure('Units','centimeters','Position',[0,0,19,15]);

% HISTROGRAMS FOR 2D PRIOR PARAMETERS
ax=subplot(6,7,17);
hist(Prior2D(:,5),30)
xlim([-2,2])
ylim([0,Num2D])
title('s')
set(gca,'XTick',[-1,1],'YTick',[])
ax.YAxis.Exponent = 0;

ax=subplot(6,7,15);
hist(Prior2D(:,4),30)
xlim([min(Prior2D(:,4)),max(Prior2D(:,4))])
title('l [m]')
ax.YAxis.Exponent = 0;
set(gca,'YTick',[],'XTick',[250 350])

ax=subplot(6,7,16);
hist(Prior2D(:,3),30)
xlim([min(Prior2D(:,3)),max(Prior2D(:,3))])
title('t [m]')
ax.YAxis.Exponent = 0;
set(gca,'YTick',[],'XTick',[-100,100])

ax=subplot(6,7,8);
hist(Prior2D(:,1),30)
xlim([min(Prior2D(:,1)),max(Prior2D(:,1))])
title('z [m]')
ax.YAxis.Exponent = 0;
set(gca,'YTick',[],'XTick',[5 10 15])

ax=subplot(6,7,9);
hist(Prior2D(:,2),30)
xlim([min(Prior2D(:,2)),max(Prior2D(:,2))])
title('\theta [^o]','Interpreter','tex')
ax.YAxis.Exponent = 0;
set(gca,'YTick',[],'XTick',[4 6 8])

% HISTROGRAMS FOR 1D PRIOR PARAMETERS
ax=subplot(16,5,[9,14,19,24]+10);
hist(Prior1D(1,:),30)
xlim([min(Prior1D(1,:)),max(Prior1D(1,:))])
title('Z [m]')
ax.YAxis.Exponent = 0;
set(gca,'YTick',[],'XTick',[5 10 15])

ax=subplot(16,5,[10,15,20,25]+10);
hist(Prior1D(2,:),30)
ylim([0 8E5*10])
xlim([0,max(Prior1D(2,:))])
title('T [m]')
ax.YAxis.Exponent = 0;
ax.YAxisLocation = "right";
set(gca,'YTick',[],'XTick',[0,20,40,60])

%Calculate Marginal Distribution for P(Wedge|depth) and 2D Prior Distribution
for i = 1:Num1D/100
    Cur1D = Prior1D(:,i);

    LogicalBelowTop = Cur1D(1) < MarginDepths;
    LogicalAboveBottom = (Cur1D(1)+Cur1D(2)) > MarginDepths;

    indices = and(LogicalBelowTop,LogicalAboveBottom);

    MarginalDistribution(indices) = MarginalDistribution(indices)+1;
end

MarginalDistribution = MarginalDistribution/(Num1D/100);

P1D = zeros(NVERT,N1D2D);

for i = 1:N1D2D
    P1D(:,i) = MarginalDistribution;
end

MaxP = round(100*max(P2D(:)))/100;

subplot(10,16,[97:101,113:117,129:133,145:149])
imagesc(XS,MarginDepths,P2D)
xlabel('X [m]')
ylabel('Depth [m]')
clim([0,MaxP])

%
subplot(10,16,[107:111,123:127,139:143,155:159])
imagesc(XS,MarginDepths,P1D)
xlabel('X [m]')
ylabel('Depth [m]')
clim([0,MaxP])
%
MarginalDistribution2D = sum(P2D,2)/N1D2D;
%
subplot(10,16,[102,118,134,150])
fill(MarginalDistribution2D,MarginDepths,[0,0,0.6])
set(gca,'ydir','reverse')
set(gca,'YTick',[],'XTick',[])
xlim([0,0.2])
yyaxis right
ylabel('Marginal PDF - F(Depth)')
set(gca,'YColor','k')
set(gca,'YTick',[])
%
subplot(10,16,[112,128,144,160])
fill(MarginalDistribution,MarginDepths,[0,0,0.6])
set(gca,'ydir','reverse')
set(gca,'YTick',[],'XTick',[])
xlim([0,0.2])
yyaxis right
ylabel('Marginal PDF - F(Depth)')
set(gca,'YColor','k')
set(gca,'YTick',[])

subplot(20,16,[177:181])
imagesc(0:256,0:1,(1:256))
set(gca,'XAxisLocation','Top')
set(gca,'xtick',[0,256],'ytick',[],'xticklabel',{'0',num2str(MaxP)})
title('Pixelwise P(Wedge)','VerticalAlignment','top')

subplot(20,16,[187:191])
imagesc(0:256,0:1,1:256)
set(gca,'XAxisLocation','Top')
set(gca,'xtick',[0,256],'ytick',[],'xticklabel',{'0',num2str(MaxP)})
title('Pixelwise P(Wedge)','VerticalAlignment','top')

%Draw Bounding Boxes
annotation(Fig3b,"rectangle",[0.515 0.025 0.47 0.85],'LineWidth',2)
annotation(Fig3b,"rectangle",[0.015 0.025 0.47 0.85],'LineWidth',2)

annotation(Fig3b,"rectangle",[0.03 0.04 0.44 0.465],'LineWidth',1)
annotation(Fig3b,"rectangle",[0.03 0.51 0.44 0.35],'LineWidth',1)

annotation(Fig3b,"rectangle",[0.53 0.04 0.44 0.465],'LineWidth',1)
annotation(Fig3b,"rectangle",[0.53 0.51 0.44 0.35],'LineWidth',1)

annotation(Fig3b,'textbox',[0.02 0.9 0.46 0.1],'string',"\rho_{2D}(\bf{m}\rm_{2D})",'FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none')
annotation(Fig3b,'textbox',[0.52 0.9 0.46 0.1],'string',"\rho_{1D}(\bf{m}\rm_{1D})",'FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none')
annotation(Fig3b,'textbox',[0.52 0.77 0.46 0.1],'string',"\bf{m}\rm_{1D}",'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none')
annotation(Fig3b,'textbox',[0.02 0.77 0.46 0.1],'string',"\bf{m}\rm_{2D}",'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none')

%Axis labels
annotation(Fig3b,'textbox',[0.03 0.78 0.2 0.1],'string',"a)",'FontSize',16,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none')
annotation(Fig3b,'textbox',[0.53 0.78 0.2 0.1],'string',"c)",'FontSize',16,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none')
annotation(Fig3b,'textbox',[0.03 0.42 0.2 0.1],'string',"b)",'FontSize',16,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none')
annotation(Fig3b,'textbox',[0.53 0.42 0.2 0.1],'string',"d)",'FontSize',16,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none')

exportgraphics(Fig3b,'Fig3Vectorformat.eps','Resolution',300)
end