function [TD,OD,TDin,ODin,APPROXIMATE,APPROXIMATE1D] = Primary(SNR,Nreals,noisemode,switchnoise,rundate,redo_priors,noisemode0,type,NN)

N = 1; %number of variations
starts = 1;
ends = 100;

% Whats the type?
switch type
    case 'seis'
        disp('Settíng up synthetic seismic inverse problem')
        disp('')
    case 'TEM'
        disp('Settíng up synthetic TEM inverse problem')
        disp('')
end

% CODE PARAMETERS FOR REJECTION SAMPLING OF SEISMIC
true_depth = 10;
true_dip = atand(0.1);
true_rhos = [1500 1800 1500];
true_vps = [2000 2000 2000];
true_res = [70 10 70];
true_x0 = -100;
true_length = 400;
profile_length = 200;            %meters
resolution = 100;                 %Distance between traces, meters
true_inc = 1;                   %true slope direction, 1 for positive, -1 for negative.

true_params = [true_depth true_dip true_x0 true_length, true_inc];

% Set parameters for wedge realizations %
wedge_depth = [1 20]; %Depth interval of top of wedge in meters, format: [minimum maximum]
wedge_inclination = [3 9]; %wedge dip angle interval in degrees, format: [minimum maximum]
x0 = [-150 150];
length = [200 400];
incs = [1 -1]; %slopes

%N = profile_length/resolution+1;
rho1 = [1400 1600];
rho2 = [1800 1800];
rho3 = [1400 1600];
rho = [rho1;rho2;rho3];

vp1 = [2000 2000];
vp2 = [2000 2000];
vp3 = [2000 2000];
vp = [vp1;vp2;vp3];

res1 = [70,70];
res2 = [5,14];
res3 = [70,70];
ress = [res1;res2;res3];

switch type
    case 'TEM'
        intervals = [wedge_depth;wedge_inclination;x0;length;incs;ress;ress];
    case 'seis'
        intervals = [wedge_depth;wedge_inclination;x0;length;incs;rho;vp];
end
%% SETUP FOLDER STRUCTURE

switch type
    case 'seis'
        typefolder = '\seis';
        true_mod = [true_depth true_dip true_rhos true_vps true_x0 true_length profile_length resolution true_inc];

    case 'TEM'
        typefolder = '\TEM';
        true_mod = [true_depth true_dip true_res true_res true_x0 true_length profile_length resolution true_inc];

end

%Check if folders exist
Name1D = [pwd,typefolder,'\Results\',rundate,'\',noisemode,'\1D_Prior\'];
Check1D = isfolder(Name1D);

Name2D = [pwd,typefolder,'\Results\',rundate,'\',noisemode,'\2D_Prior\'];
Check2D = isfolder(Name2D);

if ~Check1D
    mkdir(Name1D)
end

if ~Check2D
    mkdir(Name2D)
end

%%
Vals = linspace(starts,ends,N);
MA = [];
Means = [];

NoiseType = 1; %1 for correlated noise, 0 for independent noise

%Produce data for 2D
[covm,noise,~,cleandata]=GetSyntheticData(SNR,true_mod,NoiseType,type);
cleandata = cleandata(:);

cov_original = covm;

if switchnoise == 0
    new_noise = chol(diag(diag(covm)))'*randn(78,1);
    new_noise = reshape(new_noise,[26,3]);

    noise = new_noise;

    CA = chol(diag(diag(covm)));
    CA = CA';
else
    CA = chol(covm);
    CA = CA';
end

if noisemode0 == 1
    covinv = inv(covm);
else
    covinv = inv(diag(diag(covm)));
end

shap = size(noise);
noise = noise(:);

noises = [];

for i = 1:Nreals
    noises = [noises CA*randn(prod(shap),1)];
end

%% PRODUCE 2D LOOKUP TABLE
UseExistingTEM = 1;
NLookUpTEM = 1E4;

UseExistingSEIS = 1;
NLookUpSEIS = 1E4;

switch type
    case 'TEM'
        switch UseExistingTEM
            case 1

                name = 'LookupTableTEM2D_N10000.mat';
                disp('Loading 2D TEM Lookup Table')
                %This choice loads a lookuptable
                %The lookuptable contains some number of forward responses
                %It has one column for each defining parameter
                %It has one column for each of the 33 gates
                TS = load(name);
                T = TS.LookupTable;
                T = single(T);
                clear TS

                disp('2D TEM Lookup Table Loaded')
            case 0
                %This choice produces a lookuptable
                T = WriteLookupTable(NLookUpTEM,intervals,profile_length,resolution,type);
                T = single(T);

        end
    case 'seis'
        switch UseExistingSEIS
            case 1

                name = 'LookupTableseis2D_N10000.mat';
                disp('Loading 2D SEIS Lookup Table')
                %This choice loads a lookuptable
                %The lookuptable contains some number of forward responses
                %It has one column for each defining parameter
                %It has one column for each of the 33 gates
                TS = load(name);
                T = TS.LookupTable;
                T = single(T);
                clear TS

                disp('2D SEIS Lookup Table Loaded')


            case 0
                %This choice produces a lookuptable
                T = WriteLookupTable(NLookUpSEIS,intervals,profile_length,resolution,type);
                T = single(T);

        end
end

%% SECTION FOR PRODUCING THE 1D PRIOR, THIS IS DERIVED FROM "Num2D" REALIZATIONS OF THE 2D PRIOR
%The 1D prior lookuptable is saved as a .mat file called "mar" which is subsequently
%used for sampling.

if redo_priors == 1
    Num2D = 2000;
    Num1Dper2D = 5;
    Num1D = Num2D*Num1Dper2D;

    Prior1D = zeros(Num1D,8);
    Prior2D = zeros(Num2D,11);

    %Get 2D realizations
    switch type
        case 'seis'
            for i = 1:Num2D
                [~, rhos, vps, thick2, theta2, x0b, len2, inc2] = get_seismic_wedge(x0,wedge_depth,wedge_inclination,profile_length,resolution,rho,vp,length,incs);
                Prior2D(i,:) = [thick2(1), theta2, x0b, len2, inc2, rhos', vps'];
            end
        case 'TEM'
            for i = 1:Num2D
                [~, res, resdummy, thick2, theta2, x0b, len2, inc2] = get_TDEM_wedge(x0,wedge_depth,wedge_inclination,profile_length,resolution,ress,ress,length,incs);
                Prior2D(i,:) = [thick2(1), theta2, x0b, len2, inc2, res', resdummy'];
            end
    end

    %Translate each 2D realization to "Num1Dper2D" realizations
    ODT=translate2Dto1D(Prior2D,Prior1D,type);

    clear Prior1D
    clear Prior2D
end



%% Sample with the 2D prior (wedge)
MA = zeros(Nreals,15);
Nnoise = NN;
MA1 = cell(1,Nnoise);

ACC_RATE = zeros(1,Nnoise);
APPROXIMATE = cell(1,Nnoise);
NumWorkers2D = 8;

parfor (i = 1:Nnoise, NumWorkers2D)
%for i = 1:Nnoise
    switch type
        case 'TEM'
            [Mt,PROP,isapprox] = MainSampler(Nreals/Nnoise,true_mod,covinv,cleandata+noises(:,i),cleandata,T,NLookUpTEM);
        case 'seis'
            [Mt,PROP,isapprox] = MainSampler(Nreals/Nnoise,true_mod,covinv,cleandata+noises(:,i),cleandata,T,NLookUpSEIS);
    end
    APPROXIMATE{i} = isapprox{1};

    MA1{i} = Mt;

    ACC_RATE(i) = Nreals/PROP;

    display(['Currently at..',num2str(i*(Nreals/Nnoise)),'...samples out of...',num2str(Nreals)])
end
%%

for i = 1:Nnoise
    ind1 = (i-1)*(Nreals/Nnoise)+1;
    ind2 = (Nreals/Nnoise)*i;

    arr = MA1{i};
    if numel(arr)
        MA(ind1:ind2,:) = arr;
    end
end
%% PLOT 2D DATA
D2fig1 = figure();
plot(cleandata(:),'-r','DisplayName','Without Noise')
hold on
grid on
for i = 1:Nnoise
    if i == 1
        plot(cleandata(:)+noises(:,i),'.-k','DisplayName','With Noise')
    else
        plot(cleandata(:)+noises(:,i),'.-k','HandleVisibility','off')
    end
end
s = size(cleandata,1);
V = 1:s;

errorbar(V,cleandata,2*sqrt(diag(covm)),'.r','displayname','2\sigma')
title(['2D Data and noise, SNR:', num2str(SNR)])
legend()
switch type
    case 'TEM'
        ylim([-25,0])
    case 'seis'
        ylim([-0.6,0.6])
end
exportgraphics(D2fig1,[Name2D,'Data2D_SNR_',num2str(SNR),'.png'])

%% PLOT P VALUES
Pfig1 = figure();

for i = 1:Nnoise
    nexttile
    curma = MA1{i};
    curP = curma(:,9);
    hist(curP)
    title(['Noise Real ',num2str(i)])
    grid on
    xlabel('P')
end

nexttile
curP = MA(:,9);
hist(curP)
title('All Noise Reals')
grid on
xlabel('P')

exportgraphics(Pfig1,[Name2D,'SNR_',num2str(SNR),'_Phist.png'])

Pfig2 = figure();

for i = 1:Nnoise
    nexttile
    curma = MA1{i};
    curP = curma(:,9);
    plot(sort(curP))
    title(['Noise Real ',num2str(i)])
    grid on
    ylabel('P')
end

nexttile
curP = MA(:,9);
plot(sort(curP))
title('All Noise Reals')
grid on
ylabel('P')

exportgraphics(Pfig2,[Name2D,'SNR_',num2str(SNR),'_Psort.png'])

%%

NUM = 100;
NSeeds = Nnoise;
wedgethicks = linspace(0,10,NUM);
cov1D = cov_original(1:(size(cov_original,1)/3),1:(size(cov_original,1)/3));


%% 1D: Assume uncorrelated or correlated noise?
if switchnoise == 0
    if noisemode0 == 1
        covinv1D = inv(cov1D);
        cov1D = diag(diag(cov1D));

    else
        covinv1D = inv(diag(diag(cov1D)));
        cov1D = diag(diag(cov1D));

    end
else
    if noisemode0 == 1
        covinv1D = inv(cov1D);

    else
        covinv1D = inv(diag(diag(cov1D)));
    end
end
%% Produce 1D reference data
close all

switch type
    case 'TEM'
        data1D = zeros(33,NUM);

        %Produce data for 1D
        thick = [true_depth*ones(NUM,1), wedgethicks'];

        Fmode = 3;
        [tim,fwr] = get_forwards(1,1,1,1,1,1,Fmode,1,thick',true_res);
        da1D = fwr(:,:,2);

        data1D = da1D';

    case 'seis'
        for i = 1:NUM
            %Produce data for 1D
            thick = [true_depth wedgethicks(i)];

            [tim, da1D] = seismic_forward(thick',true_rhos',true_vps',4);

            data1D(:,i) = da1D;
        end
end

times = tim;

%Prepare noise for 1D inversion
C=chol(cov1D);
C = C';
noise1Da = C*randn(size(noise,1)/3,NSeeds);

%% PRINT CERTAIN FIGURES
assumed_model = inv(covinv1D);
true_model = cov1D;

noisemod_fig=figure();
imagesc(assumed_model)
clim([0,max(assumed_model(:))])
cbar=colorbar;
cbar.Ticks = linspace(0,max(assumed_model(:)),5);
title('Assumed Covariance Model')
exportgraphics(noisemod_fig,[Name1D,num2str(SNR),'_noisemodel.png'])

truemod_fig=figure();
imagesc(true_model)
clim([0,max(true_model(:))])
cbar=colorbar;
cbar.Ticks = linspace(0,max(assumed_model(:)),5);
title('True Covariance Model')
exportgraphics(truemod_fig,[Name1D,num2str(SNR),'_noisetrue.png'])
APPROXIMATE1D = cell(size(APPROXIMATE));

%% PLOT 1D DATA
D2fig1 = figure();

plot(data1D(:,70),'-r','DisplayName','Without Noise')
hold on
grid on
for i = 1:Nnoise
    if i == 1
        plot(data1D(:,70)+noise1Da(:,i),'.-k','DisplayName','With Noise')
    else
        plot(data1D(:,70)+noise1Da(:,i),'.-k','HandleVisibility','off')
    end
end
V2 = size(data1D,1);

errorbar(1:V2,data1D(:,70),2*sqrt(diag(cov1D)),'.r','displayname','2\sigma')
title(['1D Data and noise, SNR:', num2str(SNR)])
legend()
switch type
    case 'TEM'
        ylim([-25,0])
    case 'seis'
        ylim([-0.6,0.6])
end
exportgraphics(D2fig1,[Name1D,'Data1D_SNR_',num2str(SNR),'.png'])
%% PRODUCE 1D LOOKUP TABLE
UseExistingTEM1D = 1;
UseExistingSEIS1D = 1;
A = tic();

N_use = 1E4;
switch type
    case 'TEM'
        switch UseExistingTEM1D
            case 1

                name = 'LookupTableTEM1D_N10000.mat';
                disp('Loading 1D TEM Lookup Table')
                %This choice loads a lookuptable
                %The lookuptable contains some number of forward responses
                %It has one column for each defining parameter
                %It has one column for each of the 33 gates
                TS_1D = load(name);
                LookupTable1D = TS_1D.LookupTable1D;
                clear TS_1D

                Ntotal = size(LookupTable1D,1);
                Randvector = randperm(Ntotal,N_use);
                
                LookupTable1D = LookupTable1D(Randvector,:);
                LookupTable1D = single(LookupTable1D);
                clear randvector

                disp('1D TEM Lookup Table Loaded')
            case 0
                A = tic();
                %This choice produces a 1D lookuptable from ODT
                disp('Calculating Forward Responses')

                Th = ODT(1:2,:)';
                Re = ODT(3:5,:)';
                N1D = size(ODT,2);

                Fmode = 100;
                [~,fwr] = get_forwards(N1D,1,1,1,1,1,Fmode,1,Th,Re);
                FWR = fwr(:,:,2);

                LookupTable1D = [ODT',FWR];
                LookupTable1D = single(LookupTable1D);
                
                LookupTableName = ['LookupTable',type,'1D_N',num2str(size(LookupTable1D,1)),'.mat'];
                save(LookupTableName,"LookupTable1D",'-v7.3')
                clear ODT
        end

    case 'seis'
        switch UseExistingSEIS1D
            case 1

                name = 'LookupTableseis1D_N10000.mat';
                disp('Loading 1D SEIS Lookup Table')
                %This choice loads a lookuptable
                %The lookuptable contains some number of forward responses
                %It has one column for each defining parameter
                %It has one column for each of the 33 gates
                TS_1D = load(name);

                LookupTable1D = TS_1D.LookupTable1D;
                clear TS_1D

                Ntotal = size(LookupTable1D,1);
                Randvector = randperm(Ntotal,N_use);

                LookupTable1D = LookupTable1D(Randvector,:);
                LookupTable1D = single(LookupTable1D);
                clear Randvector

                disp('1D SEIS Lookup Table Loaded')
            case 0
                
                %This choice produces a 1D lookuptable from ODT
                disp('Calculating SEIS Forward Responses')

                N1D = size(ODT,2);
                Z = tic();
                FWR = NaN(N1D,26);

                G = waitbar(0,['Calculating ',num2str(N1D), ' forward seismic responses.']);

                for ii = 1:N1D


                    if mod(ii,5000) == 0


                        Progress = ii;
                        target = N1D;
                        TimeSpent = toc(Z);
                        TimePerProgress = TimeSpent/Progress;
                        RemainingProgress = N1D-Progress;
                        RemainingTime = RemainingProgress*TimePerProgress;

                        PercProg = Progress/target;

                        waitbar(PercProg,G,['Calculating ',num2str(N1D), ' forward seismic responses. (',num2str(RemainingTime/3600),'H left)']);

                    end

                    
                    cur_rhos = ODT(3:5,ii);
                    cur_thicknesses = ODT(1:2,ii);
                    cur_vps = [2000 2000 2000]'; %VP is hardcoded to be 2E3 m/s

                    sampf = 4;

                    [~,fwr] = seismic_forward(cur_thicknesses,cur_rhos,cur_vps,sampf);
                    FWR(ii,:) = fwr;
        
                end

                delete(G)

                LookupTable1D = [ODT',FWR];
                LookupTable1D = single(LookupTable1D);
                LookupTableName = ['LookupTable',type,'1D_N',num2str(size(LookupTable1D,1)),'.mat'];
                save(LookupTableName,"LookupTable1D",'-v7.3')
                clear ODT
        end


end

B = toc(A);

disp(['Lookuptable completed. Time Elapsed: ',num2str(B),' seconds'])
clear T T2

%% INVERT 1D DATA

parfor i = 1:Nnoise

    curdata = data1D+repmat(noise1Da(:,i),[1,NUM]);
    
    switch type
        case 'TEM'
            [M,~,isapprox1D] = MainSampler1D(Nreals/Nnoise,covinv1D,curdata,LookupTable1D,N_use);
        case 'seis'
            [M,~,isapprox1D] = MainSampler1D(Nreals/Nnoise,covinv1D,curdata,LookupTable1D,N_use);
    end


    APPROXIMATE1D{i} = isapprox1D;

    MA0{i} = M;

    display(['Currently at..',num2str(i*(Nreals/Nnoise)),'...samples out of...',num2str(Nreals)])
end
close all

MA02 = [];
for i = 1:Nnoise
    MA02 = [MA02;MA0{i}];
end

% switch type
%     case 'seis'
%         if SNR == 5
%             keyboard
%         end
% 
%     case 'TEM'
%         if SNR == 500
%             keyboard
%         end
% end
%% CALC RES LIMITS
% ISSUE BEFORE THIS
PlotIndividualResLimits = 0;

if PlotIndividualResLimits == 1
    %Resolution Limits for each Noise Realization Separately
    for i = 1:Nnoise
        %Find the reals for noise model "i" in 1D results
        ODin(i)=resolutionplot4(MA0{i},NUM,SNR,noisemode,true_depth,rundate,i,type);
        TDin(i)=resolutionplot5(MA1{i},SNR,noisemode,true_params,rundate,i,type);
        close all
    end
else
    ODin(i) = NaN;
    TDin(i) = NaN;
end
%Combined Resolution Limits
OD=resolutionplot4(MA02,NUM,SNR,noisemode,true_depth,rundate,0,type);
TD=resolutionplot5(MA,SNR,noisemode,true_params,rundate,0,type);
close all
end