function LookupTable = WriteLookupTable(N,intervals,profile_length,resolution,type)  

% Set parameters for wedge realizations %
switch type
    case 'TEM'
        wedge_depth = intervals(1,:); %Depth interval of top of wedge in meters, format: [minimum maximum]
        wedge_inclination = intervals(2,:); %wedge dip angle interval in degrees, format: [minimum maximum]
        x0 = intervals(3,:);
        length = intervals(4,:);
        incs = intervals(5,:); %slopes
        ress = intervals(6:8,:);
        dummyress = intervals(9:11,:);
    case 'seis'
        wedge_depth = intervals(1,:); %Depth interval of top of wedge in meters, format: [minimum maximum]
        wedge_inclination = intervals(2,:); %wedge dip angle interval in degrees, format: [minimum maximum]
        x0 = intervals(3,:);
        length = intervals(4,:);
        incs = intervals(5,:); %slopes
        rhos = intervals(6:8,:);
        vps = intervals(9:11,:);        
end
%Lookuptable will be N*profile_length/resolution+1 long
n = (profile_length/resolution+1);
Nrows = N*n;

%Lookuptable will have [RealID,traceID,ress,thick,theta,x0,length,inc]
Ncols = 11;

LookupTable = NaN(Nrows,Ncols);

count = 0;
disp('Producing Lookuptable')
A = tic();

switch type
    case 'TEM'
        for i = 1:N
            [pos, rezz, ~, thick, theta, xnull, lengths, inc] = get_TDEM_wedge(x0,wedge_depth,wedge_inclination,profile_length,resolution,ress,dummyress,length,incs);

            for j = 1:n
                count = count+1;

                LookupTable(count,1) = i;
                LookupTable(count,2) = j;
                LookupTable(count,3:5) = rezz;
                LookupTable(count,6:7) = thick(j,:);
                LookupTable(count,8) = theta;
                LookupTable(count,9) = xnull;
                LookupTable(count,10) = lengths;
                LookupTable(count,11) = inc;

                
            end
        end
    case 'seis'
        for i = 1:N
            [pos, Rhos, ~, thick, theta, xnull, lengths, inc] = get_seismic_wedge(x0,wedge_depth,wedge_inclination,profile_length,resolution,rhos,vps,length,incs);

            for j = 1:n
                count = count+1;

                LookupTable(count,1) = i;
                LookupTable(count,2) = j;
                LookupTable(count,3:5) = [Rhos(1),Rhos(2),Rhos(1)];
                LookupTable(count,6:7) = thick(j,:);
                LookupTable(count,8) = theta;
                LookupTable(count,9) = xnull;
                LookupTable(count,10) = lengths;
                LookupTable(count,11) = inc;
            end
        end
end

disp('Calculating Forward Responses')

switch type
    case 'TEM'
        T = LookupTable(:,6:7);
        R = LookupTable(:,3:5);
        
        Fmode = 100;

        [~,fwr] = get_forwards(count,1,1,1,1,1,Fmode,1,T,R);
        FWR = fwr(:,:,2);

        LookupTable = [LookupTable,FWR];

    case 'seis'


        sampf = 4;

        FWR = NaN(Nrows,26);
        count = 0;
        G = waitbar(0,['Calculating ',num2str(N), ' forward seismic responses.']);
        Z = tic();

        for ii = 1:N


            if mod(ii,5000) == 0


                Progress = ii;
                target = N;
                TimeSpent = toc(Z);
                TimePerProgress = TimeSpent/Progress;
                RemainingProgress = N-Progress;
                RemainingTime = RemainingProgress*TimePerProgress;

                PercProg = Progress/target;

                waitbar(PercProg,G,['Calculating ',num2str(N*n), ' forward seismic responses. (',num2str(RemainingTime/3600),'H left)']);

            end

            for jj = 1:n
                count = count+1;
                cur_rhos = LookupTable(count,3:5)';
                cur_thicknesses = LookupTable(count,6:7)';
                cur_vps = [2000 2000 2000]';

                [~,fwr] = seismic_forward(cur_thicknesses,cur_rhos,cur_vps,sampf);
                FWR(count,:) = fwr;
            end
        end

        LookupTable = [LookupTable,FWR];
end
B = toc(A);
disp(['Lookuptable completed. Time Elapsed: ',num2str(B),' seconds'])

LookupTableName = ['LookupTable',type,'2D_N',num2str(N),'.mat'];
save(LookupTableName,"LookupTable",'-v7.3')
end