function [M,PROP,isapprox]=MainSampler(Nreals,trues,covinv,data,cleandata,T,NLU)
profile_length = trues(11);            %meters
resolution = trues(12);                 %Distance between traces, meters

cleandata_fmax = DataMisfit(cleandata,data,covinv);

N = profile_length/resolution+1;

NLookUp = NLU;

M = zeros(1,NLookUp);

for i = 1:NLookUp
    First = (i-1)*N+1;
    Last = i*N;

    curfwr = T(First:Last,12:end);
    curfwr = curfwr';
    cfwr = curfwr(:);

    MF = DataMisfit(cfwr,data,covinv);

    if isnan(MF)
        MF = -1E9;
    end

    M(i) = MF;
end

fmaxM = max(M);
fmax = max(fmaxM,cleandata_fmax);

% SAMPLE POSTERIOR DISTRIBUTION WITH 2D PRIOR
display('Sampling Posterior with 2D prior')
NREAL = Nreals;
m_acc = zeros(NREAL,15);

%Get Acc_P_values
Acc_P_values = exp(M-fmaxM);
[SortedP,iS] = sort(Acc_P_values);
CumulativeP = cumsum(SortedP)/sum(SortedP);

%Draw NREAL random numbers
real_rands = rand(1,NREAL);

%Find out where the cumulative P and real_rands intersect
AcceptanceMatrix = CumulativeP>real_rands';

%Get accepted models
SI = sum(AcceptanceMatrix==0,2)+1;

%Find accepted models original index
Original_Inds_Accepted = iS(SI);

m_acc(:,14:15) = [real_rands',M(Original_Inds_Accepted)'];

Original_Inds_Accepted = 3*Original_Inds_Accepted;

%format: [res1,res2,res1,res1,res2,res1,thick1,theta,x0,lengths,inc]
m_acc(:,1:13) =[T(Original_Inds_Accepted,3:5) T(Original_Inds_Accepted,3:6) T(Original_Inds_Accepted,8) Acc_P_values(Original_Inds_Accepted/3)' T(Original_Inds_Accepted,9) M(Original_Inds_Accepted/3)' T(Original_Inds_Accepted,10:11)];

isapprox = {Original_Inds_Accepted};
PROP = NLookUp;
M = m_acc;
end