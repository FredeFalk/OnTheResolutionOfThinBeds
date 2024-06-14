function [M,PROP,isapprox]=MainSampler1D(Nreals,covinv,data,T,NLU)
%Number of soundings per prior realization:
N = 1;

NLookUp = NLU;
NSounding = size(data,2);
m_acc = [];

isapprox = zeros(Nreals,NSounding);
for j = 1:NSounding
    M = zeros(1,NLookUp);
    curdata = data(:,j);

    for i = 1:NLookUp
        First = (i-1)*N+1;
        Last = i*N;

        curfwr = T(First:Last,9:end);
        MF = DataMisfit(curfwr,curdata,covinv);
        M(i) = MF;
    end

fmaxM = max(M);

% SAMPLE POSTERIOR DISTRIBUTION WITH 2D PRIOR
disp(['Sampling Posterior ',num2str(j),' with 1D prior'])
NREAL = Nreals;
m_acc_temp = NaN(NREAL,4);

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

m_acc_temp(:,1) = j;
m_acc_temp(:,2) = T(Original_Inds_Accepted,1);
m_acc_temp(:,3) = T(Original_Inds_Accepted,2);
m_acc_temp(:,4) = NREAL;

m_acc = [m_acc;m_acc_temp];
isapprox(:,j) = Original_Inds_Accepted;
end


PROP = NLookUp;
M = m_acc;
end