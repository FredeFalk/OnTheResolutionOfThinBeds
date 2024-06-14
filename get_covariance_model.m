function [varmodel C_fin] = get_covariance_model(fwr,y,dy,NT)

%A gaussian semivariogram model is chosen due to its parabolic nature near the origin,
%reflecting the smooth seismic wavelet
gaus = @(x,a,S) S*exp(-3*x.^2/a.^2);

%perform some gridsearch for variables Sill (S) and Practical Range (a)
as = 9:0.05:13;
Ss = 0.003:0.0002:0.0085;

%as = 20.5:0.02:23;

likelihoodgrid = zeros(length(as),length(Ss));

datz = fwr;

data = fwr;
rem = [];
for i = 1:size(data,2)
    i = size(data,2)+1-i;
    test = std(data(:,i));
    if test == 0
        k = i;
        rem = [rem k];
    end
end

for i = rem
    datz(:,i) = [];
end


rem = [];
for i = 1:size(data,1)
    i = size(data,1)+1-i;
    test = std(data(i,:));
    if test == 0
        k = i;
        rem = [rem k];
    end
end

for i = rem
    datz(i,:) = [];
end

for i = 1:length(as)
    for j = 1:length(Ss)
        
        model = @(x) gaus(x,as(i),Ss(j));
        C = get_covariance_matrix(datz,model,4,0,1);
        
        D = exp(-0.5*sum(log(2*pi*eig((C))))-0.5*(datz(:)'/(C)*datz(:)));
        likelihoodgrid(i,j)=D;
    end
end

[A S] = find(max(max(likelihoodgrid))==likelihoodgrid);
aval = as(A);
sval = Ss(S);

display('effective range is:')
aval
display('Sill is:')
sval


%translate model to covariance model
varmodel = @(x) gaus(x,aval,sval);
C_fin = get_covariance_matrix(datz,varmodel,2,0,NT);

end