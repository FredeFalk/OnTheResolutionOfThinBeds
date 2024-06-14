function [noise,data,COV,cleand] = get_synthetic_TEM(x0,depth,theta,dx,resolution,res,length,inc,SNR,NT)

%Get exact model from prior generator
[pos, ress, ~, thick] = get_TDEM_wedge(x0,depth,theta,dx,resolution,res,res,length,inc);

N = size(pos,2);
data = [];


%Calculate forward response of model
Fmode = 3; %Ref model

[times,fwr] = get_forwards(1,1,1,1,1,1,Fmode,1,thick',ress);

for i = 1:N
    data = [data,fwr(i,:,2)];
end

cleand = data;

%Add noise
cov_import = load('COVAR.mat');
cov = cov_import.cov1D;

[data, noise, COV] = generate_noise_TEM(data,cov,SNR);

end