function [noise,data,cov,cleand] = get_synthetic_seismic(x0,depth,theta,dx,resolution,rho,vp,length,inc,SNR,NT)

%Get exact model from prior generator
[pos rhos vps thick] = get_seismic_wedge(x0,depth,theta,dx,resolution,rho,vp,length,inc);

N = size(pos,2);
data = [];

%Calculate forward response of model
for i = 1:N
    [times,fwr] = seismic_forward(thick(i,:)',rhos,vps,0.4);
    data = [data fwr];
end

%Do covariance analysis
[y,x,dy,dx]=get_semivariogram(data,resolution,0.4);

%Calculate coarser forward response of model

data = [];
for i = 1:N
    [times,fwr] = seismic_forward(thick(i,:)',rhos,vps,4);
    data = [data fwr];
end

[covarmodel,C] = get_covariance_model(data,y,dy,NT);

cleand = data;

%Add noise
cov = get_covariance_matrix(data,covarmodel,4,1,NT);
cov = cov/SNR;

[data, noise] = generate_noise(data,cov);
end