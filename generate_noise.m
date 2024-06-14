function [out, noise2] = generate_noise(fwr,cov)
%Function for adding noise to forward response, fwr.
%std is the standard deviation of the uncorrelated noise
%cov is the covariance matrix for the correlated noise

noise1 = randn(size(fwr));

L = chol((cov));
L = L';

sized = size(noise1);
noise2 = L*noise1(:);

out = reshape(fwr(:)+noise2,sized);
end