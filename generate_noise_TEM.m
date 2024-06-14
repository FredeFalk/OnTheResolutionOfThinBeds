function [out, noise2, COV] = generate_noise_TEM(fwr,cov,SNR)
%Function for adding noise to forward response, fwr.
%std is the standard deviation of the uncorrelated noise
%cov is the covariance matrix for the correlated noise

noise1 = randn(size(fwr));

N = numel(fwr)/size(cov,1);

CovCell = repmat({cov}, 1, N);
BigCov = blkdiag(CovCell{:});

% ENSURE THAT THE MEAN SNR IS 1 %
ScalingFactor= 1/mean(power(fwr,2)./diag(BigCov)');
BigCov = BigCov./ScalingFactor;

BigCov = BigCov/SNR;

L = chol((BigCov));
L = L';

sized = size(noise1);
noise2 = L*noise1(:);

COV = BigCov;

out = reshape(fwr(:)+noise2,sized);
end