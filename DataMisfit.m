function loglike=DataMisfit(d,d_obs,covinv)

A = d(:)-d_obs(:);

loglike = -0.5*(A'*covinv)*A;
end