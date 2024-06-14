function cov=get_covariance_matrix(fwr,covarmodel,dy,a,NT)
%The covariance matrix is n-by-n for n samples, and
%is going to be constructed given an input covariance model (gaussian)
rows = size(fwr,1);
cols = size(fwr,2);
hold off
%covarmodel is the covariance model
N = length(fwr(:));
mainmat = zeros(N,N);
cov = zeros(N,N);

%set horizontal spacing, inf when determining noise model, other number
%when using the model to produce noise.

if a==0;
dx = 1E20; %This can be modified and doesnt have to be the true distance, to simulate anisotrophy.
else
dx = 3;
end

for i = 1:N
    for j = 1:N
        %Translates from linear index in cov matrix to subscript indexes in
        %the spatial model space
        
        [MRow1,MCol1] = ind2sub([rows cols],i);
        [MRow2,MCol2] = ind2sub([rows cols],j);
        

        dist = sqrt(power(dy*(MRow2-MRow1),2)+power(dx*(MCol2-MCol1),2));
        
        if NT == 0
            dist = inf;
        end
        if i == j
           dist = 0; 
        end
        cov(i,j) = covarmodel(dist);
    end
end
end