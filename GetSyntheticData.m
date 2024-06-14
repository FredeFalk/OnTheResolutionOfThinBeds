function [cov,noise,data,cleandata]=GetSyntheticData(SNR,trues,NT,type)
        
true_depth = trues(1);
true_dip = trues(2);
true_x0 = trues(9);
true_length = trues(10);
profile_length = trues(11);            %meters
resolution = trues(12);                 %Distance between traces, meters
true_inc = trues(13); 

switch type
    case 'seis'
        true_rhos = trues(3:5);
        true_vps = trues(6:8);

        %Fetch synthetic data
        [noise, data, cov, cleandata] = get_synthetic_seismic([true_x0 true_x0],[true_depth true_depth],[true_dip true_dip],profile_length,resolution,[true_rhos(1) true_rhos(1);true_rhos(2) true_rhos(2);true_rhos(3) true_rhos(3)],[true_vps(1) true_vps(1);true_vps(2) true_vps(2);true_vps(3) true_vps(3)],[true_length true_length],[true_inc true_inc],SNR,NT);
        
        noise = reshape(noise,size(data));
    case 'TEM'
        true_res = trues(6:8);

        %Fetch synthetic data
        [noise, data, cov, cleandata] = get_synthetic_TEM([true_x0 true_x0],[true_depth true_depth],[true_dip true_dip],profile_length,resolution,[true_res(1) true_res(1);true_res(2) true_res(2);true_res(3) true_res(3)],[true_length true_length],[true_inc true_inc],SNR,NT);
        
        noise = reshape(noise,size(data));
end
end