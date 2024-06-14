function [times, Fwr] = seismic_forward(thicknesses,rhos,vps,sampf)

samp = 0.1; %Initial sampling rate in ms (for forward response accuracy)
 %sampf final sampling rate (output)
time = 100; %Sampling time in ms
N = time/samp+1; %number of samples

samps = linspace(0,time,N);


TWT = cumsum(2000*thicknesses./vps(1:end-1)); %Calculate two-way-times to layer boundaries
IMPs = rhos.*vps;
IMPC = diff(IMPs);
IMPS = [IMPs(1)+IMPs(2) IMPs(2)+IMPs(3)]';
INDX = round(TWT/samp);
List = zeros(1,N);
List2 = zeros(1,N);

if INDX(1) == INDX(2)
    INDX = INDX(2);
    IMPC = sum(IMPC);
    IMPS = IMPs(1)+IMPs(3);
end

RS = IMPC./IMPS;

List(INDX) = RS;

R=conv(ricker(50,40/(samp)+1,samp/1000),List); %Reflection profile by convolution with a Ricker Wavelet

Fwr = R;
Fwr = Fwr';

times = samps;
Fwr = Fwr(20/samp+1:end-20/samp);

times = times(1:sampf/samp:end);
Fwr = Fwr(1:sampf/samp:end);
end