function [mods, output] = get_forwards(N,RES1,RES2,THICK1,THICK2,s,Fmode,res,wedge_t,wedge_r)
% N is the number of models, s is the chunksize
if Fmode == 1
   N = power(res,2); 
   s = N;
end

if Fmode == 2
   N = power(res,2); 
   s = N;
end

if Fmode == 3
    N = size(wedge_t,2)*N; 
    s = N;
end

if Fmode == 4
    N = N;
    s = N;
end

if Fmode == 100
    N = N;
    s = 2000;
end

G = waitbar(0,['Calculating ',num2str(N), ' forward responses.']);

remainder = mod(N,s);

J = remainder == 0;

if remainder == 0
    M = floor(N/s);
else
    M = ceil(N/s);
end

output = zeros(N,33,2);
mods = zeros(N,4);
Z = tic();

for k = 1:M
Progress = (k-1)*s;
target = N;
TimeSpent = toc(Z);
TimePerProgress = TimeSpent/Progress;
RemainingProgress = N-Progress;
RemainingTime = RemainingProgress*TimePerProgress;

PercProg = Progress/target;

waitbar(PercProg,G,['Calculating ',num2str(N), ' forward responses. (',num2str(RemainingTime/3600),'H left)']);

%Start by writing the model file
if Fmode == 100
    if k < M
        curchunksize = s;
        mods((k-1)*s+1:k*s,:) = write_modfile(curchunksize,RES1,RES2,THICK1,THICK2,Fmode,wedge_t((k-1)*s+1:k*s,:),wedge_r((k-1)*s+1:k*s,:));
    else
        if J
            curchunksize = s;
            mods((k-1)*s+1:k*s,:) = write_modfile(curchunksize,RES1,RES2,THICK1,THICK2,Fmode,wedge_t((k-1)*s+1:k*s,:),wedge_r((k-1)*s+1:k*s,:));
        else
            curchunksize = remainder;
            mods((end-remainder+1):end,:) = write_modfile(curchunksize,RES1,RES2,THICK1,THICK2,Fmode,wedge_t((end-remainder+1):end,:),wedge_r((end-remainder+1):end,:));
        end
    end

else
    mods((k-1)*s+1:k*s,:) = write_modfile(s,RES1,RES2,THICK1,THICK2,Fmode,wedge_t,wedge_r);   
end

%Run Forward Response on the model file
system('"AarhusInv64.exe" "tempo.mod"')

if Fmode == 100
    s = curchunksize;
end

for j = 1:s
%Read FWR file
ex = '00000';
ex = ex(1:end-length(num2str(j)));

filename = strcat('tempo',ex,num2str(j),'.fwr');
fileID2 = fopen(filename,'r');
dat = fscanf(fileID2,'%s');

Indcs = [];
Indcs2 = [];
len = 32;
for i = 0:len
    Indcs = [Indcs (2227+i*31):(2238+i*31)];
    Indcs2 = [Indcs2 (2218+i*31):(2226+i*31)];
end

fwr = dat(Indcs);
fwr = reshape(fwr,[12,len+1])';
fwr = str2num(fwr);

%% LOG TRANSFORM
fwr = real(log10(fwr));
%%

timesa = dat(Indcs2);
timesa = reshape(timesa,[9,len+1])';
timesa = str2num(timesa);

output((k-1)*s+j,:,1) = timesa;
output((k-1)*s+j,:,2) = fwr;

fclose(fileID2);
delete(filename)

end
end
delete(G)
end