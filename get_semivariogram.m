function [yvals,xvals,ydistances,xdistances,XNS,YNS]=get_semivariogram(datz,dx,dy)
data = datz;
dimension = min(size(datz)); %is it 1D or 2D?


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

%% NEXT STEP IS TO DO THE SEMIVARIOGRAM ON DATZ AND NOT DATA
data = datz;
kx = (1:(size(data,2))-1)*dx; %Set of distances in horizontal direction
ky = (1:(size(data,1))-1)*dy; %Set of distances in vertical direction

yvals = [];
xvals = [];

xdistances = [];
ydistances = [];

XNS = [];
YNS = [];

for i = 1:length(ky)
    
   distance = ky(i);
   
   cols = size(data,2);
   N = (length(ky)-i+1)*cols;
   
   YNS = [YNS N];
   
   list = [];
   
   %make list of residuals between parameter pairs with the given distance
    for j = 1:N/cols
        L = data(j,:)-data(j+i,:);
        list = [list;L];
    end
    
   list = list.^2;
   val = sum(sum(list))/(2*N);
   
   %vals contains semivariogram values
   yvals = [yvals val];
   ydistances = [ydistances distance]; 
end

data = data';

for i = 1:length(kx)
    
   distance = kx(i);
   
   cols = size(data,2);
   N = (length(kx)-i+1)*cols;
   
   XNS = [XNS N];
   
   list = [];
   
   %make list of residuals between parameter pairs with the given distance
    for j = 1:N/cols
        L = data(j,:)-data(j+i,:);
        list = [list;L];
    end
    
   list = list.^2;
   val = sum(sum(list))/(2*N);
   
   %vals contains semivariogram values
   xvals = [xvals val];
   xdistances = [xdistances distance]; 
end

end