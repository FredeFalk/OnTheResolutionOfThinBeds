function [pos, ress, resdummys, thick, theta, x0, length, inc] = get_TDEM_wedge(x0,depth,theta,dx,resolution,res,resdummy,length,inc) %takes intervals
xsize = dx/resolution+1;

        theta = drawrandom(theta);
        depth = drawrandom(depth);
        x0 = drawrandom(x0);
        length = drawrandom(length);

        res1 = res(1,:);
        res2 = res(2,:);
        
        DRO = drawrandom(res1);
        ress = [DRO;drawrandom(res2);DRO];


A = randperm(2);
ind = A(1);
inc = inc(ind);

thick = [];
for i = 1:xsize
    
    a = depth;
    
    if inc == 1
        if (i-1)*resolution > x0
            b = tand(theta)*((i-1)*resolution-x0);
            
            if (i-1)*resolution >x0+length
                b = 0;
            end      
        else
            b = 0;
        end
    end
    
    if inc == -1
        if (i-1)*resolution < x0
            b = tand(theta)*(x0-(i-1)*resolution);
            
            if (i-1)*resolution <x0-length
                b = 0;
            end
            
        else
            b = 0;
        end
    end
    thick = [thick;[a b]];
end

pos = linspace(0,dx-resolution,xsize);

resdummys = ress;
end