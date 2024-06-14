function [pos rhos vps thick theta x0 length inc] = get_seismic_wedge(x0,depth,theta,dx,resolution,rho,vp,length,inc) %takes intervals
xsize = dx/resolution+1;

        theta = drawrandom(theta);
        depth = drawrandom(depth);
        x0 = drawrandom(x0);
        length = drawrandom(length);

        rho1 = rho(1,:);
        rho2 = rho(2,:);

        rhos = [drawrandom(rho1);drawrandom(rho2);drawrandom(rho1)];
        
        vp1 = vp(1,:);
        vp2 = vp(2,:);
        
        vps = [drawrandom(vp1);drawrandom(vp2);drawrandom(vp1)];

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

pos = linspace(0,dx,xsize);

end