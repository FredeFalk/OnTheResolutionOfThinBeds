function ra = drawrandom(interval)
ra = interval(:,1)+rand(size(interval,1),1).*(interval(:,2)-interval(:,1));
end