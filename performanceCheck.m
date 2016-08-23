%performance checks here

Ncomp = 100000;

tic
zeroMean = zeros(1, 8);
unitCov = eye(8);
for i = 1:Ncomp
   
    A = mvnrnd(zeroMean, unitCov);
    
end
t1 = toc



tic
for i = 1:Ncomp
   
    B = normrnd(0, 1, 1, 8);
    
end
t2 = toc
