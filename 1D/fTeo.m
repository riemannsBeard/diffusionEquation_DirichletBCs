function [ uTeo ] = fTeo( x, t, N )
    
uTeo = 0;

for n = 1:1:N
    uTeo = uTeo + (8/pi^2)*(1/n^2)*sin(0.5*n*pi)*sin(n*pi*x)*exp(-n*n*pi*pi*t);
end

end

