function [ value ] = PearsonVII( x,amp,x0,gamma,m )
    value = amp * gamma^(2*m) ./ ((((2^(1/m))-1)*(x-x0).^2 + gamma^2).^m);
end