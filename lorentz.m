function [ value ] = lorentz( x,amp,x0,gamma )
    value = amp * gamma^2 ./ ((x-x0).^2 + gamma^2);
end