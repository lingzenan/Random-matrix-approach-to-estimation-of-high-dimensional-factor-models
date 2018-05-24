% ASD
% arg: c = N/T
%        b: the global mean-reversion
%rtn: pr: asymptotic spectral distribiton 
function pr = ASD_AR(b,c)
a = sqrt(1-b^2);
syms m
t = 0.01:0.1:9.01; %supp
num = length(t);
d = zeros(4,num);
for k = 1:num
    z = t(k)+0.0001i;
    d(:,k) = solve(a^4*c^2*m^4+2*a^2*c*(-(1+b^2)*z+a^2*c)*m^3+((1-b^2)^2*z^2-2*a^2*c*(1+b^2)*z+(c^2-1)*a^4)*m^2-2*a^4*m-a^4,'m');
    d = double(d);
end
pr1 = zeros(4,num);
prr = zeros(4,num);
for j = 1:num
    for k = 1:4
        if (-1/pi*imag((d(k,j)+1)/(t(j)+0.0001i))) > 0.00001 % choose root
            prr(k,j) =  d(k,j);
            pr1(k,j) = -1/pi*imag((d(k,j)+1)/(t(j)+0.0001i));
        end
    end
end
pr = sum(pr1,1);
end
