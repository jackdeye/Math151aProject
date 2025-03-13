f1 = @(x) 1./(1+1000*(x-0.5).^2); %[0,1] =  Analytical solution of 
f2 = @(x) 1./x.*sqrt((1+x)./(1-x)).*log((2*x.^2+2*x+1)./(2*x.^2-2*x+1)); %[-1,1] = 4pi arccot sqrt(phi) obviously
f3 = @(x) sqrt(abs(x)); %[-1,1] = 4/3
f4 = @(x) 1./(1+x.^2); %[0,1] = pi/4
f5 = @(x) cos(x); % [-pi/2,pi/2] = 2
e1 = atan(5*sqrt(10))/(5*sqrt(10));
e2 = 4*pi*acot(sqrt(( 1 + sqrt(5) ) / 2));
e3 = 4/3;
e4 = pi/4;
e5 = 2;
eval = @(n) gauss_legendre_quadrature(f5,-pi/2,pi/2,n);

num = 80;
skip = 2;
begin = 2;
x = 1:skip:num;
values = zeros(num/skip,1);
for n = begin:skip:num
    diff = e5-eval(n);
    values(n/skip) = diff;
    disp(n)
    if abs(diff) < eps
        break;
    end
end
figure;
values = values(:);
plot(x, values(:), 'Color', 'magenta');
figure;
plot(log(x), log(values(:)), 'Color', 'green');