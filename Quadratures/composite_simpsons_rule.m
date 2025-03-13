function I = composite_simpsons_rule(f, a, b, n)
    % COMPOSITE_SIMPSONS_RULE Calculates integral approximation based off of
    % the composite simpson's rule, given function handle f, lower bound a,
    % upper bound b, subintervals n (must be even)
if mod(n, 2) ~= 0
    error("n Must be even");
else
    h = (b - a) / n;
    x = a:h:b;

    odds = 2:2:n;
    evens = 3:2:n-1;

    I = h / 3 * (f(x(1)) + f(x(n+1)) + 4 * sum(f(x(odds))) + 2 * sum(f(x(evens))));
end

