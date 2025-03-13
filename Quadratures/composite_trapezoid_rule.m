function I = composite_trapezoid_rule(f, a, b, n)
    % COMPOSITE_TRAPEZOID_RULE Calculates integral approximation based off of
    % the composite trapezoid rule, given function handle f, lower bound a,
    % upper bound b, subintervals n
h = (b - a) / n;
x = a:h:b;
I = h / 2 * (f(a) + 2*sum(f(x(2:end-1))) + f(b));
end

