function [nodes, weights] = naive_nodes_weights(n)
    syms x;
    P = legendre_poly(n);
    nodes = find_roots(P, n);
    P_prime = diff(P, x);
    weights = zeros(size(nodes));
    for i = 1:length(nodes)
        weights(i) = 2/((1 - nodes(i)^2)*(double(subs(P_prime, x, nodes(i))))^2);
    end
end

function P = legendre_poly(n)
    %computes symbolic Legrendre Polynomial 
    syms x;
    P0 = 1;
    P1 = x;
    if n == 0
        P = P0;
    elseif n == 1
        P = P1;
    else
        for k = 1:n-1
            P = ((2*k + 1) * x * P1 - k * P0) / (k + 1);
            P0 = P1;
            P1 = P;
        end
    end
end

function roots = find_roots(P, n)
    syms x;
    P_prime = diff(P, x); 
    roots = zeros(n, 1);
    initial = cos((2 * (1:n) - 1) * pi / (2 * n));

    % Newton's method to find roots
    tol = 1e-10; 
    max_iter = 100; 
    for i = 1:n
        x0 = initial(i); 
        for iter = 1:max_iter
            f = double(subs(P, x, x0)); 
            f_prime = double(subs(P_prime, x, x0)); 
            x1 = x0 - f / f_prime; 
            if abs(x1 - x0) < tol
                break;
            end
            x0 = x1;
        end
        roots(i) = x1;
    end
    roots = sort(roots);
end