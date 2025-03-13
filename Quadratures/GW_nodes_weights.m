function [nodes, weights] = GW_nodes_weights(n)
   % Golub-Welsch algorithm
   % GLINT  Gauss-Legendre numerical integration.
   % Fundamentals of Numerical Computation, By Tobin A. Driscoll and Richard J. Braun
   % https://github.com/tobydriscoll/fnc-extras/blob/master/fnc/glint.m
    beta = 0.5 ./ sqrt(1-(2*(1:n-1)).^(-2));
    T = diag(beta, 1) + diag(beta, -1);         
    [V, D] = eig(T);                                
    nodes = diag(D);                            
    weights = 2 * (V(1, :).^2)';                
    [nodes, idx] = sort(nodes);
    weights = weights(idx);
end
