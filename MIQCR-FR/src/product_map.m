function T = product_map(n, i, j)
% vec(Y)=e00+T*[x;y], with arrow equations and symmetric products.
% Each x_i fills Y_0i, Y_i0 and Y_ii; each y_ij fills Y_ij and Y_ji.
d = n+1;
k = numel(i);
a = (1:n)';
T = sparse([1+d*a; 1+a; 1+a+d*a; i+1+d*j; j+1+d*i], ...
    [a; a; a; n+(1:k)'; n+(1:k)'], 1, d*d, n+k);
end
