function [I,r] = gausscomp_my(a,b,N,n,fun)
%
%   [I,r] = gausscomp_my(a,b,N,n,fun)
%
%GAUSSCOMP Implementa la formula di quadratura di Gauss-Legendre composita
%a (n+1) punti (= ordine n)
%
%Input:
%   a → Estremo di integrazione sinistro
%   b → Estremo di integrazione destro
%   N → Numero di sottointervalli in cui si vuole suddividere [a,b]
%   n → ordine di quadratura (intero ≥ 0)
%   fun → Funzione da integrare
%
%Output:
%   I → Valore approssimato dell'integrale
%   r → Grado di esattezza della formula

H = (b - a) / N;

if ~isscalar(n) || n < 0 || n ~= floor(n)
    fprintf("ERRORE: n deve essere un intero >= 0\n")
    I = NaN; r = NaN;
    return
end

% Golub-Welsch: nodi e pesi su [-1,1] con (n+1) punti
m = n + 1;
if m == 1
    y_ = 0;
    alpha_ = 2;
else
    k = (1:m-1)';
    beta = k ./ sqrt(4*k.^2 - 1);     % coeff. ricorrenza Legendre
    J = diag(beta,1) + diag(beta,-1); % matrice di Jacobi simmetrica
    [V,D] = eig(J);
    [y_, idx] = sort(diag(D));        % nodi = autovalori
    alpha_ = 2 * (V(1,idx).^2)';      % pesi
end

x = a:H:b;                 % estremi di tutti gli intervalli
v1 = x(1:N) + x(2:N+1);    % tutti gli a + b
v2 = x(2:N+1) - x(1:N);    % tutti i b - a
y = v1/2 + v2/2 .* y_;
alpha = v2/2 .* alpha_;
I = sum(sum(alpha.*fun(y)));
r = 2 * n + 1;