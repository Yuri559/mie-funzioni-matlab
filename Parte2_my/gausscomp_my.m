function [I,r] = gausscomp_my(a,b,N,n,fun)
%
%   [I,r] = gausscomp_my(a,b,N,n,fun)
%
%GAUSSCOMP Implementa la formula di quadratura di Gauss-Legendre composita
%a n punti (= ordine n)
%
%Input:
%   a → Estremo di integrazione sinistro
%   b → Estremo di integrazione destro
%   N → Numero di sottointervalli in cui si vuole suddividere [a,b]
%   n → ordine di quadratura (da 0 a 4)
%   fun → Funzione da integrare
%
%Output:
%   I → Valore approssimato dell'integrale
%   r → Grado di esattezza della formula

H =  (b - a) / N;

if n == 0
    y_ = 0;
    alpha_ = 2;
end
if n == 1
    y_ = [-1 / sqrt(3), 1/sqrt(3)]';
    alpha_ = [1, 1]';
end
if n == 2
    y_ = [-sqrt(15)/5, 0, sqrt(15)/5]';
    alpha_ = [5/9, 8/9, 5/9]';
end
if n == 3
    y_ = [-0.861136, -0.339981, 0.339981, 0.861136]';
    alpha_ = [0.3478, 0.652155, 0.652155, 0.3478]';
end
if n ~= 0 && n ~= 1 && n ~=2 && n ~= 3
    fprintf("ERRORE: n inserito non adatto, scegliere n = [0,1,2,3]\n")
end

x = a:H:b; % Vettore contenente gli estremi di tutti gli intervalli
v1 = x(1:N) + x(2:N+1); % Vettore contenente tutti gli a + b
v2 = x(2:N+1) - x(1:N); % Vettore contenente tutti i b - a

y = v1/2 + v2/2 .* y_;
alpha = v2/2 .* alpha_;

I = sum(sum(alpha.*fun(y)));

r = 2 * n + 1;

% Oss: questa funzione la ho scritta io, claude dice che è corretta è
% analiticamente equivalente a quella del prof, e per certi aspetti più
% efficiente