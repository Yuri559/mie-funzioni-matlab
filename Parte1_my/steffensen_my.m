function [x, N] = steffensen_my(f, x0, tol, max_iter)
%
% [x, N] = steffensen_my(f, x0, tol, max_iter)
%
% Trova lo zero di una funzione usando il metodo di Steffensen.
% Convergenza quadratica (come Newton) ma senza calcolare la derivata.
% Formula: x_{k+1} = x_k - f(x_k)^2 / (f(x_k + f(x_k)) - f(x_k))
%
% Parametri di ingresso:
%   f        - handle della funzione
%   x0       - punto di partenza iniziale
%   tol      - tolleranza (criterio arresto: |f(x)| < tol)
%   max_iter - numero massimo di iterazioni
%
% Parametri di uscita:
%   x        - approssimazione dello zero
%   N        - numero di iterazioni effettuate
%

x = x0;
N = 0;

while abs(f(x)) >= tol && N < max_iter
    fx          = f(x);
    denominatore = f(x + fx) - fx;

    if denominatore == 0
        error('Denominatore nullo: impossibile continuare l''iterazione.');
    end

    x = x - (fx^2) / denominatore;
    N = N + 1;
end

if N == max_iter && abs(f(x)) >= tol
    warning('Raggiunto il numero massimo di iterazioni (%d) senza convergenza.', max_iter);
else
    fprintf(' Convergenza al passo k : %d \n', N);
    fprintf(' Radice calcolata       : %-12.8f \n', x);
end

end
