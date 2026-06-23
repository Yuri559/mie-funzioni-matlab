function [x, X_storia, iter] = metodo_corde_universale_my(F, C, x0, tol, nmax)
%
% [x, X_storia, iter] = metodo_corde_universale_my(F, C, x0, tol, nmax)
%
% Metodo delle corde per sistemi non lineari: x = x - C^{-1}*F(x)
% C resta costante (es. Jacobiana valutata in x0) -> si fattorizza UNA
% volta sola con LU fuori dal ciclo, per efficienza.
%
% Parametri di ingresso:
%   F     - funzione anonima del sistema (es. @(x) A*x - b + ...)
%   C     - matrice di iterazione costante (es. Jacobiana in x0)
%   x0    - vettore di partenza (colonna)
%   tol   - tolleranza sull'incremento ||delta||
%   nmax  - numero massimo di iterazioni
%
% Parametri di uscita:
%   x        - soluzione finale
%   X_storia - matrice con tutte le iterate [x0, x1, ..., xN]
%   iter     - numero di iterazioni effettuate
%

n = length(x0);

% preallocazione: salva tutte le iterate
X_storia = zeros(n, nmax + 1);
X_storia(:, 1) = x0;
x = x0;

% fattorizzo C una sola volta fuori dal ciclo (cuore dell'efficienza)
[L, U, P] = lu(C);

for iter = 1:nmax
    % valuto F nel punto corrente
    Fx = F(x);

    % risolvo C*delta = Fx sfruttando LU (fwsub + bksub)
    delta = U \ (L \ (P * Fx));

    % aggiorno l'iterata
    x = x - delta;
    X_storia(:, iter + 1) = x;

    % criterio di arresto sulla norma dell'incremento
    if norm(delta) < tol
        X_storia = X_storia(:, 1:iter+1);   % pulisco colonne vuote
        return;
    end
end

end
