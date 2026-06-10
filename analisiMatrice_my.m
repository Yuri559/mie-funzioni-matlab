function analisiMatrice_my(A, P, al, k)
%% ============================================
%  ANALISI PROPRIETÀ MATRICE A
%  ============================================
%
%   analisiMatrice_my(A, P, al, k)
%
%  INPUT:
%    A  - matrice da analizzare (obbligatorio)
%    P  - precondizionatore per Richardson (opzionale)
%    al - alpha del metodo di Richardson (opzionale)
%    k  - numero di iterazioni da stampare (opzionale)
%
%  Esempi d'uso:
%    analisiMatrice(A)              % solo proprietà della matrice
%    analisiMatrice(A, P, al, k)    % anche raggio spettrale Richardson

n = size(A,1);
I = eye(n);   % matrice identità

%% CONDIZIONE NECESSARIA E SUFFICIENTE (per LU senza pivoting)
if det(A) == 0
    fprintf('matrice singolare!\n')
end

con_nc = 1;
ii = 1;
while ii <= n && con_nc == 1
    if det(A(1:ii,1:ii)) == 0
        con_nc = 0;          % minore principale nullo
    end
    ii = ii + 1;
end

if con_nc == 1
    disp('condizione necessaria e sufficiente soddisfatta')
else
    disp('condizione necessaria e sufficiente non soddisfatta')
end

%% SDP - Simmetrica Definita Positiva
% verifico simmetria
if isequal(A, A')
    disp('matrice simmetrica')
else
    disp('matrice non simmetrica')
end

% verifico definita positiva
aut = eig(A);                % autovalori
if all(aut > 0)
    disp('matrice definita positiva')
else
    disp('matrice non definita positiva')
end

%% DOMINANZA STRETTA PER RIGHE
domr = 1;
ii = 1;
while domr == 1 && ii <= n
    if abs(A(ii,ii)) <= sum(abs(A(ii,:))) - abs(A(ii,ii))
        domr = 0;
    end
    ii = ii + 1;
end

if domr == 1
    disp('dominanza stretta per righe')
else
    disp('non dominanza stretta per righe')
end

%% DOMINANZA STRETTA PER COLONNE
domc = 1;
ii = 1;
while domc == 1 && ii <= n
    if abs(A(ii,ii)) <= sum(abs(A(:,ii))) - abs(A(ii,ii))
        domc = 0;
    end
    ii = ii + 1;
end

if domc == 1
    disp('dominanza stretta per colonne')
else
    disp('non dominanza stretta per colonne')
end

%% RAGGIO SPETTRALE - RICHARDSON CON ALPHA
% eseguito solo se vengono passati P e al
if nargin >= 3
    rho = max(abs(eig(I - al*(P\A))));   % raggio spettrale matrice di iterazione

    fprintf('alpha = %g \n', al);
    fprintf('rho   = %g \n', rho);

    if nargin >= 4
        fprintf('#iterazioni: %g \n', k);
    end

    % criterio di convergenza
    if rho < 1
        disp('il metodo converge (rho < 1)')
    else
        disp('il metodo NON converge (rho >= 1)')
    end
end

end