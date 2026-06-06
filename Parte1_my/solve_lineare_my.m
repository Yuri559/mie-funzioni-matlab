function [x, resrel, errstim, K] = solve_lineare_my(A, b, metodo)
%
% [x, resrel, errstim, K] = solve_lineare_my(A, b, metodo)
%
% Risolve il sistema lineare A*x = b con il metodo diretto scelto e
% restituisce residuo normalizzato ed errore relativo STIMATO a posteriori.
% Pensata per i temi d'esame Parte 1 (sistemi lineari, metodi diretti).
%
% INGRESSO:
%   A       matrice n x n (non singolare)
%   b       termine noto (vettore colonna)
%   metodo  intero che seleziona il metodo:
%             1 -> A\b           (backslash di Matlab, robusto)
%             2 -> fwsub         (sost. avanti: A triang. INFERIORE)
%             3 -> bksub         (sost. indietro: A triang. SUPERIORE)
%             4 -> lugauss       (LU senza pivoting)
%             5 -> lu            (LU con pivoting per RIGA: P*A = L*U)
%             6 -> lu totale     (pivoting TOTALE righe+colonne: P*A*Q = L*U)
%             7 -> thomas        (A tridiagonale)
%             8 -> chol          (Cholesky: A simmetrica def. positiva)
%
% USCITA:
%   x        soluzione numerica
%   resrel   residuo normalizzato  = ||b - A*x|| / ||b||
%   errstim  errore relativo stimato = K(A) * resrel       (stima a posteriori)
%   K        numero di condizionamento in norma 2 (cond(A))
%
% NOTA: errstim e' la stima  errrel <= K2(A)*resrel  (valida se dA ~ 0,
%       cioe' metodi con pivoting / Cholesky). E' una MAGGIORAZIONE.
%
% ----------------------------------------------------------------------

% --- controllo input base ---
b = b(:);                       % forzo vettore colonna
n = length(b);
if size(A,1) ~= n || size(A,2) ~= n
    error('Dimensioni A e b incompatibili')
end

% --- risoluzione in base al metodo ---
switch metodo

    case 1   % backslash
        x = A \ b;

    case 2   % fwsub (A triangolare inferiore)
        x = fwsub(A, b);

    case 3   % bksub (A triangolare superiore)
        x = bksub(A, b);

    case 4   % LU senza pivoting (lugauss)
        [L, U] = lugauss(A);
        y = fwsub(L, b);        % L*y = b
        x = bksub(U, y);        % U*x = y

    case 5   % LU con pivoting per riga (lu di Matlab): P*A = L*U
        [L, U, P] = lu(A);
        y = fwsub(L, P*b);      % L*y = P*b
        x = bksub(U, y);        % U*x = y

    case 6   % LU con pivoting TOTALE (righe + colonne): P*A*Q = L*U
        % NOTA: la sintassi a 4 output richiede matrice SPARSA -> sparse(A)
        [L, U, P, Q] = lu(sparse(A), [1.0 1.0]);  % thresh=1 -> pivot. parziale
        % P*A*Q = L*U  =>  x = Q*z,  L*U*z = P*b
        y = fwsub(L, P*b);      % L*y = P*b
        z = bksub(U, y);        % U*z = y
        x = Q * z;              % ripristino ordine colonne

    case 7   % Thomas (A tridiagonale)
        [~, ~, x] = thomas(A, b);

    case 8   % Cholesky (A SPD): A = R'*R
        R = chol(A);            % R triangolare superiore
        y = fwsub(R', b);       % R'*y = b
        x = bksub(R, y);        % R*x = y

    otherwise
        error('metodo deve essere un intero da 1 a 8')
end

x = full(x(:));   % garantisco vettore colonna pieno

% --- diagnostica: residuo, condizionamento, errore stimato ---
resrel  = norm(b - A*x) / norm(b);     % residuo normalizzato
K       = cond(A);                     % K2(A) = cond(A) in norma 2
errstim = K * resrel;                  % stima a posteriori: err <= K*resrel

% --- stampa riassuntiva ---
fprintf('--- solve_lineare_my (metodo %d) ---\n', metodo);
fprintf('  resrel  (||b-Ax||/||b||) = %.4e\n', resrel);
fprintf('  K(A)    (cond norma 2)   = %.4e\n', K);
fprintf('  errstim (K*resrel)       = %.4e\n', errstim);

end