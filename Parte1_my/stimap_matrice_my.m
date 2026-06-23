function [p, c] = stimap_matrice_my(xvect, xex)
%
% [p, c] = stimap_matrice_my(xvect, xex)
%
% Estensione di stimap (Politecnico di Milano, 04/04/2024) per iterate
% vettoriali. Accetta una matrice in cui ogni colonna e' un'iterata
% vettoriale x^(k) e calcola ordine e fattore di convergenza.
%
% Rispetto a stimap aggiunge:
%   - accetta xvect come matrice n x K (non solo vettore scalare)
%   - calcola l'errore come norma della differenza tra iterate successive
%     se xex non e' fornita, oppure come norma dell'errore vero ||x^(k)-xex||
%     se xex e' fornita
%
% Parametri di ingresso:
%   xvect  - matrice n x K (n incognite, K iterate)
%   xex    - (OPZIONALE) soluzione esatta vettoriale (n x 1)
%            se assente: usa differenze successive come stima dell'errore
%
% Parametri di uscita:
%   p      - vettore delle stime dell'ordine di convergenza
%   c      - vettore delle stime del fattore di abbattimento
%
% Uso tipico:
%   [p, c] = stimap_matrice_my(xvect);          % senza soluzione esatta
%   [p, c] = stimap_matrice_my(xvect, x_true);  % con soluzione esatta
%

if nargin < 2
    xex = [];
end

[~, K] = size(xvect);

% pre-allocazione
p   = zeros(1, K);
c   = zeros(1, K);
err = zeros(1, K);

% --- 1. Calcolo sequenza errori ---
if isempty(xex)
    % senza soluzione esatta: norma diff. tra iterate successive
    for j = 1:K-1
        err(j) = norm(xvect(:, j+1) - xvect(:, j));
    end
    limite = K - 1;
else
    % con soluzione esatta: errore vero
    xex = xex(:);
    for j = 1:K
        err(j) = norm(xvect(:, j) - xex);
    end
    limite = K;
end

% --- 2. Calcolo ordine p e fattore c ---
% servono almeno 3 errori -> ciclo parte da i=3
for i = 3:limite
    e_k1 = err(i);      % errore passo corrente
    e_k  = err(i-1);    % errore passo precedente
    e_km = err(i-2);    % errore due passi fa

    if (e_k1 == 0 || e_k == 0 || e_km == 0)
        disp('Attenzione: convergenza esatta raggiunta (errore = 0).');
        break;
    else
        num  = log(e_k1 / e_k);
        den  = log(e_k  / e_km);
        p(i) = num / den;
        c(i) = e_k1 / (e_k^p(i));
    end
end

% --- 3. Stampa ---
indici = find(p ~= 0);
if ~isempty(indici)
    ultimo = indici(end);
    fprintf(' Ordine stimato       : %12.8f \n', p(ultimo));
    fprintf(' Fattore di riduzione : %12.8f \n', c(ultimo));
else
    disp(' Numero di iterazioni insufficiente per stimare l''ordine!');
end

end