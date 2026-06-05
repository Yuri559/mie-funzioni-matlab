function [s, lambda_min] = scelta_shift_my(A, k, perturb)
%
% [s, lambda_min] = scelta_shift_my(A, k, perturb)
%
% Calcola gli SHIFT da usare con invpowershift.m per approssimare i k
% autovalori di MODULO MINIMO di A (reali o complessi coniugati).
%
% Ogni shift viene posto VICINO all'autovalore target ma NON uguale, cosi':
%   - A - s*I resta invertibile (la LU dentro invpowershift non esplode);
%   - l'autovalore target e' quello UNIVOCAMENTE piu' vicino a s
%     (condizione di convergenza del metodo delle potenze inverse con shift).
%
% Lo scostamento e' una piccola frazione del "gap" tra l'autovalore target
% e il suo vicino piu' prossimo -> garantito corretto in ogni caso.
%
% INPUT:
%   A        matrice quadrata (n x n)
%   k        quanti autovalori di modulo minimo cercare (default 1)
%   perturb  frazione del gap usata per scostare s dall'autovalore,
%            deve stare in (0, 0.5)  (default 1e-2 = 1%);
%            piu' piccola -> s piu' vicino al target
%
% OUTPUT:
%   s           vettore (k x 1) degli SHIFT -> da passare a invpowershift
%   lambda_min  vettore (k x 1) degli autovalori VERI di modulo minimo
%               (solo riferimento per giustificare la scelta all'esame:
%                NON usarli come shift, A - lambda*I sarebbe singolare!)
%
% NOTA: matrice reale con autovalori minimi complessi -> gli shift escono
%       complessi (uno per ogni coniugato). E' corretto, deve essere cosi'.
%
% ESEMPIO D'USO:
%   [s, lm] = scelta_shift(A, 2);          % 2 autovalori minimi in modulo
%   lambda = zeros(numel(s),1);
%   for j = 1:numel(s)
%       lambda(j) = invpowershift(A, s(j), 1e-6, 1000, ones(size(A,1),1));
%   end
%   % lambda contiene i 2 autovalori di modulo minimo approssimati
%
%                                                              versione _MY
%

    if nargin < 2 || isempty(k),       k = 1;          end
    if nargin < 3 || isempty(perturb), perturb = 1e-2; end

    [n, m] = size(A);
    if n ~= m,                  error('Solo per matrici quadrate'); end
    if k > n,                   error('k non puo'' superare n'); end
    if perturb <= 0 || perturb >= 0.5
        error('perturb deve stare in (0, 0.5)');
    end

    % 1) autovalori VERI, ordinati per modulo crescente
    lambda = eig(A);
    [~, idx] = sort(abs(lambda));
    lambda = lambda(idx);

    lambda_min = lambda(1:k);   % i k piu' piccoli in modulo (valori veri)

    % 2) costruisco uno shift per ciascun target
    s = zeros(k, 1);
    for j = 1:k
        lt = lambda(j);
        % distanza dal vicino piu' prossimo (escludo se stesso)
        dist = abs(lambda - lt);
        dist(j) = Inf;
        d = min(dist);
        if d == 0                          % autovalore ripetuto: fallback
            d = perturb * max(abs(lambda));
            if d == 0, d = perturb; end    % matrice nulla: fallback estremo
        end
        % scosto di una piccola frazione del gap: s univocamente piu' vicino a lt
        s(j) = lt + perturb * d;
    end

end