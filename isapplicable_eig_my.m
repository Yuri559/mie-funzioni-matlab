function isapplicable_eig_my(A, s)
%
% isapplicable_eig_my(A)
% isapplicable_eig_my(A, s)
%
% Verifica l'applicabilita' dei metodi delle potenze a partire dallo
% spettro di A (calcolato con eig). Stampa SI/NO per:
%   - Potenze dirette (eigpower)      -> lambda_1 (modulo MAX)
%   - Potenze inverse (invpower)      -> lambda_n (modulo MIN)
%   - Potenze inverse con shift s     -> autovalore piu' vicino a s
%                                        (solo se 's' viene passato)
%
% Parametri di ingresso:
% A   Matrice quadrata (n x n)
% s   (opzionale) shift per il metodo delle potenze inverse con shift
%
% NOTA: usa una tolleranza relativa per dichiarare "uguali" due moduli,
%       cosi' coglie anche coppie complesse coniugate / autovalori +-lambda.
%
%                                              versione _my (open book)
%

lambda = eig(A);
n = length(lambda);
m = sort(abs(lambda), 'descend');   % moduli ordinati: m(1) >= ... >= m(n)

% tolleranza relativa per confrontare due moduli
tol = 1e-10 * max(m(1), 1);

fprintf('\n--- Spettro di A ---\n');
disp(lambda);
fprintf('Moduli ordinati (decrescente):\n');
disp(m);

% ---- POTENZE DIRETTE: serve |lambda_1| > |lambda_2| stretto ----
fprintf('\n[Potenze DIRETTE - eigpower] lambda max modulo:\n');
if abs(m(1) - m(2)) > tol
    fprintf('  APPLICABILE  (|l1|=%.4g > |l2|=%.4g)\n', m(1), m(2));
else
    fprintf('  NON applicabile (|l1|=|l2|=%.4g: modulo max NON unico)\n', m(1));
end

% ---- POTENZE INVERSE: serve |lambda_{n-1}| > |lambda_n| stretto ----
fprintf('[Potenze INVERSE - invpower]  lambda min modulo:\n');
if abs(m(n-1) - m(n)) > tol
    fprintf('  APPLICABILE  (|l_{n-1}|=%.4g > |l_n|=%.4g)\n', m(n-1), m(n));
else
    fprintf('  NON applicabile (|l_{n-1}|=|l_n|=%.4g: modulo min NON unico)\n', m(n));
end

% ---- POTENZE INVERSE CON SHIFT: serve |lk-s| < |li-s| stretto ----
if nargin == 2
    d = sort(abs(lambda - s), 'ascend');   % distanze da s
    tols = 1e-10 * max(d(end), 1);
    fprintf('[Potenze INVERSE con SHIFT s=%.4g] autovalore piu'' vicino a s:\n', s);
    if d(1) < tols
        fprintf('  NON applicabile (s coincide con un autovalore: A-sI singolare)\n');
    elseif abs(d(2) - d(1)) > tols
        % trova quale autovalore e' il piu' vicino
        [~, idx] = min(abs(lambda - s));
        fprintf('  APPLICABILE  -> approssima lambda = %s\n', num2str(lambda(idx)));
        fprintf('               (dist min=%.4g < seconda dist=%.4g)\n', d(1), d(2));
    else
        fprintf('  NON applicabile (due autovalori equidistanti da s: dist=%.4g)\n', d(1));
    end
end

fprintf('\n');

return