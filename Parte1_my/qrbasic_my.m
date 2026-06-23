function [D, lambda2_hist, rate] = qrbasic_my(A, tol, nmax)
%
%   [D, lambda2_hist, rate] = qrbasic_my(A, tol, nmax)
%
% Estensione di qrbasic (Politecnico di Milano, 04/04/2024).
% Stesse iterazioni QR, stessi parametri di ingresso, stesso D in uscita.
%
% Rispetto a qrbasic aggiunge:
%   lambda2_hist  - storia delle approssimazioni di T(2,2) ad ogni
%                   iterazione (cambia (2,2) con (k,k) per il k-esimo
%                   autovalore). Utile per stimap e analisi convergenza.
%   rate          - stima dell'ordine di convergenza per lambda_2,
%                   calcolata con regressione log-log sulle differenze
%                   successive (mediana per robustezza).
%
% Parametri di ingresso:
%   A     - matrice quadrata reale (n x n)
%   tol   - tolleranza sul criterio d'arresto
%   nmax  - numero massimo di iterazioni
%
% Parametri di uscita:
%   D            - vettore n x 1 degli autovalori approssimati
%   lambda2_hist - storia approssimazioni secondo autovalore
%   rate         - stima ordine di convergenza (NaN se non stimabile)
%

    [n, m] = size(A);
    if n ~= m, error('La matrice deve essere quadrata'); end

    % --- Inizializzazione ---
    k            = 1;
    lambda2_hist = NaN(nmax, 1);
    T            = A;
    niter        = 0;
    test         = max(max(abs(tril(T, -1))));

    % --- Iterazioni QR: T = R*Q ---
    while niter < nmax && test > tol
        [Q, R] = qr(T);
        T      = R * Q;
        niter  = niter + 1;
        test   = max(max(abs(tril(T, -1))));

        % salva approssimazione lambda_2 (cambia (2,2) per altri autovalori)
        if n >= 2
            lambda2_hist(k) = T(2,2);
        else
            lambda2_hist(k) = NaN;
        end
        k = k + 1;
    end

    % --- Tronca storia ---
    lambda2_hist = lambda2_hist(1:k-1);

    if niter >= nmax
        fprintf('Il metodo non converge nel massimo numero di iterazioni.\n');
    else
        fprintf('Il metodo converge in %d iterazioni.\n', niter);
    end

    D = diag(T);

    % --- Stima ordine convergenza (log-log su differenze successive) ---
    % modello: errore_{k+1} ~ C * errore_k^p
    rate = NaN;
    if length(lambda2_hist) >= 4
        lambda_true_est = lambda2_hist(end);
        e   = abs(lambda2_hist - lambda_true_est);
        idx = find(e > 0);
        if length(idx) >= 3
            ek    = e(idx);
            ek1   = ek(2:end);
            ek0   = ek(1:end-1);
            p_est = log(ek1) ./ log(ek0);    % stime locali di p
            rate  = median(p_est(isfinite(p_est)));   % mediana per robustezza
        end
    end
end