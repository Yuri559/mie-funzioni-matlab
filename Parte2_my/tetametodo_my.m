function [t_h,u_h,iter_pf] = tetametodo_my(f,t_max,y0,h,teta,A)
%
% [t_h,u_h,iter_pf] = tetametodo_my(f,t_max,y0,h,teta)
% [t_h,u_h,iter_pf] = tetametodo_my(f,t_max,y0,h,teta,A)
%
% Risolve il problema di Cauchy (anche sistemi di EDO del 1° ordine)
% col Teta-Metodo.
%
% DUE MODALITA':
%  - PROBLEMA LINEARE  y' = A*y + g(t):  passa anche A (6° input)
%       -> implementazione EFFICIENTE: (I - teta*h*A) fattorizzata UNA
%          volta con LU, poi solo sostituzioni avanti/indietro nel ciclo.
%          (in questo caso f deve essere f = @(t,y) A*y + g(t) )
%  - PROBLEMA NON LINEARE: NON passare A
%       -> iterazioni di punto fisso ad ogni passo (come versione base).
%
% Input:
%   f       -> funzione f(t,y), y colonna (N_eq x 1)
%   t_max   -> estremo dell'intervallo
%   y0      -> dato iniziale (colonna N_eq x 1)
%   h       -> passo temporale
%   teta    -> 0 -> EA ; 1 -> EI ; 0.5 -> CN
%   A       -> (OPZIONALE) matrice del sistema lineare. Se presente,
%              attiva il ramo lineare efficiente con LU.
% Output:
%   t_h     -> istanti temporali discreti
%   u_h     -> soluzione approssimata (N_eq x N_istanti)
%   iter_pf -> n° iterazioni di punto fisso ad ogni passo
%              (vettore di zeri nel ramo lineare: non si itera)

t0  = 0;
t_h = t0:h:t_max;

N_istanti = length(t_h);
y0   = y0(:);             % forza colonna
N_eq = length(y0);

u_h = zeros(N_eq,N_istanti);
u_h(:,1) = y0;

iter_pf = zeros(1,N_istanti);

% ===== scelta del ramo =====
if nargin >= 6 && ~isempty(A)
    % ---------- RAMO LINEARE EFFICIENTE (LU) ----------
    % g(t) = f(t,0) perche' f = A*y + g(t)  ->  f(t,0) = g(t)
    g = @(t) f(t, zeros(N_eq,1));

    % matrice di iterazione, COSTANTE: fattorizzo UNA sola volta
    B = eye(N_eq) - teta*h*A;
    [L,U,P] = lu(B);

    for n = 1:(N_istanti-1)
        % termine noto: u^n + h[(1-teta)(A u^n + g(t_n)) + teta g(t_{n+1})]
        b = u_h(:,n) + h*( (1-teta)*( A*u_h(:,n) + g(t_h(n)) ) ...
                         +    teta * g(t_h(n+1)) );
        y = fwsub(L, P*b);            % sostituzione in avanti
        u_h(:,n+1) = bksub(U, y);     % sostituzione all'indietro
    end

else
    % ---------- RAMO NON LINEARE (punto fisso) ----------
    nmax = 100;
    toll = 1e-5;

    for n = 1:(N_istanti-1)

        phi = @(y) u_h(:,n) + h*( (1-teta)*f( t_h(n),   u_h(:,n) ) ...
                                +    teta *f( t_h(n+1), y        ) );

        xv  = u_h(:,n);
        err = 1 + toll;
        it  = 0;
        while (it < nmax && err > toll)
            xn  = phi(xv);
            xn  = xn(:);
            err = norm(xn - xv);
            it  = it + 1;
            xv  = xn;
        end

        u_h(:,n+1)   = xv;
        iter_pf(n+1) = it;
    end
end

end