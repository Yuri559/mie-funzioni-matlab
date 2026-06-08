function u = uj_diff_trasp_my(mu, eta, h, u0, uNp1, N)
%
%   u = uj_diff_trasp_my(mu, eta, h, u0, uNp1, N)
%
% DIFF_FINITE_DT_MY  Differenze finite centrate per diffusione-trasporto 1D.
%   Risolve  -mu*u'' + eta*u' = 0  su (0,1) con condizioni di Dirichlet
%   u(0)=u0, u(N+1)=uNp1, usando lo schema (9.26) pag. 209.
%
%   Funziona sia con valori numerici sia con mu, eta, h, u0, uNp1 passati
%   come variabili simboliche (syms): in tal caso restituisce le uj come
%   espressioni simboliche in funzione dei parametri.
%
%   INPUT:
%     mu    -> coeff. diffusione (davanti a -u'')
%     eta   -> coeff. trasporto  (davanti a  u')
%     h     -> passo di griglia  (h = 1/(N+1))
%     u0    -> valore al bordo sinistro  u(x0)   = u(0)
%     uNp1  -> valore al bordo destro    u(x_N+1)= u(1)
%     N     -> numero di nodi INTERNI (incognite). (N+1) sottointervalli
%
%   OUTPUT:
%     u -> vettore [u0; u1; ...; uN; uN+1]  (tutti i nodi, bordi inclusi)
%
%   Stencil (9.26), riga j-esima:
%     (-mu/h^2 - eta/(2*h)) u_{j-1}
%   + ( 2*mu/h^2 )          u_{j}
%   + (-mu/h^2 + eta/(2*h)) u_{j+1}  = 0

    % coefficienti dello stencil
    a = -mu/h^2 - eta/(2*h);   % sotto-diagonale (u_{j-1})
    d =  2*mu/h^2;             % diagonale       (u_{j})
    c = -mu/h^2 + eta/(2*h);   % sopra-diagonale (u_{j+1})

    % matrice tridiagonale A (N x N) e termine noto b
    A = d*eye(N) + c*diag(ones(N-1,1),1) + a*diag(ones(N-1,1),-1);
    b = zeros(N,1);

    % se qualche input e' simbolico, lavoro in simbolico
    if isa(mu,'sym') || isa(eta,'sym') || isa(h,'sym') || ...
       isa(u0,'sym') || isa(uNp1,'sym')
        A = sym(A);
        b = sym(b);
    end

    % inserisco le condizioni al contorno nel termine noto
    b(1)   = b(1)   - a*u0;     % il nodo j=1 usa u0
    b(end) = b(end) - c*uNp1;   % il nodo j=N usa u_{N+1}

    % risolvo il sistema interno
    uint = A\b;

    if isa(A,'sym')
        uint = simplify(uint);
    end

    % ricompongo il vettore completo coi bordi
    u = [u0; uint; uNp1];
end