function [u, x] = risolvi_bvp(a, b, h, alfa_fun, beta_fun, gamma_fun, f_fun, tipo_sx, val_sx, tipo_dx, val_dx)
%
%   [u, x] = risolvi_bvp(a, b, h, alfa_fun, beta_fun, gamma_fun, f_fun, tipo_sx, val_sx, tipo_dx, val_dx)
%
% Risolve un problema ai limiti (BVP) lineare del secondo ordine.
% nell'intervallo [a, b] utilizzando il metodo delle differenze 
% finite centrate del secondo ordine.
% Il modello matematico generale risolto è:
%       -alfa(x)*u''(x) + beta(x)*u'(x) + gamma(x)*u(x) = f(x)
%
% INPUT:
%   a, b     : Scalari. Estremi dell'intervallo di integrazione [a, b].
%   h        : Scalare. Passo di discretizzazione spaziale (ampiezza della griglia).
%   alfa_fun : Function handle del coefficiente della derivata seconda, alfa(x).
%   beta_fun : Function handle del coefficiente della derivata prima, beta(x).
%   gamma_fun: Function handle del coefficiente del termine noto lineare, gamma(x).
%   f_fun    : Function handle del termine sorgente/noto, f(x).
%   tipo_sx  : Stringa. Tipo di condizione al contorno a sinistra (x = a).
%              Valori ammessi: 'Dirichlet' o 'Neumann'.
%   val_sx   : Scalare. Valore numerico della condizione al contorno a sinistra.
%              Se Dirichlet: u(a) = val_sx. Se Neumann: u'(a) = val_sx.
%   tipo_dx  : Stringa. Tipo di condizione al contorno a destra (x = b).
%              Valori ammessi: 'Dirichlet' o 'Neumann'.
%   val_dx   : Scalare. Valore numerico della condizione al contorno a destra.
%              Se Dirichlet: u(b) = val_dx. Se Neumann: u'(b) = val_dx.
%
% OUTPUT:
%   u        : Vettore colonna (M x 1) contenente la soluzione numerica 
%              approssimata nei nodi della griglia.
%   x        : Vettore colonna (M x 1) dei nodi della griglia spaziale.
%

    % 1. Definizione del dominio nell'intervallo [a, b]
    x = (a:h:b)'; 
    M = length(x); 
    
    A = zeros(M, M);
    b_vettore = zeros(M, 1);

    % 2. Costruzione equazioni per i nodi interni (Universale)
    for i = 2:M-1
        xi = x(i);
        
        % Valutazione delle funzioni nel nodo corrente x_i
        alfa  = alfa_fun(xi);
        beta  = beta_fun(xi);
        gamma = gamma_fun(xi);
        f     = f_fun(xi);
        
        % Assemblaggio della riga i-esima usando le differenze centrate
        A(i, i-1) = -alfa/h^2 - beta/(2*h);       % Peso per u_{i-1}
        A(i, i)   = 2*alfa/h^2 + gamma;           % Peso per u_i
        A(i, i+1) = -alfa/h^2 + beta/(2*h);       % Peso per u_{i+1}
        
        b_vettore(i) = f;                         % Termine noto f(x_i)
    end

    % 3. Condizione al Contorno Sinistra (x = a, indice 1)
    if strcmpi(tipo_sx, 'Dirichlet')
        A(1, 1) = 1;
        b_vettore(1) = val_sx;
    elseif strcmpi(tipo_sx, 'Neumann')
        A(1, 1) = -1/h;
        A(1, 2) = 1/h;
        b_vettore(1) = val_sx;
    else
        error('Tipo di condizione a sinistra non valido.');
    end

    % 4. Condizione al Contorno Destra (x = b, indice M)
    if strcmpi(tipo_dx, 'Dirichlet')
        A(M, M) = 1;
        b_vettore(M) = val_dx;
    elseif strcmpi(tipo_dx, 'Neumann')
        A(M, M-1) = -1/h;
        A(M, M)   = 1/h;
        b_vettore(M) = val_dx;
    else
        error('Tipo di condizione a destra non valido.');
    end

    % 5. Risoluzione
    u = A \ b_vettore;
end