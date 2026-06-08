function Pe = peclet_globale_my(beta, alfa, L)
%
%   Pe = peclet_globale_my(beta, alfa, L)
%
% Calcola il numero di Peclet GLOBALE per un problema di
% diffusione-trasporto:  -alfa*u'' + beta*u' + gamma*u = f
%
% INPUT:
%   beta : coefficiente del trasporto (modulo)
%   alfa : coefficiente di diffusione (viscosita')
%   L    : lunghezza dell'intervallo  L = b - a
%
% OUTPUT:
%   Pe   : numero di Peclet globale
%          Pe > 1  ->  trasporto dominante

    Pe = abs(beta) * L / (2 * alfa);

end
