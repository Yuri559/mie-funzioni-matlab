function Pe = peclet_globale_my(eta, mu, L)
%
%   Pe = peclet_globale_my(eta, mu, L)
%
% Calcola il numero di Peclet GLOBALE per un problema di
% diffusione-trasporto:  -mu*u'' + eta*u' + sigma*u = f
%
% INPUT:
%   beta : coefficiente del trasporto (modulo)
%   alfa : coefficiente di diffusione (viscosita')
%   L    : lunghezza dell'intervallo  L = b - a
%
% OUTPUT:
%   Pe   : numero di Peclet globale
%          Pe > 1  ->  trasporto dominante

    Pe = abs(eta) * L / (2 * mu);

end
