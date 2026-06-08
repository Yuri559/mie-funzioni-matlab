function Pe = peclet_locale_my(eta, mu, h)
%
%   Pe = peclet_locale_my(eta, mu, h)
%
% Calcola il numero di Peclet LOCALE per un problema di
% diffusione-trasporto:  -mu*u'' + eta*u' + sigma*u = f
%
% INPUT:
%   beta : coefficiente del trasporto (modulo)
%   alfa : coefficiente di diffusione (viscosita')
%   h    : passo di discretizzazione della griglia
%
% OUTPUT:
%   Pe   : numero di Peclet locale
%          Pe < 1  ->  centrate stabili
%          Pe > 1  ->  centrate oscillano: usa upwind o Scharfetter-Gummel

    Pe = abs(eta) * h / (2 * mu);

end
