function Pe = peclet_locale_my(beta, alfa, h)
%
%   Pe = peclet_locale_my(beta, alfa, h)
%
% Calcola il numero di Peclet LOCALE per un problema di
% diffusione-trasporto:  -alfa*u'' + beta*u' + gamma*u = f
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

    Pe = abs(beta) * h / (2 * alfa);

end
