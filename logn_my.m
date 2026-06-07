function y = logn_my(x, n)
%
%   y = logn_my(x, n)
%
% LOGN  Logaritmo in base n di x
%   Usa il cambio di base: log_n(x) = log(x)/log(n)
%
%   Input:
%     x - argomento del logaritmo (>0)
%     n - base del logaritmo (>0, ≠1)
%   Output:
%     y - log in base n di x

    % 🔹 controlli rapidi
    if any(x <= 0)
        error('x deve essere > 0');
    end
    if n <= 0 || n == 1
        error('La base n deve essere > 0 e diversa da 1');
    end

    % 🔹 proprietà del cambio di base
    y = log(x) ./ log(n);
end