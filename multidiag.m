function [A] = multidiag(v,n)
% multiDiag: costruisce matrici multidiagonali (tri- o penta-diagonale)
% Uso:
%   A = multidiag(v,n)
% Input:
%   v : vettore dei coefficienti lungo le diagonali, ordinati simmetricamente
%       - se numel(v)==3: v = [sotto, diag, sopra]  (tridiagonale)
%       - se numel(v)==5: v = [sottosotto, sotto, diag, sopra, soprasopra] (pentadiagonale)
%   n : dimensione della matrice quadrata A (n x n)
% Output:
%   A : matrice n x n con i valori di v posti sulle rispettive diagonali
% Note:
%   I vettori per le diagonali sono replicati con la lunghezza appropriata.
A=zeros(n);
if(numel(v)==3)
    A= A+ diag(v(2)*ones(n,1)) + diag(v(1)*ones(n-1,1),-1) + diag(v(3)*ones(n-1,1),+1);
end
if(numel(v)==5)
        A= A+ diag(v(3)*ones(n,1)) + diag(v(2)*ones(n-1,1),-1) + diag(v(4)*ones(n-1,1),+1) + diag(v(1)*ones(n-2,1), -2) + diag(v(5)*ones(n-2,1),+2);
end
