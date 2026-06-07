function [v_A] = normA_my(v,A)
%
%   [v_A] = normA_my(v,A)
%
%NORMA_MY restituisce la norma A (norma energia) di un vettore rispetto a
%una matrice SDP
%
%   ||v||_A = sqrt(v' * A * v)
%
%INPUT
%   v   →   Vettore
%   A   →   Matrice (Simmetrica e Definita positiva) 
%
%OUTPUT
%   v_A → ||v||_A

[R, p] = chol(A);

if issymmetric(A) && p == 0
    v_A = sqrt(v' * A * v);
else
    fprintf("La matrice A non è SDP\n");
end