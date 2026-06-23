function [radiic] = gershcircles(A)
% 
% gershcircles(A)
%
% Rappresenta i cerchi di Gershgorin nel piano complesso associati alla
% matrice A quadrata ed a valori reali
%
% Parametri di ingresso:
% A         Matrice quadrata (n x n) a valori reali
%
%                                         Politecnico di Milano, 04/04/2024
%

% Verifica che A sia quadrata, altrimenti errore
n = size(A);
if n(1) ~= n(2)
    error('Solo matrici quadrate');
else
    n = n(1);
    % Prealloca matrici per punti dei cerchi (201 campioni)
    circler = zeros(n,201);
    circlec = circler;
end
% Calcola centri e raggi Gershgorin per righe e colonne
center = diag(A);
radiic = sum(abs(A-diag(center)));
radiir = sum(abs(A'-diag(center)));
% Vettori campione e angoli (201 punti, ciclo chiuso)
one = ones(1,201);
cosisin = exp(i*[0:pi/100:2*pi]);
% Finestra per cerchi riga
h1 = figure;
title('Cerchi riga');
xlabel('Re'); ylabel('Im');
% Finestra per cerchi colonna
h2 = figure;
title('Cerchi colonna');
xlabel('Re'); ylabel('Im');
% Per ogni riga/colonna costruisci e disegna i cerchi
for k = 1:n
   circlec(k,:) = center(k)*one + radiic(k)*cosisin;
   circler(k,:) = center(k)*one + radiir(k)*cosisin;
   figure(h1);
   patch(real(circler(k,:)),imag(circler(k,:)),...
   'red');
   hold on
   plot(real(circler(k,:)),imag(circler(k,:)),'k-',...
   real(center(k)),imag(center(k)),'kx');
   figure(h2);
   patch(real(circlec(k,:)),imag(circlec(k,:)),...
   'blue');
   hold on
   plot(real(circlec(k,:)),imag(circlec(k,:)),'k-',...
   real(center(k)),imag(center(k)),'kx');
end
% Ridisegna contorni per assicurare tracciati neri sovrapposti
for k = 1:n
   figure(h1);
   plot(real(circler(k,:)),imag(circler(k,:)),'k-',...
   real(center(k)),imag(center(k)),'kx');
   figure(h2);
   plot(real(circlec(k,:)),imag(circlec(k,:)),'k-',...
   real(center(k)),imag(center(k)),'kx');
end
figure(h1);
axis image;
hold off
figure(h2);
axis image;
hold off
% disegno tutti i cerchi nello stesso grafico
h3 = figure;
for k = 1:n
   figure(h3);
   title('Cerchi di Gershgorin')
   xlabel('Re'); ylabel('Im');
   plot(real(circler(k,:)),imag(circler(k,:)),'g','LineWidth', 2);
   hold on
   plot(real(circlec(k,:)),imag(circlec(k,:)),'k-');
   plot(real(center(k)),imag(center(k)),'rx','MarkerSize', 8);
end
legend('cerchi riga', 'cerchi colonna', 'centro cerchi','Location','NorthOutside' );
figure(h3);
axis image;
hold off
return
