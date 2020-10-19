%% Project 2

clear MinS Filter type tol P NewIm Counter error T
clc
load wbarb

M               = rand(size(X)) > 0.925;
X_Cor           = X.*M;

MinS            = 2; 
Filter          = 'Daubechies';     % Daubechies/Haar
type            = 8;                % Daubechies Type, must be even and greater then 4
tol             = 10^-7;
P               = 2.2;                % Penalty on the cost function



figure(1)
clf

subplot(1,3,1)
imagesc(X)
colormap gray
axis square

subplot(1,3,2)
imagesc(X_Cor)
colormap gray
axis square

subplot(1,3,3)
colormap gray
axis square


drawnow

[NewIm,Counter] = MyCGImageInterp(X_Cor, MinS, Filter,P, tol,type);
relative_error = (abs(norm(NewIm)-norm(X))/norm(X))*100;


figure(2)
imagesc(NewIm)
colormap gray
axis square
title(['Filter=', Filter, ', P=', num2str(P),', Min Level=', num2str(MinS),', Itterations=', num2str(Counter),' Relative Error=', num2str(round(relative_error,2)),'%']);

relative_error = (abs(norm(NewIm)-norm(X))/norm(X))*100
