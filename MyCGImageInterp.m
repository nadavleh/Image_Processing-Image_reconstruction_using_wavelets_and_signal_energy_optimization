function [ xj,Counter ] = MyCGImageInterp( Im, MinS, Filter, P, tol,type )
% This MATLAB function takes an image (Im) with missing data and preforms
% minimization of a cost function, defined by the minimun level of wavelets
% transform (MinS), wavelet's type ('type' = 'Haar' or 'Daubechies' for
% instance) and a penalty (P). The function runs as long as the residual's
% norm is larger then 'tol'.
% 
% This practicly can result in an efficient interpulation if there is a
% large amount of missing data.


if nargin < 5
    tol = 10^-1;
    Filter = 'Haar';
elseif nargin < 6
    Filter = 'Haar';
end

n               = size(Im,1);
qmf             = MakeONFilter(Filter,type);
MaxS            = log2(n);

D = zeros(n);
for s = MaxS:-1:MinS
    D(1:2^s,1:2^s) = 2^(P*s);
end

Const    = (Im~=0); 

xj          = Im;

temp        = FWT2_PO(xj,MinS,qmf);
temp        = temp.*D;
temp        = IWT2_PO(temp,MinS,qmf);

temp(logical(Const)) = Im(logical(Const));
rj          = 0 - temp;
pj          = rj; 
norm_rj     = norm(rj);
Counter     = 0;

while norm_rj > tol
    
    temp        = FWT2_PO(pj,MinS,qmf);
    temp        = temp.*D;
    temp        = IWT2_PO(temp,MinS,qmf);
    
    pj(logical(Const)) = 0;
    
    alpha_jp1 = (norm_rj^2)/(pj(:)'*temp(:)); 
    
    xjp1 = xj + alpha_jp1*pj;
    rjp1 = rj - alpha_jp1*temp;
    
    rjp1(Const) = 0;
    
    norm_rjp1 = norm(rjp1);
    beta_jp1 = (norm_rjp1/norm_rj)^2;
    pjp1 = rjp1 + beta_jp1*pj;             
    
    xj = xjp1;
    rj = rjp1;
    pj = pjp1;
    norm_rj = norm(rj(:));
        
    Counter = Counter + 1;
    
    if ~mod(Counter,150)
        imagesc(xj)
        axis square
        drawnow
        Counter
        norm_rj
    end
    
end

end
