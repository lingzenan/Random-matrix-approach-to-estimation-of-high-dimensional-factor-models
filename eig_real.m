%arg: R : the input data
%        p: the number of factors
%rtn:   the eigenvalue of the covariance matrix of the residuals
function E = eig_real( R,p)
[~,T] = size(R);
[F,~] = eigs(R'*R,p); % obtain the p largest eigenvalues and corresponding eigenvectors
F=F';
L =  R*F'/(F*F');
U = R-L*F;
U = mapstd(U);
C = U*U'/T;
E= eig(C);
end