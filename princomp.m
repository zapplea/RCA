function [V, lambda] = princomp(X,mode);
% [V, lambda] = princomp(X, mode);
% computes principal components of
% the data in X (N rows, m cols)
% returns a matrix of eigenvectors V
% and a vector of eigenvalues, lambda.  
% The cols of V correspond to the rows 
% of lambda, which are order by size
% mode = 0 (default) uses the correlation
% matrix, mode = 1 uses the covariance.
% the program plots C1 vs C2, C2 vs C3, 
% and a 3D plot of C1, C2, and C3
% Written by Gerry Middleton, Feb 1997.
if nargin < 2, mode = 0; end
[N,m] = size(X);
if mode == 0
   R = corrcoef(X); 
else 
   R = cov(X);
end; 
[V,D] = eig(R);
lambda = diag(D);  % extract eigenvalues
[lambda,i] = sort(-lambda);  % sort in decreasing order
lambda = -lambda;
V = V(:,i);        % and sort the eigenvectors
Xd = X - ones(N,1)*mean(X); % deviations from mean
if mode == 0
   Z = Xd ./ (ones(N,1)*std(X)); % z values
else
   Z = Xd;
end;
c1 = Z*V(:,1);   % compute the first two components
c2 = Z*V(:,2);
c3 = Z*V(:,3);
plot(c1,c2,'o');  % normally use this
axis('equal');
% plot(c1(1:50),c2(1:50), 'o',c1(51:100),c2(51:100),'x',...
%   'LineWidth',1);  % use this only for iris data
xlabel('C1 scores'), ylabel('C2 scores');
figure;
plot(c2,c3,'o');
axis('equal');
xlabel('C2'), ylabel('C3');
figure;
plot3(c1,c2,c3,'o');
xlabel('C1'), ylabel('C2'),zlabel('C3')