% RANKCHECK - Checks the rank of the augmented system matrix after pruning for
%       non-zero entries in A.
% 
%  R = RankCheck(A) - Returns the rank of Au= diag{A,A,...,A} after pruning the
%       rows corresponding to the non-zero entries in A.

function R= RankCheck(alpha)

[N,L]= size(alpha);
alpha=randn(size(alpha)).*alpha;

R= 0;   % --> initialize rank of the augmented sparse system
nz= 0;  % --> initialize total # of zeros in the augmented system

for l= 1:L,
    v= find(alpha(:,l)==0);
    nz= nz+length(v);
    R= R+ rank(alpha(v,:));
end

%fprintf(1,'# of columns in alpha= %i',L);
%fprintf(1,'\n# of zeros in alpha= %i',nz);
%fprintf(1,'\nrank of augmented system after pruning = %i',R);

if (nz < L*(L-1))
    fprintf(1,'\n\nWarning: # of zeros in alpha < L*(L-1)\n');
else
    if (R < L*(L-1))
        fprintf(1,'\n\nWarning: rank(Au) < # L*(L-1)\n');
    end
end