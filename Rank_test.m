% rank_test - Check the contribution of each trascription factor to the total rank in A
% 
% Usage:
% [res,d_rank]=Rank_test(A)
%           
% res contains the number of TF and its rank in each row
% d_rank contains the number of TF and its rank with rank problem 


function [res,d_rank]=rank_test(A)
[N,L]= size(A);
A=randn(size(A)).*A;
nz= 0;  % --> initialize total # of zeros in the augmented system
res=[];
for l= 1:L,
    v= find(A(:,l)==0);
    nz= nz+length(v);
    R=rank(A(v,:));
    res=[res;l R];
end
d_rank=[];
for k=1:1:L
    if(res(k,2)~=(L-1))
        d_rank=[d_rank;res(k,:)];
    end
end
