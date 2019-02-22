% This code is for normalization of CS(a) and TFA(p) matrix after NCA run
% [s,a,p]=NCA(E,A,0);
% L : the number of transcription factor
% N : the number of genes
% Usage: function [nor_a, nor_p]=normalization(a,p)
% nor_a: normalized a(CS); nor_p: normalized p(TFA)
% written by Young Lyeol Yang


function [nor_a, nor_p]=normalization(a,p)
[N,L]=size(a);
res=[];
NF=[];
for(k=1:1:L)
    for(kk=1:1:N)
        if (a(kk,k)~=0)
            temp=abs(a(kk,k));
            res=[res;temp];
        end
    end
    temp1=mean(res,1);
    NF=[NF;temp1];
    res=[];
end

% The following codes are for calculating the normalization factor of CS columns
NF1=NF.^-1;

% The following codes are for calculating the normalized CS(a) and TFA(p)
for(k=1:1:L)
   nor_a(:,k)=a(:,k)*NF1(k,1);
   nor_p(k,:)=p(k,:)*NF(k,1);
end
   