function consensus = scaconsensus(a,kstart,kend,nloop)
%
% This software and its documentation are copyright 2004 by the
% Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever.
% Neither the Broad Institute nor MIT can not be responsible for its use, misuse,
% or functionality.
%
% Model selection for SCA
%
% a (n,m) : N (genes) x M (samples) original matrix
%           numerical data only.
%
% kstart, kend : range of values of k to test consensus for.
%
% nloop : number of initial conditions per k
%         (start with 10 or 20, check with more)
%
% consensus : 3d array of consensus matrices
%             dimensions : kend x M x M
%             Values in consensus(1:kstart-1,:,:) should be ignored
%

[n,m]=size(a);

consensus=zeros(kend,m,m);
conn=zeros(m,m);

sW = 0.5; %0.8; %define sparseness
sH = []; %0.3; %define sparseness
iteration = 500;
 
for j=kstart:kend
    connac=zeros(m,m);
    [kmeansInx,h0] = kmeans(a, j);
    for iloop=1:nloop;
        %[kmeansInx,h0] = kmeans(a, j);
        [w,h,objhist] = nmfscHL1( a, h0, j, sW, sH, [], iteration);
        %[w,h,objhist] = nmfscHL2( a, h0, j, [], sH, [], iteration);
        %
        % compute consensus matrix
        %
        conn=nmfconnectivity(h);
        connac=connac+conn; % accumulate connectivity matrices
    end
    consensus(j, :, :)=connac/nloop; %average
end

function conn= nmfconnectivity(h)
%
% Jean-Philippe Brunet
% Cancer Genomics 6/10/03
%
mm=size(h);
k=mm(1);
m=mm(2);


% compute m x m matrix which is 1 if samples are together, 0 elsewhere

% determine sample assignment by its largest metagene expresion value
[y,index]=max(h,[],1);

mat1=repmat(index,m,1); % spread index down
mat2=repmat(index',1,m); % spread index right

conn=mat1==mat2; % 1 when for pair of samples with same assignement



