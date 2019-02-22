figure (6)
YY= pdist(XX','correlation');
ZZ= linkage(YY,'average');
[H,T] = dendrogram(ZZ,0,'colorthreshold','default','labels', sampleID);
xticklabel_rotate;