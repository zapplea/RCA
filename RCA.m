function returnVal = RCA(arraydataFileName, TFbinddataFileName, zMode, nmfFlag, keepOriTFFlag)
% this program is Li Huai wrote before he went to China
% (1) Open TF binding data file

fid = fopen(TFbinddataFileName,'r');            
Block = 1;
while (~feof(fid)) 
    tline = fgetl(fid);
    expr = '([\w\-\/\.\$\,\(\)\:\+]*)';
    if (isempty(strfind(tline, '_UNKNOWN')))
        [tokens{Block,1}] = regexp(tline, expr, 'tokens');
        Block = Block+1;
    end  
end
fclose(fid);

TFSym = cell(length(tokens),1);
for i = 1:length(tokens)
     TFSym(i,1) = tokens{i,1}{1,1};
end

% (2) Load original data and gene Symbol ID files

arrayDataStr = importdata(arraydataFileName);

% Get data
subDataX = arrayDataStr.data;

[N, M] = size(subDataX);  % N: number of genes. 
                        % M: number of experiments.
                        
                        

% Get gene symbol of array data
geneSymID = cell(N,1);
RowID = cell(N,1);
for i = 1:N
     geneSymID{i,1} = arrayDataStr.textdata{i+1,2};
     RowID{i,1} = arrayDataStr.textdata{i+1,1};
end
geneSymID = upper(geneSymID);

% Get sample label of array data
sampleID = cell(1,M);
for i = 1:M
     sampleID{1,i} = arrayDataStr.textdata{1,i+2};
end
sampleID = upper(sampleID);

% construct connectionTF sparse matrix from conserved motif information
numTF = length(TFSym);
connectionTF = [];

for ii = 1:numTF
      tmpTFTargetSym = cell(length(tokens{ii,1})-1,1);
      for jj = 1:length(tokens{ii,1})-1
          tmpTFTargetSym(jj) = strtok(tokens{ii,1}{1,jj+1});
      end
      tmpTFTargetSym = upper(tmpTFTargetSym);     
      tf = ismember(geneSymID, tmpTFTargetSym); 
      connectionTF = [connectionTF tf];
      clear tf;
      clear tmpTFTargetSym;
end


% get TF subdataset from subdataset subDataX
subTFDataX = [];
if strcmpi(zMode, 'pls') == 1
    xmean = mean(connectionTF);
    xstd = std(connectionTF);
    ymean = mean(subDataX);
    ystd = std(subDataX);
    connectionTFcentered = (connectionTF - xmean(ones(N,1),:))./xstd(ones(N,1),:);
    subDataXcentered = (subDataX - ymean(ones(N,1),:))./ystd(ones(N,1),:);
    [Xloadings,Yloadings,Xscores,Yscores,Weights] = simpls(connectionTFcentered,subDataXcentered,numTF);
    subTFDataX = Xloadings*Yloadings';
else
    for iiii = 1:numTF
        inx_targetgene = find(connectionTF(:,iiii)==1);
        if isempty(inx_targetgene),
            oneRowTFDataX = subDataX(1,:); %randomly assign
        elseif length(inx_targetgene) == 1,
            oneRowTFDataX = subDataX(inx_targetgene,:);
        else
            switch lower(zMode)
                case 'median'
                    oneRowTFDataX = median(subDataX(inx_targetgene,:));
                case 'kmean'
                    [kmeansInx,Centriod] = kmeans(subDataX(inx_targetgene,:),2, 'distance','correlation','replicates', 10);
                    if std(Centriod(1,:)) > std(Centriod(2,:))
                        oneRowTFDataX = Centriod(1,:);
                    else
                        oneRowTFDataX = Centriod(2,:);
                    end
                case 'pca'
                    [U,S,V] = svd(subDataX(inx_targetgene,:));
                    X2 = U*S;
                    oneRowTFDataX = X2(1,:);
            end
        end
        subTFDataX = [subTFDataX; oneRowTFDataX];
    end
end

% output Z file
str_file_name1111 = [strtok(arraydataFileName,'\.') '_' strtok(TFbinddataFileName,'\.') '_' zMode '_' 'ZMat_' num2str(numTF) '.txt'];
fid1111 = fopen(str_file_name1111, 'wt');
fprintf(fid1111,'%s\t','TF_NAME');
for ii = 1:M-1,
    fprintf(fid1111,'%s\t',sampleID{1,ii});
end
fprintf(fid1111,'%s\n',sampleID{1,M}); %finish firsr row

for i = 1:numTF,
    fprintf(fid1111,'%s\t',TFSym{i});
    for j = 1:M-1,
        fprintf(fid1111,'%f\t',subTFDataX(i,j));
    end
    fprintf(fid1111,'%f\n',subTFDataX(i,M)); %finish one inter row        
end
fclose(fid1111);

% rank standard derivation of TF profile across conditions and select top
% TFs with higher std
stdsubTFDataX = std(subTFDataX, 0, 2);
[sortedStdsubTFDataX,sortedInxsubTFDataX] = sort(stdsubTFDataX,'descend');
if numTF < 100,
    numTopTFs = numTF;
else 
    numTopTFs = 100;
end
sortedsubTFYSYM = TFSym(sortedInxsubTFDataX(1:numTopTFs));
sortedsubTFDataX = subTFDataX(sortedInxsubTFDataX(1:numTopTFs),:);
sortedconnectionTF = connectionTF(:,sortedInxsubTFDataX(1:numTopTFs));
% sortedsubTFYSYM = TFSym(1:numTopTFs);
% sortedsubTFDataX = subTFDataX(1:numTopTFs,:);
% sortedconnectionTF = connectionTF(:,1:numTopTFs);

% (3) Using SCA to infer regulatory network 

 Z0 = sortedsubTFDataX;
 L = numTopTFs;
 sW = 0.86; %0.967; %0.96; %define sparseness
 iteration = 1000;
 fname = ['results_SCA_' num2str(L) '.mat']; 
 %[Y,Z,objhist] = nmfscfixedH( subDataX, Z, L, sW, fname, iteration);
 [Y,Z,objhist] = nmfscfixedHInitW( subDataX, sortedconnectionTF, Z0, L, sW, nmfFlag, keepOriTFFlag, fname, iteration);
 
 
 % plot(objhist);
 
% construct clusters based on the Y  matrix nonzero elements

[gene_inx,tf_inx_repeat,val_weight] = find(Y);

MCell = cell(L,1); % MCell is a cell array in which each element contains a subcluster
clusterInx = [];
clusterInxStart = 0;
for ii = 1:L
   temp_gene_inx_array = [];
   for jj = 1:length(gene_inx)
       if tf_inx_repeat(jj) == ii
           temp_gene_inx_array = [temp_gene_inx_array gene_inx(jj)];
           clusterInx = [clusterInx;clusterInxStart];
       end
   end
   if length(temp_gene_inx_array) ~= 0
       clusterInxStart = clusterInxStart + 1;
   end
   MCell{ii} = temp_gene_inx_array;
   clear temp_gene_inx_array;
end


geneSYM = geneSymID(gene_inx);


% output cluster file
str_file_name3 = [strtok(arraydataFileName,'\.') '_' strtok(TFbinddataFileName,'\.') '_' zMode '_' 'NNTopology_' num2str(L) '.txt']; 
fid3 = fopen(str_file_name3, 'wt');
fprintf(fid3,'%s\t%s\t%s\t%s\t%s\n','Source','Target','IntType','TrueNNVeri','Weight');
for i = 1:length(geneSYM),
    if Y(gene_inx(i),tf_inx_repeat(i))>0
        fprintf(fid3,'%s\t%s\t%s\t%d\t%f\n',sortedsubTFYSYM{tf_inx_repeat(i)},geneSYM{i},'Activate', sortedconnectionTF(gene_inx(i),tf_inx_repeat(i)), Y(gene_inx(i),tf_inx_repeat(i)));
    elseif Y(gene_inx(i),tf_inx_repeat(i))< 0
        fprintf(fid3,'%s\t%s\t%s\t%d\t%f\n',sortedsubTFYSYM{tf_inx_repeat(i)},geneSYM{i},'Inhabit', sortedconnectionTF(gene_inx(i),tf_inx_repeat(i)), Y(gene_inx(i),tf_inx_repeat(i)));
    else
        fprintf(fid3,'%s\t%s\t%s\t%d\t%f\n',sortedsubTFYSYM{tf_inx_repeat(i)},geneSYM{i},'Noeffect', sortedconnectionTF(gene_inx(i),tf_inx_repeat(i)), Y(gene_inx(i),tf_inx_repeat(i)));
    end
end
fclose(fid3);

% Load prior total TFs in human
[humanTFRowID, humanTFGeneIDCell]=xlsread('humanTFListNew.xls');
humanTFGeneID = humanTFGeneIDCell(:,2);
[overlapTFListSYM, inxTFListInsubYSYM, inxInhumanTF]= intersect(geneSymID, humanTFGeneID);

% output node attribute file
meansubDataX = mean(subDataX, 2);
typeTFinsubDataX = zeros(size(geneSymID));
typeTFinsubDataX(inxTFListInsubYSYM) = 1;
str_file_name1 = [strtok(arraydataFileName,'\.') '_' strtok(TFbinddataFileName,'\.') '_' zMode '_' 'NodeAtt_' num2str(L) '.txt']; 
fid1 = fopen(str_file_name1, 'wt');
fprintf(fid1,'%s\t%s\t%s\t%s\n','RowID','Symbol','TYPE','AveExVal');
for i = 1:length(geneSymID),
    fprintf(fid1,'%s\t%s\t%d\t%f\n',RowID{i},geneSymID{i},typeTFinsubDataX(i),meansubDataX(i));
end
fclose(fid1);

% output Y file
fprintf('=====================')
str_file_name11 = [strtok(arraydataFileName,'\.') '_' strtok(TFbinddataFileName,'\.') '_' zMode '_' 'TFCluY_' num2str(L) '.dat'];
fprintf(str_file_name11);

fid11 = fopen(str_file_name11, 'wt');
fprintf(fid11,'%s\t%s\t','UNIQID','NAME');
for ii = 1:L-1,
    fprintf(fid11,'%s\t',sortedsubTFYSYM{ii});
end
fprintf(fid11,'%s\n',sortedsubTFYSYM{L}); %finish firsr row

for i = 1:N,
    fprintf(fid11,'%s\t%s\t',RowID{i},geneSymID{i});
    for j = 1:L-1,
        fprintf(fid11,'%f\t',Y(i,j));
    end
    fprintf(fid11,'%f\n',Y(i,L)); %finish one inter row        
end
fclose(fid11);

returnVal = 0;

