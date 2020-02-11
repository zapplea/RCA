function [W,H,objhistory] = nmfscfixedHInitW( V, W0, H, rdim, sW, nmfFlag, keepOriTFFlag, fname, iteration)
% nmfsc - non-negative matrix factorization with sparseness constraints
%
% SYNTAX:
% [W,H] = nmfsc( V, rdim, sW, sH, fname, showflag );
%
% INPUTS:
% V          - data matrix
% rdim       - number of components (inner dimension of factorization)
% sW         - sparseness of W, in [0,1]. (give [] if no constraint)
% sH         - sparseness of H, in [0,1]. (give [] if no constraint)
% fname      - name of file to write results into
% showflag   - binary flag. if set then graphically show progress
%
% Note: Sparseness is measured on the scale [0,1] where 0 means
% completely distributed and 1 means ultimate sparseness.
%
% NOTE: There is NO CONVERGENCE CRITERION. The estimation never ends,
% but rather has to be terminated manually. See README file of code
% package for details.
%


% Globally rescale data to avoid potential overflow/underflow
V = V/max(V(:));

% Data dimensions
vdim = size(V,1);
samples = size(V,2);

% Create initial matrices
if ~isempty(W0),
    W = W0;
else
    W = randn(vdim,rdim);
end
H = H./(sqrt(sum(H.^2,2))*ones(1,samples));
% Huai compute sparseness of W0 by column
 

    for i=1:rdim, 
        L1 = 0;
        L2 = 0;
        for j=1:vdim, 
            L1 = L1 + abs(W(j,i));
            L2 = L2 + W(j,i)^2;
        end
        sWs(i) = (sqrt(vdim)- L1/sqrt(L2))/(sqrt(vdim)-1);
    end

% % Make initial matrices have correct sparseness
if ~isempty(sW),
    L1a = sqrt(vdim)-(sqrt(vdim)-1)*sW;
%     for i=1:rdim, 
%         if (sWs(i) < sW)
%             W(:,i) = projfunc(W(:,i),L1a,1,nmfFlag); 
%         end
%     end
end

% Calculate initial objective
objhistory = 0.5*sum(sum((V-W*H).^2));

% Initial stepsizes
stepsizeW = 1;

timestarted = clock;

% Start iteration
for iter=1:iteration,
    fprintf('iter: ',iter)
    iter1 = iter-1;
    % Show progress
    fprintf('[%d]: %.5f\n',iter1,objhistory(end));

    % Save every once in a while
    if rem(iter1,5)==0,
        elapsed = etime(clock,timestarted);
        fprintf('Saving...');
        save(fname,'W','H','sW','iter1','objhistory','elapsed');
%        csvwrite('/home/yibing/Documents/code/bioinfo/RCA/RCAresult_Y.txt',W)
%        csvwrite('/home/yibing/Documents/code/bioinfo/RCA/RCAresult_sW.txt',sW)
%        csvwrite('/home/yibing/Documents/code/bioinfo/RCA/RCAresult_H.txt',H)
%        csvwrite('/home/yibing/Documents/code/bioinfo/RCA/RCAresult_iter1.txt',iter1)
%        csvwrite('/home/yibing/Documents/code/bioinfo/RCA/RCAresult_objhistory.txt',objhistory)
%        csvwrite('/home/yibing/Documents/code/bioinfo/RCA/RCAresult_elapsed.txt',elapsed)
        fprintf('Done!\n');
    end
 
    % ----- Update W ---------------------------------------

    if ~isempty(sW),

        % Gradient for W
        dW = (W*H-V)*H';
        begobj = 0.5*sum(sum((V-W*H).^2));

        % Make sure we decrease the objective!
        while 1,

            % Take step in direction of negative gradient, and project
            Wnew = W - stepsizeW*dW;
            norms = sqrt(sum(Wnew.^2));
            for i=1:rdim,
                Wnew(:,i) = projfunc(Wnew(:,i),L1a*norms(i),(norms(i)^2),nmfFlag);
            end
            if (nmfFlag)
                Wnew(find(Wnew<0)) = 0;
            end
            % Calculate new objective
            newobj = 0.5*sum(sum((V-Wnew*H).^2));

            % If the objective decreased, we can continue...
            if newobj<=begobj,
                break;
            end

            % ...else decrease stepsize and try again
            stepsizeW = stepsizeW/2;
            fprintf(',');
            if stepsizeW<1e-200,
                fprintf('Algorithm converged.\n');
                return;
            end

        end

        % Slightly increase the stepsize
        stepsizeW = stepsizeW*1.2;
        
        %Huai modified
        W = Wnew;
        if keepOriTFFlag == 1,
            min_a = 1.0; 
            max_b = 5.0;
            RandW = min_a + (max_b-min_a)*rand(size(W0));     
            % commented by Yibing Liu based on the requirements of Dr. Yan
            % Bin.
            %W(intersect(find(abs(W)<1.0),find(W0~=0)))=RandW(intersect(find(abs(W)<1.0),find(W0~=0)));
            %W(intersect(find(W==0),find(W0~=0)))=W0(intersect(find(W==0),find(W0~=0)));
            %W(intersect(find(W~=0),find(W0==0)))=W0(intersect(find(W~=0),find(W0==0)));
        end
    end

    % Calculate objective
    newobj = 0.5*sum(sum((V-W*H).^2));
    objhistory = [objhistory newobj];
    if newobj < 1.5
        break;
    end
end
huai = 0;
