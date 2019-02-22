%------------------------------------------------------------------------------
%SIMPLS Basic SIMPLS.  Performs no error checking.
function [Xloadings,Yloadings,Xscores,Yscores,Weights] = simpls(X0,Y0,ncomp)

[n,dx] = size(X0);
dy = size(Y0,2);

% Preallocate outputs
outClass = superiorfloat(X0,Y0);
Xloadings = zeros(dx,ncomp,outClass);
Yloadings = zeros(dy,ncomp,outClass);
if nargout > 2
    Xscores = zeros(n,ncomp,outClass);
    Yscores = zeros(n,ncomp,outClass);
    if nargout > 4
        Weights = zeros(dx,ncomp,outClass);
    end
end

% An orthonormal basis for the span of the X loadings, to make the successive
% deflation X0'*Y0 simple - each new basis vector can be removed from Cov
% separately.
V = zeros(dx,ncomp);

Cov = X0'*Y0;
for i = 1:ncomp
    % Find unit length ti=X0*ri and ui=Y0*ci whose covariance, ri'*X0'*Y0*ci, is
    % jointly maximized, subject to ti'*tj=0 for j=1:(i-1).
    [ri,si,ci] = svd(Cov,'econ'); ri = ri(:,1); ci = ci(:,1); si = si(1);
    ti = X0*ri;
    normti = norm(ti); ti = ti ./ normti; % ti'*ti == 1
    Xloadings(:,i) = X0'*ti;
    
    qi = si*ci/normti; % = Y0'*ti
    Yloadings(:,i) = qi;
    
    if nargout > 2
        Xscores(:,i) = ti;
        Yscores(:,i) = Y0*qi; % = Y0*(Y0'*ti), and proportional to Y0*ci
        if nargout > 4
            Weights(:,i) = ri ./ normti; % rescaled to make ri'*X0'*X0*ri == ti'*ti == 1
        end
    end

    % Update the orthonormal basis with modified Gram Schmidt (more stable),
    % repeated twice (ditto).
    vi = Xloadings(:,i);
    for repeat = 1:2
        for j = 1:i-1
            vj = V(:,j);
            vi = vi - (vi'*vj)*vj;
        end
    end
    vi = vi ./ norm(vi);
    V(:,i) = vi;

    % Deflate Cov, i.e. project onto the ortho-complement of the X loadings.
    % First remove projections along the current basis vector, then remove any
    % component along previous basis vectors that's crept in as noise from
    % previous deflations.
    Cov = Cov - vi*(vi'*Cov);
    Vi = V(:,1:i);
    Cov = Cov - Vi*(Vi'*Cov);
end

if nargout > 2
    % By convention, orthogonalize the Y scores w.r.t. the preceding Xscores,
    % i.e. XSCORES'*YSCORES will be lower triangular.  This gives, in effect, only
    % the "new" contribution to the Y scores for each PLS component.  It is also
    % consistent with the PLS-1/PLS-2 algorithms, where the Y scores are computed
    % as linear combinations of a successively-deflated Y0.  Use modified
    % Gram-Schmidt, repeated twice.
    for i = 1:ncomp
        ui = Yscores(:,i);
        for repeat = 1:2
            for j = 1:i-1
                tj = Xscores(:,j);
                ui = ui - (ui'*tj)*tj;
            end
        end
        Yscores(:,i) = ui;
    end
end