function X=grid(P, N)
%
% Produce a regular grid of this set.
% <p>
% Algorithm:
% <ol>
% <li>Compute outer bounding hypercube</li>
% <li>Grid the hypercube</li>
% <li>Test each point for inclusion in the set, discarding those outside</li>
% </ol>
%
% @param N Number of points in each dimension of the grid
% @return Matrix of n grid points in <eq>\mathbb{R}^{n\times d}</eq>
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS=mptopt;
end

narginchk(2, 2);
error(P.rejectArray());

if P.isEmptySet
    error('Empty set, there is nothing to be gridded here.');
end

if ~P.isBounded
    error('Can only grid bounded sets.'); 
end

% get bounding box
bbox = P.outerApprox;

lb = bbox.Internal.lb;
ub = bbox.Internal.ub;


% grid the state-space into equidistantly placed points
dimB = size(lb(:),1);
Xpoints = zeros(N, dimB);
for ii=1:dimB
    Xpoints(:,ii) = linspace(lb(ii),ub(ii),N)';
end

% generate all possible combinations of states
% one could use kron() here, but that one fails for high number of elements
n_states = dimB;
ZZ=[];
ZZ{n_states}=Xpoints(:,n_states);
for ii=n_states-1:-1:1,
    Zd=[];
    for jj=1:size(Xpoints,1),
        Zd=[Zd; repmat(Xpoints(jj,ii),length(ZZ{ii+1}),1)];
    end
    ZZ{ii}=Zd;
end
for ii=2:n_states,
    ZZ{ii}=repmat(ZZ{ii},length(ZZ{ii-1})/length(ZZ{ii}),1);
end
datapoints=[];
for ii=1:n_states,
    datapoints(:,ii) = ZZ{ii};
end
npoints = size(datapoints,1);
isin = false(npoints, 1);
for i = 1:npoints,
    isin(i) = P.contains(datapoints(i,:)');
end
X = datapoints(isin, :);



end
