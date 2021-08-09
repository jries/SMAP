function [ MC ] = ELSclique( A )
%{

THIS IS A MATLAB FUNCTION.
AUTHOR: YANHAO WEI
DATE: AUG, 2014
BASED ON: maximalCliques() by Jeffrey Wildman, 2011

   It finds maximal cliques using the Bron-Kerbosch algorithm with
   both pivoting and degeneracy ordering. Degeneracy ordering speeds up the
   algorithm especially when the graph is large & sparse.

   Given a graph's boolean adjacency matrix, A, find all maximal cliques 
   on A using the Bron-Kerbosch algorithm in a recursive manner.  The 
   graph is required to be undirected and must contain no self-edges. A can
   be full or sparse, or even logical.

   The output is a sparse matrix where each column indicates a clique. You
   can change that into a full matrix by using the full() Matlab function.

   Part of the code is based on maximalCliques() by Jeffrey Wildman, 2011.

   Algorithm Ref: Eppstein, Loffler, and Strash "Listing All...", 2010

%}


% first, some input checking

if size(A,1) ~= size(A,2)
    error('MATLAB:maximalCliques', 'Adjacency matrix is not square.');
elseif ~all(all((A==1) | (A==0)))
    error('MATLAB:maximalCliques', 'Adjacency matrix is not boolean (zero-one valued).')
elseif ~all(all(A==A.'))
    error('MATLAB:maximalCliques', 'Adjacency matrix is not undirected (symmetric).')
elseif trace(abs(A)) ~= 0
    error('MATLAB:maximalCliques', 'Adjacency matrix contains self-edges (check your diagonal).');
end


% second, set up some variables

n = size(A,2);      % number of vertices
SC = nan(1e4,2);    % storage of the sparse indices for the matrix of maximal cliques, i.e. output
nc = 0;             % number of found cliques
ns = 0;             % used length of SC
nSC = 1e4;          % preallocated length of SC

% step 2.5: degeneracy ordering

O = nan(1,n);       % to store a degeneracy ordering of the vertices
A0 = double(A);
for i=1:n
    Nb = sum(A0,2);
    [~,j] = min(Nb);
    A0(j,:) = 0; A0(:,j) = 0; A0(j,j) = inf;
    O(i) = j;
end

% third, run the algorithm!

for i = 1:n                     % the outer layer of the algorithm by Eppstein et al. (2010)
    v = O(i);
    R = false(1,n); R(v) = 1;
    Nv = A(v,:)==1;
    P = Nv; P(O(1:i)) = 0;
    X = Nv; X(O(i:end)) = 0;
    BKpivot(R, P, X);           % the inner layer is the Bron-Kerbosch algorithm with pivoting
end

SC = SC(1:ns,:);
MC = sparse(SC(:,1),SC(:,2),1,n,nc);

    function [] = BKpivot ( R, P, X )

        if ~any(P) && ~any(X)                   % report R as a maximal clique
            nc = nc + 1;
            nr = sum(R);
            if ns + nr > nSC
                SC = [SC; nan(1e4,2)];
                nSC = nSC + 1e4;
            end
            SC(ns+1:ns+nr,1) = find(R');
            SC(ns+1:ns+nr,2) = nc;
            ns = ns + nr;
        else
            ppivots = find(or(P,X));            % potential pivots          
            pcounts = A(ppivots,:)*double(P)';  % cardinalities of the sets of neighbors of each ppivots intersected with P
            [~,ind] = max(pcounts);
            u_p = ppivots(ind);                 % select one of the ppivots with the largest count
            
            for u = find(and(A(u_p,:)==0,P))    % all prospective nodes who are not neighbors of the pivot
                Rnew = R; Rnew(u)=1;
                Nu = A(u,:)==1;
                Pnew = and(P,Nu);
                Xnew = and(X,Nu);
                BKpivot(Rnew, Pnew, Xnew);
                P(u) = 0;
                X(u) = 1;
            end
        end
        
    end % BKpivot
       
end % ELSclique