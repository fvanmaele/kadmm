% M = MOPfilter2(n,m,F,boundaryHandling,boundary)
% M = MOPfilter2(n,m,F,boundaryHandling,i,j)
%
% Creates a matrix operator M which applies an FIR-filter
% to the data u:
%   y = M*u
% u and y represent the two-dimensional input and output data, respectively,
% by stacking together the columns of their matrix representation.
% The function expects the data size n and m of u (n rows, m columns)
% and the two-dimensional (FIR) filter mask F (in matrix form) and
% returns a sparse matrix M.
% By default, the filter is applied to all positions within the n x m-data
% where we have enough information, i.e. no extrapolation is required.
% If a 2x2-vector 'boundary' is provided, the filter is applied also to
% additional positions and implements extrapolation of missing values.
% One of the following methods can be chosen for 'boundaryHandling':
%  'zero'                   all outlying are assumed to be zero
%  'linearExtrapolation'    the value is lineary extrapolated
%  'constantExtrapolation'  the value is constant extrapolated (extended)
%  'mirror'                 the data is mirrored at the boundary
% boundary(1,1) and boundary(1,2) specify the number of rows where the
% filter is applied additionally at the beginning and at the end of
% each column, respectively.
% boundary(2,1) and boundary(2,2) specify the number of columns where the
% filter is applied additionally at the beginning and at the end of
% each row, respectively
% Then the y represents the vector representation of an 2D-array of size
% (n - size(F,1) + 1 + boundary(1,1) + boundary(1,2)) times
% (m - size(F,2) + 1 + boundary(2,1) + boundary(2,2)) .
%
% Instead of 'boundary', two equally-sized arrays, 'i' and 'j' can be
% specified which describe the positions where the filter mask shall be
% applied, so y(k) is the result of applying the filter mask at row i(k),
% column j(k) in the data represented by u.
%
% example 1:
%  % define a Gaussian smoothing filter, 3x3
%  F = [0.0625    0.1250    0.0625;...
%       0.1250    0.2500    0.1250;...
%       0.0625    0.1250    0.0625];
%  M = MOPfilter2(10,15,F,'linearExtrapolation',[1,1;1,1]);
%
% example 2:
%  % define a Gaussian smoothing filter, 3x3
%  F = [0.0625    0.1250    0.0625;...
%       0.1250    0.2500    0.1250;...
%       0.0625    0.1250    0.0625];
%  % apply to a sub-grid
%  [i,j] = ndgrid(1:2:10,1:2:15);
%  M = MOPfilter2(10,15,F,'linearExtrapolation',i(:),j(:));

%  File:     MOPfilter2.m
%  Author:   Florian Becker, University of Heidelberg
%  Date:     2012-09-03
%  Language: Matlab
%  Version:  1.2
%  Synopsis: Creates a matrix which acts as an 2D FIR filter operator on 2D
%            data represented as a vector.
%  Description:
%    see above

function M = MOPfilter2(n,m,F,boundaryHandling,varargin)
  
  %% default parameter check
  if (nargin < 5)
    boundary = [0;0];
  end
  
  if (nargin == 5)
    boundary = varargin{1};
  end
  
  if (nargin == 6)
    autoGrid = false;
  else
    autoGrid = true;
  end
  
  if ( (nargin < 4) || isequal(boundaryHandling,[]) )
    boundaryHandling = 'linearExtrapolation';
  end
  
  %% create grid
  [Fn,Fm] = size(F);
  
  if (autoGrid)
    % automatically create grid
    
    B = boundary;
    if size(B,1) == 1
      B = [B;B];
    end
    
    if size(B,2) == 1
      B = [B,B];
    end
    
    if ~isequal(size(B),[2,2])
      error('boundary must be 1x1, 1x2, 2x1 or 2x2 in size');
    end
    
    I = (1 - B(1,1)) : (n - Fn + 1 + B(1,2));
    J = (1 - B(2,1)) : (m - Fm + 1 + B(2,2));
    [I,J] = ndgrid(I,J);
  else
    % use the given parameters as grid
    I = varargin{1};
    J = varargin{2};
  end
  I = I(:);
  J = J(:);
  
  N = numel(I);
  
  %% translate filter mask to a sparse representation (i,j,values)
  [Fi,Fj,Fv] = find(F);
  Fnz = length(Fi);
  Fi = Fi - 1;
  Fj = Fj - 1;
  
  % make sure, they are row vectors
  Fi = Fi(:);
  Fj = Fj(:);
  Fv = Fv(:);
  
  %% apply it for every position in the 2d data grid
  % row in the matrix representation of the data
  I = repmat(I,[1,Fnz]) + repmat(Fi',[N,1]);
  % column in the matrix representation of the data
  J = repmat(J,[1,Fnz]) + repmat(Fj',[N,1]);
  % index in the vector representation of the data
  K = repmat((1:(N))',[1,Fnz]);
  % filter value
  V = repmat(Fv',[N,1]);
  
  % vectorize
  I = I(:);
  J = J(:);
  K = K(:);
  V = V(:);
  
  %% boundary handling
  switch boundaryHandling
    case 'linearExtrapolation'
      % extrapolate in i
      [I,J,K,V] = linearExtrapolateCoordinates(I,J,K,V,n);
      % extrapolate in j
      [J,I,K,V] = linearExtrapolateCoordinates(J,I,K,V,m);
      
    case 'constantExtrapolation'
      % extrapolate in i
      [I,J,K,V] = constantExtrapolateCoordinates(I,J,K,V,n);
      % extrapolate in j
      [J,I,K,V] = constantExtrapolateCoordinates(J,I,K,V,m);
      
    case 'mirror'
      % extrapolate in i
      [I,J,K,V] = mirrorCoordinates(I,J,K,V,n);
      % extrapolate in j
      [J,I,K,V] = mirrorCoordinates(J,I,K,V,m);
      
    case 'zero'
      % extrapolate in i
      [I,J,K,V] = zeroPadding(I,J,K,V,n);
      % extrapolate in j
      [J,I,K,V] = zeroPadding(J,I,K,V,m);
      
    otherwise
      error('invalid value for boundaryHandling');
  end
  
  
  %% prepare parameters for "sparse"-function
  
  % map 2d coordinates (i,j) to 1d coordinates (vector
  % representation), which in turn is the j-index of the matrix operator
  J = I + (J-1) * n;
  I = K;
  
  % finally create the sparse matrix
  M = sparse(I,J,V,N,n*m);
  % please note: duplicate indices (i,j) cause the values to be added by
  % 'sparse'
  
end

%% assume values at outlying positions to be zero
function [i,j,k,v] = zeroPadding(i,j,k,v,n)
  % just remove all out-lying entries
  ix_in = find( and((i>=1),(i<=n)) );
  i = i(ix_in);
  j = j(ix_in);
  k = k(ix_in);
  v = v(ix_in);
end

%% mirror data at the boundary
function [i,j,k,v] = mirrorCoordinates(i,j,k,v,n)
  % extrapolate i<1
  ix_out = find(i<1);
  i(ix_out) = 2 - i(ix_out);
  
  % extrapolate i>n
  ix_out = find(i>n);
  i(ix_out) = 2 * n - i(ix_out);
end

%% constant extrapolation at the boundary
function [i,j,k,v] = constantExtrapolateCoordinates(i,j,k,v,n)
  % extrapolate i<1
  i(i<1) = 1;
  
  % extrapolate i>n
  i(i>n) = n;
end

%% linear extrapolation at the boundary
function [i,j,k,v] = linearExtrapolateCoordinates(i,j,k,v,n)
  % idea: identify outlying coefficients, remove them from the lists i,j,k,v
  % and append two new coefficients for the two position next to the boundary
  % which lineary extrapolate the value.
  % Example: v(i) for i < 1 is approximated by
  %   v(1) + (v(2) - v(1))*(i-1) = v(1)*(2-i) + v(2)*(i-1)
  % 'sparse' will later add overlapping coefficients
  
  % extrapolate i<1
  ix_out = find(i<1);
  ix_in  = find(i>=1);
  n_out  = length(ix_out);
  
  % remove outlying coefficients from the vectors and append coefficients to
  % lineary extrapolate the removed values from
  v = [v(ix_in); v(ix_out).*(2 - i(ix_out)); v(ix_out).*(i(ix_out) - 1)];
  i = [i(ix_in); ones(n_out,1)*1; ones(n_out,1)*2];
  j = [j(ix_in); j(ix_out); j(ix_out)];
  k = [k(ix_in); k(ix_out); k(ix_out)];
  
  % extrapolate i>n
  ix_out = find(i>n);
  ix_in  = find(i<=n);
  n_out  = length(ix_out);
  
  % same as above but on the upper bound
  v = [v(ix_in); v(ix_out) .* ((1-n) + i(ix_out)); v(ix_out).*(n - i(ix_out))];
  i = [i(ix_in); ones(n_out,1)*(n); ones(n_out,1)*(n-1)];
  j = [j(ix_in); j(ix_out); j(ix_out)];
  k = [k(ix_in); k(ix_out); k(ix_out)];
end
