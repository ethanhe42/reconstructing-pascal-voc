function [b,ndx,pos] = unique_no_sort(a)

%UNIQUE_NO_SORT Set unique unsorted.
%   B = UNIQUE_NO_SORT(A) for the array A returns a vector of the unique 
%   elements of A in the order that they appear in A, i.e. B is unsorted.  A can 
%   be a cell array of strings.  UNIQUE_NO_SORT does not have an option to 
%   operate on rows like the original UNIQUE function.
%
%   [B,I,J] = UNIQUE_NO_SORT(A) also returns index vectors I and J such
%   that B = A(I) and A = B(J).
%
%   Caitlin Bever
%   MIT: June 5, 2007


if nargin < 1
  error('MATLAB:UNIQUE:NotEnoughInputs', 'Not enough input arguments.');
elseif nargin > 1
  error('MATLAB:UNIQUE:TooManyInputs', 'Too many input arguments.');
end

% Convert A to a vector

a = a(:)';

% Use UNIQUE to find the forward and backward sorted unique sets and their
% indices.

[b_flip,ndx_flip,pos_flip] = unique(fliplr(a));

% Create the new index by flipping around the backward index.  Use the new
% index to create the unique unsorted set.

ndx = sort(numel(a)-ndx_flip+1);
b = a(ndx);

% Initialize the new position vector with the sorted position vector then
% exchange the elements that were sorted.

if nargout>2

    [b_noflip,ndx_noflip,pos_noflip] = unique(a);    
    pos = pos_noflip;
    
    [y,ind] = sort(b);

    switched_indeces = ind - [1:numel(b)];

    [b_pos,ndx_pos,pos_pos] = unique(pos);
    
    pos_tmp = b_pos + switched_indeces;
    pos = pos_tmp(pos_pos);
    
end