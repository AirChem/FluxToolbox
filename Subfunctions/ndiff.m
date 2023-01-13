function poo = ndiff(vector,n,dim)
%ndiff takes the difference vector of a vector and tacks on a NaN at the
%end so that yer diff vector is the same size as yer startin vector.
%
% 20211220 GMW Modified to work for all n.

if nargin<3, dim = 1; end
if nargin<2, n = 1; end

poo = diff(vector,n,dim);
poo(length(poo)+n) = NaN;