 function [E,V,A,R]=mssa(X, M,method,ne,nr) 
%  Syntax: [E,V,A,R]=mssa(X, M); [E,V,A,R]=mssa(X, M, 'BK',.9,1);
%  MSSA performs a Multichannel Singular Spectrum Analysis of the data
%  matrix X, for embedding dimension M.
%
% Input:  X - Data matrix, organized such that each of the L columns of
%             X is a time series of length N.
%         M - Embedding dimension M.
%    method - (Optional) The method of calculating the covariance matrix.
%             Can be 'unbiased', 'biased', 'VG', or 'BK'.
%             The default is 'unbiased', or N-k weighting. 'Unbiased' and
%             'VG' (Vautard/Ghil) are the same. N-weighted ('biased') and
%             Broomhead/King ('BK') are also supported. See AC, BK and COVAR.
%        ne - a fourth (optional) argument can be supplied, which determines 
%             the number of EOFs and PCs to return. If n is less than 1, it 
%             is interpreted as a fractional variance (e.g. n=.9), and enough
%             EOFs and PCs are returned to account for n*100% of the variance.
%             If N is an integer, that number of EOFs and PCs will be returned.
%             The default is to return all EOFs and PCs.
%        nr - A fifth (optional) argument can be supplied, which determines 
%             the number of RCs to return, in the same way as ne functions
%             for PCs and EOFs. The default is to not return any RCs.
%          
%             Note that you must supply NE if you supply NR. You can either 
%             supply NE alone, or both NE and NR, but not NR alone.
%
%  Output:  E - eigenfunction (T-EOF) matrix 
%           V - vector containing variances (unnormalized eigenvalues)
%           A - Matrix of principal components
%           R - Matrix of reconstructed components
%
%  See e. g. Plaut and Vautard, 1994, J. Atm. Sciences 51, 210-236.
%
%  Some general advice on using MSSA:
%    MSSA sometimes seems like a way to turn a small data set into a huge
%    morass of numbers. If you take a small 100 by 100 matrix X and do an
%    MSSA with lag 10, you get: a 1000 by 1000 eigenvector matrix, 1000
%    eigenvalues, a 100 by 1000 matrix of PCs, and a 100 by 100,000 matrix
%    of RCs. To reduce this enormous amount of output, and to make larger
%    problems computationally feasible, the following approaches are possible:
%       
%       1) Rather than using the full data matrix, use a subset of its PCs.
%          These can be computed using EOF or EOFCENT.
%       2) Return EOFs and PCs sparingly if the number of channels
%          and/or the lag are large.
%       3) Don't calculate RCs unless you know you want them. Note that after
%          you have run MSSA, you can still run MRC to calculate particular
%          RCs from the EOFs and PCs.     
%
%  Written by Eric Breitenberger.     Version date 1/11/96
%  Please send comments and suggestions to eric@gi.alaska.edu   
%

[N,L]=size(X);
X=X-ones(N,1)*mean(X);  % Center the data. 

if nargin==2         % mssa(X,M) was called
  method='unbiased';
  ne=L*M;
  nr=0;
elseif nargin==3  
  if isstr(method)   % mssa(X,M,method) was called
    ne=L*M;
    nr=0;
  else               % mssa(X,M,ne) was called
    nr=0; 
    ne=method;
    method='unbiased';
  end
elseif nargin==4
  if isstr(method)    % mssa(X,M,method,ne) was called
    nr=0;
  else                % mssa(X,M,ne,nr) was called
    nr=ne;
    ne=method;
    method='unbiased';
  end
end

[E,V]=mssaeig(X,M,method);

[A]=mpc(X,E);

% figure out how many RCs to calculate:
if nr<1 & nr~=0   % if nr is in the form of fractional variance, convert to an index
  var=nr*sum(V);
  i=find(cumsum(V)>=var);
  nr=i(1);
end  

if nr~=0
  [R]=mrc(A(:,1:nr),E(:,1:nr),L);
end

% figure out how many EOFs/PCs to keep:
if ne<1 & ne~=0   % if ne is in the form of fractional variance, convert to an index
  var=ne*sum(V);
  i=find(cumsum(V)>=var);
  ne=i(1);
end  
disp(['Returning ' num2str(ne) ' EOFs and PCs, and '  num2str(nr) ' RCs.'])
E=E(:,1:ne);
A=A(:,1:ne);

