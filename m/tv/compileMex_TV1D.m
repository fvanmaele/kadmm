%{
   File            : compileMex_TV1D.m
 
   Version History : 1.0, October 11 2016
 
   Author          : Stephen Becker, PhD, University of Colorado Boulder
 
   Description     : This file compiles the mex/C files associated
                      with the algorithms described in the research paper:
 	
                     L. Condat, "A Direct Algorithm for 1D Total Variation
                     Denoising", preprint hal-00675043, 2012.
 
                     This implementation comes with no warranty: due to the
                     limited number of tests performed, there may remain
                     bugs. In case the functions would not do what they are
                     supposed to do, please email the author (contact info
                     to be found on the web).
 
                     If you use this code or parts of it for any purpose,
                     the author asks you to cite the paper above or, in 
                     that event, its published version. Please email him if 
                     the proposed algorithms were useful for one of your 
                     projects, or for any comment or suggestion.
 
   Usage rights    : Copyright Stephen Becker.
                     This file is distributed under the terms of the CeCILL
                     licence (compatible with the GNU GPL), which can be
                     found at the URL "http://www.cecill.info".
 
   This software is governed by the CeCILL license under French law and
   abiding by the rules of distribution of free software. You can  use,
   modify and or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL :
   "http://www.cecill.info".
 
   As a counterpart to the access to the source code and rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.
 
   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.
 
   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
%}

%% Compile the files
mex -v TV1D_denoise_mex.c -largeArrayDims
mex -v TV1D_denoise_tautString_mex.c -largeArrayDims
mex -v fused_lasso_mex.c -largeArrayDims

%% Simple test:
x = -((1:10)-5).^2 + .1*randn(1,10);

y   = TV1D_denoise_mex( x, 3 );
yy  = TV1D_denoise_tautString_mex( x, 3 );
y2  = fused_lasso_mex( x, 3, .3 );
figure(1); clf;
plot(x,'o-');
hold all
plot(y,'*--')
plot(y2,'s-.')
legend('data','TV','fused lasso (TV + l1)');
norm(y-yy)

%% Test what we are solving, compared to CVX solution (denoising)
N   = 10;
rng(324);
x = -((1:N)-5).^2 + .1*randn(1,N); x = x';
D   = spdiags( [[ones(N-1,1);0], -ones(N,1) ], 0:1, N, N );
% doesn't matter if it is D or -D
lambda = 3;
cvx_begin
 variable z(N)
 minimize sum_square(z-x)/2 + lambda*norm( D*z, 1 )
cvx_end
y   = TV1D_denoise_mex( x, lambda );
% [z,y,x]
norm(z-y)
obj = @(z) sum_square(z-x)/2 + lambda*norm( D*z, 1 );
obj(z) - obj(y)
% Good! These match.

%% Next test (fused lasso)
lambda = 3;
mu     = 2;
cvx_begin
 cvx_precision best
 variable z(N)
 minimize sum_square(z-x)/2 + lambda*norm( D*z, 1 ) + mu*norm(z,1)
cvx_end
y2   = fused_lasso_mex( x, lambda, mu );
% [z,y2,x]
norm(z-y2)
obj = @(z) sum_square(z-x)/2 + lambda*norm( D*z, 1 )+ mu*norm(z,1)
obj(z) - obj(y2)

%% Some timing
nList = round(logspace(5,8,5));
TIMES = zeros(3,length(nList));
for k = 1:length(nList)
    N   = nList(k);
    x   = randn(N,1);
    lambda = .2;
    mu     = .3;
    
    t1=tic;
    y   = TV1D_denoise_mex(x,lambda);
    t   = toc(t1);
    TIMES(1,k) = t;
    
    t1=tic;
    y   = TV1D_denoise_tautString_mex(x,lambda);
    t   = toc(t1);
    TIMES(2,k) = t;
    
    t1=tic;
    y   = fused_lasso_mex(x,lambda,mu);
    t   = toc(t1);
    TIMES(3,k) = t;
end
%% Plot the timing results
clf;
loglog( nList, TIMES,'o-' ); hold all
loglog( nList, nList/nList(2)*TIMES(1,2), '--');
legend('TV','TV taut string','fused lasso','linear scaling');
