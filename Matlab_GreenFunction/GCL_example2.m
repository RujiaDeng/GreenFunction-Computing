function [GF] = GCL_example2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                         %%%
%%%  INPUT: one needs to set first the geometric domain stored in system    %%%
%%%  'S': including the positions, height, radii of the top/bottom          %%%
%%%  circle and unit direction vector of N non-overlapping cylinders.       %%%
%%%  Second, one provides an array of points, at which the Green function   %%%
%%%  should be calculated.                                                  %%%
%%%    S.x is a matrix of size N x 3, with i-th row w.r.t the location      %%%
%%%        of the i-th ball                                                 %%%
%%%    S.h is a vector of size N x 1 w.r.t the height of the cylinders      %%%
%%%    S.r is a vector of size N x 1 w.r.t the radius of the cycle          %%%
%%%    S.u is a matrix of size N x 3,with i-th row w.r.t the unit direction %%%
%%%        vector                                                           %%%
%%%    x is a matrix of size ox3 of the 'o' starting points                 %%%
%%%    y  is a matrix of size px3 of the 'p' arrival points                 %%%
%%%                                                                         %%%
%%%  OUTPUT: the matrix o x p of the values of the Green function G(x,y)    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%z%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Starting computation for the example 2..\n');
S.x = [[-2 0 0]; [2 0 0]; [0 4 0]];   % positions of cylinders
S.h = [1 1 0.5]';                     % heights of cylinders
S.r = [0.8 1 0.5]';                  % radius of cylinders
S.u = [[0.6,0,0.8]; [0 1 0]; [0 0 1]];   % unit direction vector of cylinders
%S.x = [[0 0 0]; [4 0 0]; [0 3 0]];   % positions of cylinders
%S.h = [1 1 0.5]';                     % heights of cylinders
%S.r = [1 1 0.5]';                  % radius of cylinders
%S.u = [[0 0 1]; [0 1 0]; [0 1 0]];   % unit direction vector of cylinders
y = [0 0 2];              % starting point
xx = (-10:0.1:10); 
x = [xx; zeros(size(xx)); 1.5*ones(size(xx))]';   % line of arrival points

%%%  We compute the Green Function
[GF,Fai,W,A] = GCL(x, y, S);
%plot(xx,GF);
%legend('Green function');
plot(xx,Fai);
legend('Fundamental solution');