function [GF] = GCL_example1()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                         %%%
%%%  INPUT: one needs to set first the geometric domain stored in system    %%%
%%%  'S': including the positions, height, rad-overlapping cylinders.       %%%
%%%  Second, one provides an array of points, at which the Green function   %%%
%%%  should be calculated.                    ii of the top/bottom          %%%
%%%  circle and unit direction vector of N non                              %%%
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
fprintf('Starting computation for the example 1...\n');
S.x = [[0 0 0]];  % positions of cylinders
S.h = [1];        % heights of cylinders
S.r = [1];        % radius of cylinders
S.u = [[0 0 1]];  % unit direction vector of cylinders
y = [0 0 2];              % starting point
xx = (-10:0.1:10); 
x = [xx; zeros(size(xx)); 1.5*ones(size(xx))]';   % line of arrival points

%%%  We compute the Green Function
[GF,Fai,W,A] = GCL(x, y, S);
%plot(xx,GF);
%legend('Green function');
plot(xx,Fai);
legend('Fundamental solution');
%for nmax=2:5
 %   [GF] = GCL(x, y, S,nmax);
 %   plot(xx,GF);
  %  hold on;
%end
%legend('nmax=2','nmax=3','nmax=4','nmax=5');   