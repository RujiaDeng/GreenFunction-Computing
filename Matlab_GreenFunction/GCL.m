function [GF,Fai,W,A] = GCL(x, y, S, nmax, mmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GreenCylindersLaplace package by Deng Rujia                            %%%
%%% This main function evaluates the Green function G(x,y)                 %%%
%%% INPUT:                                                                 %%%
%%%    x is a matrix of size o x 3 of o points at which G is evaluated     %%%
%%%    y is a matrix of size p x 3 of p points at which G is evaluated     %%%
%%%    S is the structure of the specific domain describing the problem    %%%
%%%    nmax is the first truncation order (nmax >= 0); if omitted, the     %%%
%%%      default value nmax = 2 is used                                    %%%
%%%    mmax is the second truncation order (mmax >=0); if omitted, the     %%%
%%%      default value mmax = 3 is used                                    %%%
%%% OUTPUT:                                                                %%%
%%%    GF is the matrix of size o x p of values of the Green function      %%%
%%%      evaluated at o points x and p points y                            %%%
%%%    Fai is the matrix of size o x p of value of Fundamental Solution    %%%
%%%      evaluated at o points x and p points y                            %%%
%%%    W is the inverse matrix of size 2M x 2M of ((F,G)(F',G'))',         %%%
%%%       ((F,G)(F',G'))', with M = N*(nmax+1)*mmax. since W does not      %%%
%%%       depend on x and y points, it is recommended to compute the       %%%
%%%       matrix W only once and then supply it for next evaluations of    %%%
%%%	     the Green function in the same domain.                            %%%
%%%    A is the matrix of size 2M x p of (C,D)', i.e. the coefficients     %%%
%%%      C_{mn}^i,D_{mn}^i evaluated at p source points y;                 %%%
%%%      this matrix can be used for next valuations of the Green function %%%
%%%      in the same domain with M = N*(nmax+1)*mmax                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 4)    % default value for the two truncation sizes
    nmax = 2;
    mmax = 3;
end
if (nargin < 5)    % default value for the second truncation size
    mmax = 3;
end

a = 10000;    %set the radius large enough as 10000

% Computing W==((F,G)(F',G'))
N = length(S.h);   % number of cylinders
nn = (nmax+1)*mmax;    % number of base functions per cylinder
M = N*nn;          % dim of the matrix F,F1,G,G1
RB = GCL_rootBessel(); %root of Bessel function
F = zeros(M,M);
F1= zeros(M,M);
G = zeros(M,M);
G1= zeros(M,M);
W = zeros(2*M,2*M);
for j=1:N
    jj = (j-1)*nn;   % index shift in the matrix
    for i=1:N
         ii = (i-1)*nn;  % index shift in the matrix   
         if not(i == j)
             for n=0:nmax
                 for g=0:nmax     
                     for m=1:mmax
                         for h=1:mmax
                             F(ii+(g+1)*h,jj+(n+1)*m)  =  GCL_FG(S,a,i,j,g,h,n,m,0,RB);
                             G(ii+(g+1)*h,jj+(n+1)*m)  =  GCL_FG(S,a,i,j,g,h,n,m,1,RB);
                             F1(ii+(g+1)*h,jj+(n+1)*m) =  GCL_FG(S,a,i,j,g,h,n,m,2,RB);
                             G1(ii+(g+1)*h,jj+(n+1)*m) =  GCL_FG(S,a,i,j,g,h,n,m,3,RB);
                         end
                     end
                 end
             end
         end
     end
end
F1 = eye(size(F1)) + F1;
G = eye(size(G)) + G;
for i=1:M
    for j=1:M
        W(i,j)=F(i,j);
        W(i,j+M)=G(i,j);
        W(i+M,j)=F1(i,j);
        W(i+M,j+M)=G1(i,j);
    end
end

% Similarly, Computing Y=(H,H')'
o = size(x,1); 
p = size(y,1);                % number of starting points          
H = zeros(M,p);               % depends on the starting point y
H1 = zeros(M,p); 
Y = zeros(2*M,p);
for q=1:p
    yq = (y(q,:))';
    for j=1:N
        jj = (j-1)*nn;  % index shift in the matrix  
        for n=0:nmax
            for m=1:mmax
                H(jj+(n+1)*m,q) = GCL_H(S,yq,a,j,n,m,0,RB);   % Computing H
                H1(jj+(n+1)*m,q) = GCL_H(S,yq,a,j,n,m,1,RB);   % Computing H'
            end
        end
    end
end
for i=1:M
    Y(i)=H(i);
    Y(i+M)=H1(i);
end
%lemda=1e-5;
%Computing A=(C,D)'
A=W\Y; %Inv(W)*Y ,of size 2M x p
% +lemda.*diag(diag(ones(W)))
%Computing the patial solution gp
xn = S.x;
gp = zeros(N, o, p);  %size N x o x p
for i=1:N
    ii = (i-1)*nn;
    for q=1:o
        [rouq,faiq,sq] = GCL_cart2cyl(S,i,x(q,:) - xn(i,:));
        for n=0:nmax
            for m = 1:mmax
                gp(i,q,:) = gp(i,q,:) + ( A(ii+(n+1)*m, :).*sin(n*faiq) + A(ii+(n+1)*m+M, :).*cos(n*faiq) ).*(exp(-RB(n+1,m)*sq)*besselj(n,RB(n+1,m)*rouq));
            end
        end
    end
end

% Computing fundamental solution Fai of size o x p
Fai = zeros(o,p);
for i=1:o
  for j=1:p
    Fai(i,j) = 1/(4*pi) / (sqrt( sum( (x(i,:) - y(j,:) ).^2 ) ) );
  end
end

GF = zeros(o,p);
%Computing Green function of size o x p
for i=1:N            % we subtract partial solutions 
    for q=1:o
        for j=1:p
            GF(q,j) = Fai(q,j) - gp(i,q,j);
        end
    end
end