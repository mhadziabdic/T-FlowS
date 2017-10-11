%-------------------------------------------------------------%
%                                                             %
%  This file defines the boundary conditions                  %
%                                                             %
%  It is organised as follows:                                %
%                                                             %
%    - First line:   visc, dens                               %
%    - Second line:  number of boundary conditions (Nbound)   %
%    - All the remaining Nbound lines:                        %
%        mark  type        U     V     W                      %
%                                                             %
%   <type> represents the type of the boundary conditions.    %
%   The following types can be used:                          % 
%     1 -> solid wall, prescribed temper.                     %
%     2 -> solid adiabatic wall                               %
%     3 -> outflow                                            %
%     4 -> symetry                                            %
%                                                             %
%  Note:                                                      %
%                                                             %  
%    Each line which begins with character different from:    %
%    '0'-'9' or ' ', is a comment line                        % 
%                                                             % 
%-------------------------------------------------------------%
#--------------------------------------------------%
#            VISc       DENc    CONc        Cp 
  1
  1  Fluid   1.567E-5   1.00    2.2e-5      1.0 

  8 
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
%===========%
% WALL FLUX %
%===========%
    1         WALLFLUX       0.0    0.0      0.0     0.1      0.0  0.0004   0.0  0.0

%===========%
% WALL FLUX %
%===========%
    2         WALLFLUX       0.0    0.0      0.0     0.00      0.0  0.0004   0.0 0.0  

%===========%
% WALL FLUX %
%===========%
    3         WALLFLUX       0.0    0.0      0.0     0.00      0.0  0.0004   0.0 0.0  

%===========%
% WALL FLUX %
%===========%
    4         WALLFLUX       0.0    0.0      0.0     0.00      0.0  0.0004   0.0 0.0  

%===========%
% WALL FLUX %
%===========%
    5         WALLFLUX       0.0    0.0      0.0     0.00      0.0  0.0004   0.0 0.0  

%=========%
% OUTFLOW %
%=========%
    6        OUTFLOW      0.0  0.0  0.0  0.0  0.6  0.0822  1.0004  0.00001

%=========%
% OUTFLOW %
%=========%
    7        OUTFLOW      0.0  0.0  0.0  0.0  0.6  0.0822  1.0004  0.00001

%========%
% INFLOW %
%========%
    8         INFLOW       file    DNS_inflow.dat 


%-- initial conditions:  U, V, W, Kin, Eps, f22, vv
   1
   1  8.0  0.0  0.0  20.0  0.01  0.001  0.08  0.01
%   1  8.0  0.0  0.0  20.0  0.008  10.0  0.008  0.0
