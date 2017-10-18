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
  1
  1   FLUID    0.00033333  1.00  
    
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
%                        U        V       W    Kin        Eps       vi2 
   1

   1       WALL         0.0   0.0     0.0   0.00 0.3  0.0000  -10.0



%==========%
% SYMMETRY %
%==========%

%-- initial conditions:  U, V, W, Kin, Eps
%-- 1 - LES , 2 - RANS
  1
  1   0.00   0.0     1.0   0.03 0.003  0.008  0.0066
