%-------------------------------------------------------------%
%                                                             %
%  This file defines the boundary conditions                  %
%                                                             %
%  It is organised as follows:                                %
%                                                             %
%    - First line:   visc, dens, conduc, capac                %
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
  1   FLUID    15.11e-6  1.00  21.28e-6 1.0   
%  1   FLUID    0.001685  1.00  0.00237 1.0  
    
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
%                        U     V     W     T      Kin     Eps    v_2     f22     TT
  2
  1       WALL          0.0   0.0   0.0  15.0    0.0    0.0   0.0   0.00   0.0
  2       WALL          0.0   0.0   0.0   5.0    0.0    0.0   0.0   0.00   0.0


%==========%
% SYMMETRY %
%==========%

%-- initial conditions:  U, V, W, T, Kin, Eps, v_2, f22, TT
%-- 1 - LES , 2 - RANS
  1
  1   0.0   0.0   0.1   10.0   0.01  0.001  0.05  0.01  0.0
