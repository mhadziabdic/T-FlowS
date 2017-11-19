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
<<<<<<< HEAD:Tests/RANS/Channel_flow_Retau_590/chan.b
  1 
  FLUID FLUID  0.0001  1.00  0.0001408  1.0  
=======
  1
  FLUID   FLUID    0.0001  1.00  0.0001408  1.0  
>>>>>>> 77a0dfd04f33aeb0ac6ea1ace412064246261650:Tests/Rans/Channel_Re_Tau_590/chan.b
    
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
%                        U        V       W    Kin        Eps       vi2 
   2
<<<<<<< HEAD:Tests/RANS/Channel_flow_Retau_590/chan.b
   ONE  WALLFLUX       0.0   0.0   0.0  0.1    0.0 0.3  0.0000  -10.0
   TWO  SYMMETRY       0.0   0.0   0.0 20.0    0.0 0.3  0.0000  -10.0

%-- initial conditions:  U, V, W, Kin, Eps
%-- 1 - LES , 2 - RANS
=======
   ONE   WALLFLUX       0.0   0.0   0.0  0.1    0.0 0.3  0.0000  -10.0
   TWO   SYMMETRY       0.0   0.0   0.0 20.0    0.0 0.3  0.0000  -10.0

%----------------------------------------%
% Initial conditions:  U, V, W, Kin, Eps %
%----------------------------------------%
>>>>>>> 77a0dfd04f33aeb0ac6ea1ace412064246261650:Tests/Rans/Channel_Re_Tau_590/chan.b
  1
  1   1.00  0.0  0.0 20.0  0.01  0.001  0.1  0.0066 10.0
