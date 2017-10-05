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
%   The following types of boundary conditions can be used:   % 
%     1 -> WALL (solid wall, prescribed temper.)              %
%     2 -> WALLFLUX (solid adiabatic wall)                    %
%     3 -> INFLOW  (U,V,W,T,....  or type FILE and then name  %
%                   of the file with defined inflow profiles. % 
%                   Look at BouLoa.f90 for more info.)        %                                   %
%     4 -> OUTFLOW                                            %
%     6 -> PRESSURE                                           %
%     7 -> CONVECTIVE                                         %
%     8 -> SYMMETRY                                           %
%                                                             %
%  Note:                                                      %
%                                                             %  
%    Each line which begins with character different from:    %
%    '0'-'9' or ' ', is a comment line                        % 
%                                                             % 
%-------------------------------------------------------------%
  1   % number of domains  
  1            FLUID      0.001                 1.00      0.001408               1.0  
% mark number  mat. type  molecular viscosity   density   thermal conductivity   thermal capacity
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
%  U V W (T or q for cases with temperature)  Kin  Eps   v2 f22
%  U V W (T or q for cases with temperature)  uu vv ww uv uw vw Eps (and f22 for EBM) (for RSM) 
   2
   1 WALLFLUX  0.0 0.0 0.0  0.1 0.0 0.3  0.0000  -10.0
   2 SYMMETRY  0.0 0.0 0.0 20.0 0.0 0.3  0.0000  -10.0

%-- initial conditions:  U, V, W, T, Kin, Eps, v2, f22 
  1
  1   10.00 0.0 0.0 20.0 0.03 0.001 0.8 0.0066 10.0
