%/////////////////////%
% PHYSICAL PROPERTIES %
%/////////////////////%
  1
  1   fLUID  0.01  1.00  

%/////////////////////%
% BOUNDARY CONDITIONS %
%/////////////////////%
   3 
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
   1         wALL         0.0      0.0   0.0
%========%
% INFLOW %
%========%
   2         iNFLOW   file  parabolic-0-3.prof 
%=========%
% OUTFLOW %
%=========%
   3         oUTFLOW  file  parabolic-0-3.prof 

%////////////////////%
% INITIAL CONDITIONS %
%////////////////////%
   1
   1   1.0  0.0  0.0
