function r = tom_pointrotate(r,phi,psi,the)
%TOM_POINTROTATE rotates point (= 3d vector)
%
%   r = tom_pointrotate(r,phi,psi,the)
%
%   A vector in 3D is rotated around the origin = [0 0 0]. The puropose is
%   for example to predict the location of a point in a volume after
%   rotating it with tom_rotate3d. Take care that the coordinates are with
%   respect to the origin!
%
%PARAMETERS
%
%  INPUT
%   r                   3D vector - e.g. [1 1 1 ]
%   phi                 Euler angle - in deg.
%   psi                 Euler angle - in deg.
%   the                 Euler angle - in deg.
%  
%  OUTPUT
%   r                   ...
%
%EXAMPLE
%   r = [1 1 1]
%   r = tom_pointrotate(r,10,20,30)
%
%REFERENCES
%
%SEE ALSO
%   TOM_ROTATE3D, TOM_ROTATE2D
%
%   created by FF 08/01/03
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


phi = phi/180*pi;psi = psi/180*pi;the = the/180*pi;

matr = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0;0 0 1];

matr = matr*[1 0 0 ; 0 cos(the) -sin(the); 0 sin(the) cos(the)];

matr = matr*[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0;0 0 1];
% matr=inv(matr);
r = matr*r';
r=r';