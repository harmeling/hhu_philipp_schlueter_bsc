## Copyright (C) 2015 Philipp Schlueter
##
## This program is free software; you can redistribute it and/or modify 
## it under the terms of the GNU General Public License as published by 
## the Free Software Foundation; either version 3 of the License, or 
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but 
## WITHOUT ANY WARRANTY; without even the implied warranty of 
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

## -- Function File: [XT,XP] = Xray(IMAGESTACK, THETA, PHI, MJD)
## -- Function File: [XT,XP] = Xray(IMAGESTACK, THETA, PHI)
## -- Function File: [XT,XP] = Xray(IMAGESTACK, THETA)
## -- Function File: [XT,XP] = Xray(IMAGESTACK)
##
##     Calculates the 3D-Xray transform of the 3D matrix IMAGESTACK at
##     angles given in THETA and PHI.  To each element of THETA and PHI
##     corresponds a matrix in XT.  In MJD you can declare the date the
##     images in IMAGESTACK were made.  The variable XP represents the
##     x-axis and y-axis of the rotated coordinate.  If THETA is not
##     defined, then -90:90 is assumed.  If PHI is not defined, then
##     -90:90 is assumed.  If MJD is not defined, then it is assumed,
##     that all images are separated by one.


function [XT,xp] = Xray(Imagestack,theta,phi, mjd)

  ## Input checking
  if (nargin == 0 || nargin > 4)
    print_usage ();
  elseif (nargin == 3)
    mjd = NaN;
  elseif (nargin == 2)
    mjd = NaN;
    phi = -90:90;
  elseif (nargin == 1)
    mjd = NaN;
    phi = -90:90;
    theta = -90:90;
  endif
  
  if (!ismatrix(Imagestack) || ndims(Imagestack) != 3)
    error('First input must be a 3d "Matrix"');
  endif
	
  [xlength,ylength,zlength]=size(Imagestack);
  if (isnan(mjd))
    mjd = 1:zlength;
  end

  if not(length(mjd)==zlength)
    error('Size of mjd vector does not equal number of images.');
  end

  ## Output allocation
  b = ceil(sqrt(sum(xlength^2+ylength^2+(max(mjd)-min(mjd)+1)^2))/2+1);
  xp = [-b:b]';

  RT = zeros(2*b+1,2*b+1,length(theta),length(phi));

  th = theta'*pi/180;
  ph = phi'*pi/180;
  pkg load ndpar;

  for p = 1: length(phi)
    XT(:,:,:,p) = ndpar_arrayfun(nproc,@singleXray,Imagestack,th,ph(p),mjd,"IdxDimensions",[0,1,0,0],"CatDimensions", [3]);
  endfor
endfunction


