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

function res = findStar(Imagestack,MJD = -99999)
  [~,~,zlength]=size(Imagestack);
  if (MJD ==  -99999)
	  MJD = 1:zlength;
  endif

  if not(length(MJD)==zlength)
    error('Size of MJD vector does not equal number of images.');
  endif

MJD = MJD-mean(MJD);

  pkg load image;
  pkg load statistics;
  phi = linspace(-89.99,89.99,11);
  theta = linspace(-89.99,89.99,11);
  xOffset = 0;
  yOffset = 0;
  xLowerLimit = 1;
  yLowerLimit = 1;

  minmax = @(A) [min(A(:)) max(A(:))];

  for i = 1:15
    xOffset = xOffset+xLowerLimit-1;
    yOffset = yOffset+yLowerLimit-1;
    [xlength,ylength,~]=size(Imagestack);
    xc = (xlength+1)/2;
    yc = (ylength+1)/2;
    zc = (min(MJD)+max(MJD))/2;
    zlength = max(MJD)-min(MJD)+1;

    RT = Xray(Imagestack,theta,phi,MJD);
    [RTxlength,~,~,~] = size(RT);
    b = (RTxlength-1)/2;

    RTc = convn(convn(RT,[0.0001,1,0.0001],'same'),[0.0001;1;0.0001],'same');
    RTc = convn(convn(RTc,normpdf(1:9,5,3)','same'),normpdf(1:9,5,3),'same');
    [m,im]=max(RTc(:));
    [i1,i2,i3,i4]=ind2sub(size(RTc),im);

    th = [theta(max(i3-2,1)) theta(min(i3+2,length(theta)))]*pi/180;
    ph = [phi(max(i4-2,1)) phi(min(i4+2,length(phi)))]*pi/180;
    vx = minmax(-tan(th));
    vy = minmax(-tan(ph)'*(1./cos(th)));
    x0 = minmax((i1-(min(MJD)-zc)*sin(th)-b-1)'*(1./cos(th))+xc-1);
    [X,Y,Z] = meshgrid(1:2,1:2,1:2);
    y0 = minmax((i2+sin(ph(X(:))).*sin(th(Y(:))).*(x0(Z(:))-xc)-b-1-sin(ph(X(:))).*cos(th(Y(:)))*(min(MJD)-zc))./cos(ph(X(:)))+yc-1);

    xLowerLimit = max(1,min(floor([x0(1),x0(1)+vx*zlength]-10)));
    xUpperLimit = min(xlength,max(ceil([x0(2),x0(2)+vx*zlength]+10)));
    yLowerLimit = max(1,min(floor([y0(1),y0(1)+vy*zlength]-10)));
    yUpperLimit = min(ylength,max(ceil([y0(2),y0(2)+vy*zlength]+10)));

    Imagestack = Imagestack(xLowerLimit:xUpperLimit,yLowerLimit:yUpperLimit,:);
    theta = linspace(theta(max(i3-2,1)),theta(min(i3+2,length(theta))),7);
    phi = linspace(phi(max(i4-2,1)),phi(min(i4+2,length(phi))),7);
  endfor
 
  M = nanmax((RT(:,:,i3,i4)-RTc(:,:,i3,i4))(:));
  M2 = nanmax((-(RT(:,:,i3,i4)-RTc(:,:,i3,i4))(:)));
  res = [ sum(x0)/2+xOffset std(x0);
          sum(y0)/2+yOffset std(y0);
          sum(vx)/2 std(vx);
          sum(vy)/2 std(vy);
          M M2];
endfunction
