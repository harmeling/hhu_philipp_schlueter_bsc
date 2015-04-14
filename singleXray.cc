/*
 Copyright (C) 2015 Philipp Schlueter

 This program is free software; you can redistribute it and/or modify 
 it under the terms of the GNU General Public License as published by 
 the Free Software Foundation; either version 3 of the License, or 
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but 
 WITHOUT ANY WARRANTY; without even the implied warranty of 
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <octave/oct.h>
#include <math.h>

DEFUN_DLD (singleXray, args, , "\
 -- Function File: XT = Xray(IMAGESTACK, THETA, PHI, MJD)\n\
\n\
     Calculates the 3D-Xray transform of the 3D matrix IMAGESTACK at\n\
     angles given in THETA and PHI. In MJD you can declare the date\n\
     the images in IMAGESTACK were made. THETA and PHI are expected\n\
     to be in units of radians. \n")
{
  int nargin = args.length ();

  if (nargin != 4)
    print_usage ();
  else
    {

    NDArray A = args(0).array_value ();
    NDArray THETA = args(1).array_value();
    NDArray PHI = args(2).array_value();
    NDArray MJD = args(3).array_value();
      
    dim_vector dA = A.dims();
    dim_vector dTHETA = THETA.dims();
    dim_vector dPHI = PHI.dims();
    float zextent = MJD(dA(2)-1)-MJD(0)+1;
    int b = ceil(sqrt((dA(0)*dA(0)+dA(1)*dA(1) + zextent * zextent ) ) / 2 + 1);
    float theta = THETA.elem(0);
    float phi = PHI.elem(0);

    Matrix singleout(2*b+1,2*b+1);
    Matrix normalizer(2*b+1,2*b+1);
    for(unsigned int i = 0; i<2*b+1; i++)
      for(unsigned int j = 0; j<2*b+1; j++){
        singleout(i,j)=NAN;
        normalizer(i,j)=NAN;
      }
      
    float xc = (dA(0)-1.0)/2.0;
    float yc = (dA(1)-1.0)/2.0;
    float zc = (MJD(dA(2)-1)+MJD(0))/2.0;
    float xoffset,yoffset,xfrac,yfrac;
    int xk,yk;
     
    for(unsigned int z = 0; z<dA(2); z++)
      for(unsigned int x = 0; x<dA(0); x++)
        for(unsigned int y = 0; y<dA(1); y++){
          //Calulate Offsets
          xoffset = ((x-xc)*cos(theta)+(MJD(z)-zc)*sin(theta)+b+1);
          yoffset = (-sin(phi)*sin(theta)*(x-xc)+cos(phi)*(y-yc)+sin(phi)*cos(theta)*(MJD(z)-zc)+b+1);
          xk = floor(xoffset);
          yk = floor(yoffset);
          xfrac = xoffset-xk;
          yfrac = yoffset-yk;
			
			    //Distribute pixelvalue to four nearest grit cells
          if(singleout(xk,yk) != singleout(xk,yk) ) {
            singleout(xk,yk) = 0;
            normalizer(xk,yk) = 0;
          } 
          singleout(xk,yk) += A(x,y,z)*(1-xfrac)*(1-yfrac);
          normalizer(xk,yk) += (1-xfrac)*(1-yfrac);
		
          if(singleout(xk,yk+1) != singleout(xk,yk+1)) {
            singleout(xk,yk+1) = 0;
            normalizer(xk,yk+1) = 0;
          }
          singleout(xk,yk+1) += A(x,y,z) * (1-xfrac) * yfrac;
          normalizer(xk,yk+1) += (1-xfrac) * yfrac;
			
          if(singleout(xk+1,yk) != singleout(xk+1,yk) ) {
            singleout(xk+1,yk) = 0;
            normalizer(xk+1,yk) = 0;
          } 
          singleout(xk+1,yk) += A(x,y,z)*xfrac*(1-yfrac);
          normalizer(xk+1,yk) += xfrac*(1-yfrac);
			
          if(singleout(xk+1,yk+1) != singleout(xk+1,yk+1) ) {
          singleout(xk+1,yk+1) = 0;
	    normalizer(xk+1,yk+1) = 0;
	  } 
	  singleout(xk+1,yk+1) += A(x,y,z)*xfrac*yfrac;
	  normalizer(xk+1,yk+1) += xfrac*yfrac;	
	 }
      
	for(unsigned int i = 0; i<2*b+1; i++)
	  for(unsigned int j = 0; j<2*b+1; j++){
	    singleout(i,j)/=normalizer(i,j);
	  }
	
    if (! error_state)
      return octave_value (singleout);
    }

  return octave_value_list ();
}
