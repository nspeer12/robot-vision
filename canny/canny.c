#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#include <stdlib.h>

#define  PICSIZE 256
#define  MAXMASK 100

         int    pic[PICSIZE][PICSIZE];
         double outpic1[PICSIZE][PICSIZE];
         double outpic2[PICSIZE][PICSIZE];
         int    edgeflag[PICSIZE][PICSIZE];
         double mask[MAXMASK][MAXMASK];
         double maskx[MAXMASK][MAXMASK];
	 double masky[MAXMASK][MAXMASK];
	 double conv[PICSIZE][PICSIZE];

void main(argc,argv)
int argc;
char **argv;
{
	/* TODO:
	 x marrh.c uses second derivatives of gauss, use first derivatives instead
	 x include partial derivative in respect to y
	 * keep the convolution code, but double up on convolution for x and y
	 * delete everything below the convolution code
	 * go to sobel code, bring in the sqrt of squares code and the scaling of the output and maxival
	 */
	double e = 2.71828;

	int     i,j,p,q,s,t,mr,centx,centy;
        double  maskval,sum,sig,maxival,minival,maxval,ZEROTOL;
        FILE    *fo1, *fo2,*fp1, *fopen();
        char    *foobar;

        argc--; argv++;
        foobar = *argv;
        printf("%s\n", foobar);
	fp1=fopen(foobar,"rb");

        argc--; argv++;
        foobar = *argv;
        printf("%s\n", foobar);
	fo1=fopen(foobar,"wb");

        argc--; argv++;
        foobar = *argv;
        sig = atof(foobar);

	int rows, cols = 256;
	/* write PGM header */
    	fprintf(fo1, "P5\n256 256\n255\n", rows, cols);


        mr = (int)(sig * 3);
        centx = (MAXMASK / 2);
        centy = (MAXMASK / 2);

        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fp1);
                  //printf("%c\n", pic[i][j]);
		}
        }
	
	// probably need an x mask and y mask, the convolve each
	
        for (p=-mr;p<=mr;p++)
        {  for (q=-mr;q<=mr;q++)
           {
	      /* change this equation */
              //maskval = ((2-(((p*p)+(q*q))/(sig*sig))) * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
              
	      // p = y, q = x
	      maskval = (q) * exp(((-1*q*q)+(p*p))/(2*sig*sig));
	      (maskx[p+centy][q+centx]) = maskval;
           }
        }
	
       	for (p=-mr;p<=mr;p++)
        {  for (q=-mr;q<=mr;q++)
           {
	      /* change this equation */
              //maskval = ((2-(((p*p)+(q*q))/(sig*sig))) * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
              
	      // p = y, q = x
	      
	      maskval += (p) * exp(((-1*q*q)+(p*p))/(2*sig*sig));
	      (masky[p+centy][q+centx]) = maskval;
           }
        }

	/* compute convolutions */
        for (i=mr;i<=255-mr;i++)
        { for (j=mr;j<=255-mr;j++)
          {
             sum = 0;
             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {
                   sum += pic[i+p][j+q] * maskx[p+centy][q+centx];
		   sum += pic[i+p][j+q] * masky[p+centy][q+centx];

                }
             }
             outpic1[i][j] = sum;
             conv[i][j] = sum;
          }
	}

	// lets output the image and see what happens
	int diff;
    	maxival = 0;
    	for (i=mr;i<256-mr;i++)
   	{
		for (j=mr;j<256-mr;j++)
		{
	    		// magnitude vector
	    		//diff = sqrt((double)((outpicx[i][j]*outpicx[i][j]) + (outpicy[i][j]*outpicy[i][j])));
	    
	    		// is the magnitude of the jump big enough?            
	    		//ival[i][j] = diff > threshold && diff < hiThreshold ? 1 : 0;
    
	    		if (outpic1[i][j] > maxival)
				maxival = outpic1[i][j];

		}
    	}
 

    // normalize pixels and output image
    for (i=0;i<256;i++)
    {
	for (j=0;j<256;j++)
	{
	     pic[i][j] = (pic[i][j] / maxival) * 255;            
	     fprintf(fo1,"%c",(char)((int)(pic[i][j])));
	     
	}
    }

}

