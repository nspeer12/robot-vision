#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100

int    pic[PICSIZE][PICSIZE];
double outpicx[PICSIZE][PICSIZE];
double outpicy[PICSIZE][PICSIZE];
int    edgeflag[PICSIZE][PICSIZE];
double mask[MAXMASK][MAXMASK];
double mag[PICSIZE][PICSIZE];
double finals[PICSIZE][PICSIZE];
double peaks[PICSIZE][PICSIZE];
double ival[256][256];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
int histo[PICSIZE];

int main(argc,argv)
int argc;
char **argv;
{

	int i,j,p,q,s,t,mr,centx,centy, more;
	double  maskval,sum1,sum2,sig,maxival,minival,maxval,ZEROTOL, dx, dy, slope, HI, LO, percent, cutOff, topsArea;
	FILE    *fo1, *fo2, *fo3, *fp1, *fopen();
	char *foobar;
	// Opening and closing files for reading and writing
	// file names are taken from cmd prompt

	/* input file */
	argc--; argv++;
	foobar = *argv;
	fp1=fopen(foobar,"rb");

	// convolution output image
	foobar = "output/convolution.pgm";
	fo1=fopen(foobar,"wb");
	fprintf(fo1,"P5\n257 257\n255\n");

	
	// peaks output image
	foobar = "output/peaks.pgm";
	fo2=fopen(foobar,"wb");
	fprintf(fo2,"P5\n257 257\n255\n");

	// edge detection image
	foobar = "output/edges.pgm";
	fo3=fopen(foobar,"wb"); 
	fprintf(fo3,"P5\n257 257\n255\n");

	slope = 0.0;
	sig = 1.0;
	mr = (int)(sig * 3);

        centx = (MAXMASK / 2);
        centy = (MAXMASK / 2);

        // read input picture into memory
        for (i=0;i<256;i++)
        {
        	for (j=0;j<256;j++)
                {
                    pic[i][j]  =  getc (fp1);
                }
        }

        for (p=-mr;p<=mr;p++)
        {
            for (q=-mr;q<=mr;q++)
            {
                // create x & y masks from gaussian partial derivatives
                dx = q * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig)))));
                maskx[p+centy][q+centx] = dx;

		dy = p * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig)))));
                masky[p+centy][q+centx] = dy;
            }
        }

        // compute convolution
        for (i=mr;i<=255-mr;i++)
        {
            for (j=mr;j<=255-mr;j++)
            {
                sum1 = 0;
                sum2 = 0;

                for (p=-mr;p<=mr;p++)
                {
                    for (q=-mr;q<=mr;q++)
                    {
                        sum1 += pic[i+p][j+q] * maskx[p+centy][q+centx];
                        sum2 += pic[i+p][j+q] * masky[p+centy][q+centx];
                    }
             }

             outpicx[i][j] = sum1;
             outpicy[i][j] = sum2;

            }
        }
	

        maxival = 0;
        for (i=mr;i<256-mr;i++)
        { 
            for (j=mr;j<256-mr;j++)
            {
		// take magnitude of each convolution output
                ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                        (outpicy[i][j]*outpicy[i][j])));

             if (ival[i][j] > maxival)
                maxival = ival[i][j];

             }
        }

	// output convolution image
        for (i=0;i<256;i++)
        {
            for (j=0;j<256;j++)
            {
                mag[i][j] = ival[i][j];
                ival[i][j] = (ival[i][j]/maxival) * 255;
                fprintf(fo1,"%c",(char)((int)(ival[i][j])));
            }
            
            fprintf(fo1,"\n");
          }

	// ??? need explaination
	// compute peaks
        for (i=mr;i<256-mr;i++)
        {
            for (j=mr;j<256-mr;j++)
            {
                if(outpicx[i][j] == 0.0)
                {
                    outpicx[i][j] == 0.00001;
                }

            slope = (double)(outpicy[i][j])/(double)(outpicx[i][j]);

            if((slope <= 0.4142) && (slope > -0.4142))
            {
                if((ival[i][j] > ival[i][j-1]) && (ival[i][j] > ival[i][j+1]))
                {
                    peaks[i][j] = 255;
                }
            }
            else if((slope <= 2.4142) && (slope > 0.4142))
            {
                if((ival[i][j] > ival[i-1][j-1]) && (ival[i][j] > ival[i+1][j+1]))
                {
                    peaks[i][j] = 255;

                }
            }
            else if((slope <= -0.4142) && (slope > -2.4142))
            {
                if((ival[i][j] > ival[i+1][j-1]) && (ival[i][j] > ival[i-1][j+1]))
                {
                    peaks[i][j] = 255;

                }
            }
            else
            {
                if((ival[i][j] > ival[i-1][j]) && (ival[i][j] > ival[i+1][j]))
                {
                    peaks[i][j] = 255;
                }
            }

           }
        }

	// find maxival for peaks image
        maxival = 0;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {

              if (peaks[i][j] > maxival)
                maxival = peaks[i][j];
           }
        }
	
	// output peaks image
        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
                peaks[i][j] = (peaks[i][j]/maxival) * 255;
                fprintf(fo2,"%c",(char)((int)(peaks[i][j])));
            }
            fprintf(fo2,"\n");
          }

        /* double thresholding */
        percent = 0.05;
        topsArea = 0;

        // histogram
        for (i=mr;i<256-mr;i++)
        {
            for (j=mr;j<256-mr;j++)
            {
                histo[(int)(ival[i][j])]++;
            }
        }

        cutOff = percent * 255 * 255;
	
	// find threshold
        for (i = 255; i > 0; i--)
        {
            topsArea += histo[i];
            if(topsArea > cutOff)
                break;
        }

        HI = i;
        LO = 0.35 * HI;
	
	// output thresholds
        printf("HI = %lf, LO = %lf \n",HI, LO);

        more = 0;

	// find edge values above HI
        for (i=0;i<256;i++)
        {
            for (j=0;j<256;j++)
            {
                if(peaks[i][j] > 254)
                {
                    if(ival[i][j] > HI)
                    {
                        peaks[i][j] = 0;
                        finals[i][j] = 255;
                    }
                    else if(ival[i][j] < LO)
                    {
                        peaks[i][j] = 0;
                        finals[i][j] = 0;
                    }
                }
            }
         }

	// ??? find mid values
        more = 1;
        while(more == 1)
        {
            more = 0;
            for (i=0;i<256;i++)
            {
                for (j=0;j<256;j++)
                {
                    if(peaks[i][j] > 254)
                    {
                        for (p=-mr;p<=mr;p++)
                        {
                            for (q=-mr;q<=mr;q++)
                            {
                                if(finals[i+p][j+q] > 254)
                                {
                                    peaks[i][j] = 0;
                                    finals[i][j] = 255;
                                    more = 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        // find maxival for normalization
        maxival = 0;
        for (i=mr;i<256-mr;i++)
        {
            for (j=mr;j<256-mr;j++)
            {

                if (finals[i][j] > maxival)
                    maxival = finals[i][j];

           }
        }
        
        // normalize pixels and write image
        for (i=0;i<256;i++)
        {
            for (j=0;j<256;j++)
            {
                finals[i][j] = (finals[i][j]/maxival) * 255;
                fprintf(fo3,"%c",(char)((int)(finals[i][j])));
            }
            
            fprintf(fo3,"\n");
         }
}
