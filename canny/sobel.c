/* Sobel.c */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>

int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];

/* convolution masks */
int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
double ival[256][256],maxival;

int main(argc, argv)
int argc;
char **argv;
{
    int i,j,p,q,mr,sum1,sum2;
    double threshold;        // lo threshold
    double hiThreshold = 100;

    FILE *fo1, *fo2, *fp1, *fopen();
    char *foobar;
       
    int magnitude = 0;

    /* input file */
    argc--; argv++;
    foobar = *argv;
    fp1=fopen(foobar,"rb");

    /* output file */
	argc--; argv++;
	foobar = *argv;
	fo1=fopen(foobar,"wb");
    
    int rows = 256;
    int cols = 256;
    
    /* write PGM header */
    fprintf(fo1, "P5\n%d %d\n255\n", rows, cols);
        
    argc--; argv++;
	foobar = *argv;

	/* threshold value */
    threshold = atof(foobar);

    // scan in (optional) Hi Threshold
    argc--; argv++;
    foobar = *argv;
    
       
    if (foobar == NULL)
    {
        // set hiThreshold arbitrarily high if no input specified
        hiThreshold = INT_MAX;
    }
    else
    {
         hiThreshold = atof(foobar);
    }
 
    /* scan in image */
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            pic[i][j]  =  getc (fp1);

            /* legacy octal eight 1s bitmask for reading chars */
            pic[i][j]  &= 0377;

            magnitude += pic[i][j];
        }
    }
    
    /* convolution */
    
    mr = 1; //  mask radius

    for (i=mr;i<256-mr;i++)
    {
        for (j=mr;j<256-mr;j++)
        {
            sum1 = 0;
            sum2 = 0;
            for (p=-mr;p<=mr;p++)
            { 
                for (q=-mr;q<=mr;q++)
                {
                    sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
                    sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
                }
            }
            
            outpicx[i][j] = sum1;
            outpicy[i][j] = sum2;
        }
    }
    
    int diff;
    maxival = 0;
    for (i=mr;i<256-mr;i++)
    {
        for (j=mr;j<256-mr;j++)
        {
            // magnitude vector
            diff = sqrt((double)((outpicx[i][j]*outpicx[i][j]) + (outpicy[i][j]*outpicy[i][j])));
            
            // is the magnitude of the jump big enough?            
            ival[i][j] = diff > threshold && diff < hiThreshold ? 1 : 0;
    
            if (ival[i][j] > maxival)
                maxival = ival[i][j];

        }
    }
    
    // normalize pixels and output image
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
             ival[i][j] = (ival[i][j] / maxival) * 255;            
             fprintf(fo1,"%c",(char)((int)(ival[i][j])));
             
        }
    }
    
    printf("Magnitude: %d\n", magnitude);
    printf("Hi Threshold: %.2f\n", hiThreshold);
    printf("Lo Threshold: %.2f\n", threshold);
}
