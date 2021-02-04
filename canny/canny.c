#include <stdio.h>
#include <math.h>
#define    PICSIZE 256
#define    MAXMASK 100

int     pic[PICSIZE][PICSIZE];
double outpicx[PICSIZE][PICSIZE];
double outpicy[PICSIZE][PICSIZE];
double mask[MAXMASK][MAXMASK];
double magnitude[PICSIZE][PICSIZE];
double final[PICSIZE][PICSIZE];
double peaks[PICSIZE][PICSIZE];
double ival[256][256];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
int histogram[PICSIZE];

int main(argc,argv)
int argc;
char **argv;
{
    int i, j, p, q, s, t, mr, centx, centy, more;
    double  maskval, sum1, sum2, sig, maxival, minival, maxval, ZEROTOL, dx, dy, slope, HI, LO, percent, cutOff, hiArea;
    FILE *fo1, *fo2, *fo3, *fp1, *fopen();
    char *foobar;

    // input image
    argc--; argv++;
    foobar = *argv;
    if (foobar == NULL)
    {
        printf("Error: specify input image\n");
        return 0;
    }

    fp1=fopen(foobar,"rb");
    foobar = NULL;

    // convolution output image
    foobar = "output/convolution.pgm";
    fo1=fopen(foobar,"wb");
    fprintf(fo1,"P5\n256 256\n255\n");


    // peaks output image
    foobar = "output/peaks.pgm";
    fo2=fopen(foobar,"wb");
    fprintf(fo2,"P5\n256 256\n255\n");


    // edge detection image
    foobar = "output/edges.pgm";
    fo3=fopen(foobar,"wb");
    fprintf(fo3,"P5\n256 256\n255\n");

    slope = 0.0;
    sig = 1.0;
    mr = (int)(sig * 3);

    centx = (MAXMASK / 2);
    centy = (MAXMASK / 2);

    // get rid of input file header
    char blank[20];
    for (i=0;i<4;i++)
    {
        fscanf(fp1,"%s", blank);
    }

    // read input picture into memory
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {

             pic[i][j] = getc (fp1);
        }
    }

    // create x & y masks from gaussian partial derivatives
    for (p=-mr;p<=mr;p++)
    {
        for (q=-mr;q<=mr;q++)
        {
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
            magnitude[i][j] = ival[i][j];
            ival[i][j] = (ival[i][j]/maxival) * 255;
            fprintf(fo1,"%c",(char)((int)(ival[i][j])));
        }
    }


    // compute peaks
    for (i=mr;i<256-mr;i++)
    {
        for (j=mr;j<256-mr;j++)
        {
            if(outpicx[i][j] == 0.0)
            {
                outpicx[i][j] = 0.00001;
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

    // normalize peaks
    maxival = 0;
    for (i=mr;i<256-mr;i++)
    {
        for (j=mr;j<256-mr;j++)
        {

            if (peaks[i][j] > maxival)
                maxival = peaks[i][j];
        }
    }

    // output peaks image
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            peaks[i][j] = (peaks[i][j]/maxival) * 255;
            fprintf(fo2,"%c",(char)((int)(peaks[i][j])));
        }
    }

    /* find HI threshold */

    // hi percentage
    percent = 0.05;

    // create histogram
    for (i=mr;i<256-mr;i++)
    {
        for (j=mr;j<256-mr;j++)
        {
            histogram[(int)(ival[i][j])]++;
        }
    }

    hiArea = 0;

    cutOff = percent * 255 * 255;

    // find threshold
    for (i = 255; i > 0; i--)
    {
        hiArea += histogram[i];

        if (hiArea > cutOff)
            break;
    }

    // take percentage of threshold
    HI = i;
    LO = 0.35 * HI;

    // take edge values above HI and discard below LO
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            if(peaks[i][j] >= 255)
            {
                if(ival[i][j] > HI)
                {
                    peaks[i][j] = 0;
                    final[i][j] = 255;
                }
                else if(ival[i][j] < LO)
                {
                    peaks[i][j] = 0;
                    final[i][j] = 0;
                }
            }
        }
    }


    // find middle edges
    more = 1;
    while (more == 1)
    {
        more = 0;
        for (i=0;i<256;i++)
        {
            for (j=0;j<256;j++)
            {
                if(peaks[i][j] >= 255)
                {
                    for (p=-mr;p<=mr;p++)
                    {
                        for (q=-mr;q<=mr;q++)
                        {
                            if(final[i+p][j+q] >= 255)
                            {
                                peaks[i][j] = 0;
                                final[i][j] = 255;
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
            if (final[i][j] > maxival)
                maxival = final[i][j];
        }
  }

  // normalize pixels and write image
  for (i=0;i<256;i++)
  {
        for (j=0;j<256;j++)
        {
             final[i][j] = (final[i][j]/maxival) * 255;
             fprintf(fo3,"%c", (char)((final[i][j])));
        }
    }
}
