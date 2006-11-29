#include <volume_io.h>

#define KERN_UNDEF    0
#define KERN_GAUSSIAN 1
#define KERN_RECT     2


float normal_dist(float c, float fwhm, float mu, float x)
{
  float sigma,t1,t2,t3,f;
  
  sigma = fwhm/2.35482;
  
  if (sigma==0) {
    if (x==mu)
      f = c / (sqrt(2*PI) * sigma);      /* !!!!!!!!!!! */
    else
      f = 0;
  }
  else {
    t1 = c / (sqrt(2*PI) * sigma);
    t2 = (x-mu)*(x-mu)/(2*sigma*sigma);
    t3 = exp(-t2);
    f = t1*t3;
  }
  
  return(f);
}


float rect_dist(float c, float fwhm, float mu, float x)
{
  
  float t;
  
  t = x-mu;
  
  if ( t >= -fwhm/2  && t < fwhm/2 ) {
    return(  (float) (c/fwhm) );
  }
  else
    return ( (float) 0.0 );
  
}


void make_kernel(float *kern, float vsize, float fwhm, int size, int type)
{

  int kindex,k;
  float c,r,max;

  (void)memset(kern,(int)0,(size_t)((2*size+1)*sizeof(float)));
  
  for ( k = -size/2; k<size/2; ++k) {

    kindex = ((k + size) % size)*2 +1; 


    switch (type) {
    case KERN_GAUSSIAN: kern[kindex] = normal_dist(1.0*vsize,fwhm,0.0,(float)(vsize*k)); break;
    case KERN_RECT:     kern[kindex] = rect_dist(1.0*vsize,fwhm,0.0,(float)(vsize*k)); break;
    default: 
      {
        (void) fprintf (stderr,"Illegal kernel type = %d\n",type);
        (void) fprintf (stderr,"Impossible error in make_kernel(), line %d of %s\n",
                        __LINE__,__FILE__);
        k = size/2;
      }
    }
  }


}

char *prog_name;

int main(int argc, char *argv[])
{
  float
    sum, sum2,
    vsize,
    fwhm,
    complex_kernel[512];
  int
    i, nvoxels;



  prog_name = argv[0];

  /* Check arguments */
  if (argc != 4) {
    (void) fprintf(stderr, "Usage: %s fwhm(mm) vsize(mm) nvoxels\n",
                   argv[0]);
    exit(EXIT_FAILURE);
  }
  
  fwhm   = (float)atof(argv[1]);
  vsize  = (float)atof(argv[2]);
  nvoxels= atoi(argv[3]);

  print ("for fwhm = %f, vsize = %f, nvoxels = %d\n", fwhm, vsize, nvoxels);

  make_kernel(complex_kernel,(float)vsize,fwhm,nvoxels,KERN_GAUSSIAN);

  sum2, sum = 0.0;
  for(i=0; i<nvoxels; i++) {
    print ("%f\n",complex_kernel[i*2 + 1]);
    sum += complex_kernel[i*2 + 1];
    sum2 += complex_kernel[i*2 + 1] * complex_kernel[i*2 + 1];
  }

  print ("sum = %f, sum^2 = %f\n",sum,sum2);

  exit(EXIT_SUCCESS);
}
