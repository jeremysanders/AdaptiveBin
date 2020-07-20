#include <iostream>
#include <cstdlib>
#include <vector>

#include <FITSFile.h>
#include <FITSImage.h>

using namespace std;

unsigned get_no_bins(const CFITSImage& im)
{
  int maxbin = -1;

  for(int y=0; y<im.GetYW(); ++y)
    for(int x=0; x<im.GetXW(); ++x)
      {
	const int b = int( im.GetPixel(x, y) );
	if( b > maxbin ) maxbin = b;
      }

  return unsigned(maxbin+1);
}

void countbins(const CFITSImage& im, vector<unsigned>* counts)
{
  for(int y=0; y<im.GetYW(); ++y)
    for(int x=0; x<im.GetXW(); ++x)
      {
	const int b = int( im.GetPixel(x, y) );

	if( b >= 0 )
	    (*counts)[b] ++;
      }

}

void paintoutput(const CFITSImage& binmap, const vector<unsigned>& counts,
		 CFITSImage* outimage)
{
  const double nan = 0./0.;
  for(int y=0; y<binmap.GetYW(); ++y)
    for(int x=0; x<binmap.GetXW(); ++x)
      {
	const int bin = int( binmap.GetPixel(x, y) );
	if( bin >= 0 )
	    outimage->SetPixel(x, y, counts[bin]);
	else
	    outimage->SetPixel(x, y, nan);
      }
}

int main(int argc, char *argv[])
{
  if( argc != 3 )
    {
      cerr << "Usage: " << argv[0] << " binmap.fits countmap.fits\n";
      exit(1);
    }

  CFITSFile infile(argv[1], CFITSFile::existingro);
  const CFITSImage inimage = infile.GetImage();

  const unsigned nobins = get_no_bins(inimage);
  vector<unsigned> bincounts(nobins, 0);
  countbins(inimage, &bincounts);

  CFITSImage outimage(inimage.GetXW(), inimage.GetYW());

  paintoutput(inimage, bincounts, &outimage);

  CFITSFile outfile(argv[2], CFITSFile::create);
  outfile.SetImage(outimage);
  outfile.WriteImage();
}
