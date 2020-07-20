#include <iostream>
#include <iomanip>
#include <fstream>
#include "FITSFile.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc != 3) {
    cerr << "Usage: " << argv[0] << " infile.fits outfile.txt\n";
    return 1;
  }

  CFITSFile file(argv[1], CFITSFile::existingro);
  CFITSImage im = file.GetImage();

  std::ofstream fout(argv[2]);
  if(! fout) {
    cerr << "Unable to create " << argv[2] << '\n';
    return 1;
  }

  const unsigned xw = im.GetXW();
  const unsigned yw = im.GetYW();

  for(int y=0; y<yw; ++y) {
    for(int x=0; x<xw; ++x) {
      fout << scientific << im.GetPixel(x, y) << ' ';
    }
    fout << '\n';
  }

  return 0;
}
