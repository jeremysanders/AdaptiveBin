#include <iostream>
#include <string>
#include <cdf.h>
#include "FITSFile.h"

int main(int argc, char *argv[])
{
  if(argc != 3) {
    cerr << "Invalid number of parameters\n"
      "Usage: FITS2CDF in.fits out.cdf\n";
    exit(1);
  }

  CFITSFile file(argv[1], CFITSFile::existingro);
  CFITSImage im = file.GetImage();

  CDFid id;
  CDFstatus status;

  long dim[2];
  dim[0] = im.GetYW();
  dim[1] = im.GetXW();

  status = CDFcreate(argv[2], 2, dim, NETWORK_ENCODING,
		     ROW_MAJOR, &id);

  if(status != CDF_OK) {
    cerr << "Could not create CDF file\n";
    exit(1);
  }

  char name[128];
  strcpy(name, "data");

  long dim_vary[] = {VARY, VARY};
  status = CDFvarCreate(id, name, CDF_REAL8, 1, VARY, dim_vary,
			im.GetImageBuffer());


  status = CDFclose(id);
  if(status != CDF_OK) {
    cerr << "Could not close CDF file\n";
    exit(1);
  }

}

