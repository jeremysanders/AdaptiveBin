#ifndef FITSFILE_NEW_H
#define FITSFILE_NEW_H

// Main FITS File routines

//      FITS Image manipulation library
//      Copyright (C) 2000 Jeremy Sanders
//      Contact: jss@ast.cam.ac.uk
//               Institute of Astronomy, Madingley Road,
//               Cambridge, CB3 0HA, UK.

//      See the file COPYING for full licence details.

//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.

//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.

//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

#include "FITSGeneral.h"
#include "FITSImage.h"
#include "FITSPosn.h"

class CFITSFile
{
public:                  // public data types
  enum COpenMode {existingro, existing, create};
    // modes: existing file, read-only
    //        existing file, read-write
    //        create new file
  enum CDataType {tint, tfloat, tstring};

public:                  // public methods
  CFITSFile();
  CFITSFile(const char *filename, COpenMode openMode);
  ~CFITSFile();

  void OpenFile(const char *fileName, COpenMode openMode);
  void CloseFile();
  void CheckStatus(const char *whereMessage);

  void WriteKey(const char *key, CDataType datatype, const void *data,
		const char *comment=0);
  void UpdateKey(const char *key, CDataType datatype, const void *data,
		 const char *comment=0);
    // can set comment=0
  int ReadKey(const char *key, CDataType datatype, void *data,
	      const void *defdata=0);
    // returns 0 if key not found (and defdata is provided)

  void ReadImage();
    // (re)reads the image from the file into m_image
  void ReadImageInclNull(CFloatType nullval = CNullValue);
    // does above, but sets NANs to nullval
  void WriteImage();
    // updates the image in the file
  void WriteImageInclNull(CFloatType nullval = CNullValue);
    // write NULLs as NULLs

  void WriteHistory(const char *hist);
    // append hist as line in history
  void UpdateChecksum();
    // update fits file checksum

  void SetImage(const CFITSImage &copy);
  CFITSImage & GetImage();
  void SetPosn(const CFITSPosn &copy);
  CFITSPosn & GetPosn();

public:
  static char m_comment_char;
  
private: // private methods
  int GetFITSDataType(CDataType d);
    // convert CDataType to fitsio int type

private: // private data
  fitsfile *m_file;

  char m_fileName[256];
  int m_status;
  int m_previousWritten;

  int m_FITSMode;

  CFITSImage m_image;
  CFITSPosn m_posnImage;
};

#endif
