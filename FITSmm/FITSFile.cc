// FITS File handling

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

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "FITSFile.h"

char CFITSFile::m_comment_char = ':';

void cf_comment()
{
  printf("%c ", CFITSFile::m_comment_char);
}

CFITSFile::CFITSFile()
{
  m_file = 0;
  m_fileName[0] = 0;
  m_status = 0;
  m_previousWritten = 0;
}

CFITSFile::CFITSFile(const char *filename, COpenMode openMode)
{
  m_file = 0;
  m_fileName[0] = 0;
  m_status = 0;
  m_previousWritten = 0;

  OpenFile(filename, openMode);
}

CFITSFile::~CFITSFile()
{
  CloseFile();
}

void CFITSFile::OpenFile(const char *fileName, COpenMode openMode)
{
  if( m_file != 0 ) {
    fprintf(stderr, "*  CFITSFile::OpenFile failed: CFITSFile already open\n");
    exit(-1);
  }

  m_file = 0;
  m_status = 0;
  strcpy(m_fileName, fileName);

  switch(openMode) {

  case existingro:
    cf_comment();
    printf("Opening %s (RO)\n", m_fileName);
    m_FITSMode = READONLY;
    m_previousWritten = 1;
    fits_open_file(&m_file, m_fileName, m_FITSMode, &m_status);
    CheckStatus("Opening file (RO)");
    ReadImage();
    break;

  case existing:
    cf_comment();
    printf("Opening %s (RW)\n", m_fileName);
    m_FITSMode = READWRITE;
    m_previousWritten=  1;
    fits_open_file(&m_file, m_fileName, m_FITSMode, &m_status);
    CheckStatus("Opening file (RW)");
    ReadImage();
    break;
  case create:
    cf_comment();
    printf("Creating %s\n", m_fileName);
    m_FITSMode = READWRITE;
    m_previousWritten = 0;
    unlink(m_fileName);  // delete old file!!
    fits_create_file(&m_file, m_fileName, &m_status);
    CheckStatus("Creating file");
    break;
  }

}

void CFITSFile::CloseFile()
{
  if( m_file != 0 ) {
    cf_comment();
    printf("Closing %s\n", m_fileName);

    fits_close_file(m_file, &m_status);
    m_file = 0;
    m_fileName[0] = 0;
    m_status = 0;
  }
}

void CFITSFile::CheckStatus(const char *whereMessage)
{
  if( m_status != 0 ) {
    fprintf(stderr, "*  CFITSFile::CheckStatus failed in %s\n",
	    whereMessage);
    
    char buffer[256];
    fits_get_errstatus(m_status, buffer);
    fprintf(stderr, "*   FITS error: %s, %i\n", buffer, m_status);
    fprintf(stderr, "*   File name: %s\n", m_fileName);
    exit(-1);
  }
}

int CFITSFile::GetFITSDataType(CDataType d)
{
  switch(d) {
  case tint:
    return TINT;
  case tfloat:
    return CFITSFloatType;
  case tstring:
    return TSTRING;
  }

  return TINT;
}

void CFITSFile::WriteKey(const char *key, CDataType datatype,
			 const void *data, const char *comment)
{
  fits_write_key(m_file, GetFITSDataType(datatype), (char*)key,
		 (void*)data, (char*)comment, &m_status);

  CheckStatus("CFITSFile::WriteKey");
}

void CFITSFile::UpdateKey(const char *key, CDataType datatype,
			  const void *data, const char *comment)
{
  fits_update_key(m_file, GetFITSDataType(datatype), (char*)key,
		  (void*)data, (char*)comment, &m_status);

  CheckStatus("CFITSFile::UpdateKey");
}

int CFITSFile::ReadKey(const char *key, CDataType datatype, void
		       *data, const void *defdata)
{
//    int stat;
//    char card[128];
//    fits_read_record(m_file, 0, card, &stat);

  int retval = fits_read_key(m_file, GetFITSDataType(datatype),
		    (char*)key, data, 0, &m_status);

  //printf("%i\n", retval);

  if(retval == 0)
    return 1;
  else {
    // return default value
    if(defdata != 0) {
      switch(datatype) {
      case tint:
	*((int*)data)=*((int*)defdata);
	break;
      case tfloat:
	*((CFloatType*)data)=*((CFloatType*)defdata);
	break;
      case tstring:
	strcpy( (char*)data, (const char *)defdata );
	break;
      }
      
      fits_clear_errmsg();
      m_status = 0;
      return 0;
    }

    char buffer[128];
    sprintf(buffer, "CFITSFile::ReadKey: '%s'", key);
    CheckStatus(buffer);

    return 0;
  }
}

void CFITSFile::ReadImage()
{
  int xw, yw;

  if( ReadKey("NAXIS1", tint, &xw) == 0 ||
      ReadKey("NAXIS2", tint, &yw) == 0 ) {
    fprintf(stderr, "*   CFITSFile::ReadImage(): no image data found\n");
    exit(-1);
  }

  cf_comment();
  printf(" Reading image, size %i x %i\n", xw, yw);

  int isnull;
  CFloatType anull = CNullValue;
  m_image.Resize(xw, yw);
  fits_read_img(m_file, CFITSFloatType, 1, m_image.GetNoPixels(),
		&anull, m_image.GetImageBuffer(),
		&isnull, &m_status);
  CheckStatus("Reading image");

  m_posnImage.ReadFITSHeader(*this);
}

void CFITSFile::ReadImageInclNull(CFloatType nullval)
{
  int xw, yw;

  if( ReadKey("NAXIS1", tint, &xw) == 0 ||
      ReadKey("NAXIS2", tint, &yw) == 0 ) {
    fprintf(stderr, "*   CFITSFile::ReadImage(): no image data found\n");
    exit(-1);
  }

  cf_comment();
  printf(" Reading image, size %i x %i\n", xw, yw);

  int isnull;
  CFloatType anull = nullval;
  m_image.Resize(xw, yw);
  fits_read_img(m_file, CFITSFloatType, 1, m_image.GetNoPixels(),
		&anull, m_image.GetImageBuffer(),
		&isnull, &m_status);
  CheckStatus("Reading image");

  m_posnImage.ReadFITSHeader(*this);
}

void CFITSFile::WriteImage()
{
  int xw = m_image.GetXW();
  int yw = m_image.GetYW();
  int bp = CFITSDataFormat;

  cf_comment();
  printf(" Writing image, size %i x %i\n", xw, yw);

  if(m_previousWritten) {
    UpdateKey("NAXIS1", tint, &xw, "X Width");
    UpdateKey("NAXIS2", tint, &yw, "Y Width");
    UpdateKey("BITPIX", tint, &bp, "Number of bits per data pixel"); 
  } else {
    long axes[2];
    axes[0] = xw; axes[1] = yw;
    fits_create_img(m_file, bp, 2, axes, &m_status);
    CheckStatus("Writing image header");
  }

  //  CFloatType anull = CNullValue;
  //  fits_write_imgnull(m_file, CFITSFloatType, 1, m_image.GetNoPixels(),
  //		     m_image.GetImageBuffer(), &anull, &m_status);

  fits_write_img(m_file, CFITSFloatType, 1, m_image.GetNoPixels(),
		 m_image.GetImageBuffer(), &m_status);
  CheckStatus("Writing image");

  m_posnImage.WriteFITSHeader(*this);

  m_previousWritten = 1;
}

void CFITSFile::WriteImageInclNull(CFloatType nullval)
{
  int xw = m_image.GetXW();
  int yw = m_image.GetYW();
  int bp = CFITSDataFormat;

  cf_comment();
  printf(" Writing image, size %i x %i\n", xw, yw);

  if(m_previousWritten) {
    UpdateKey("NAXIS1", tint, &xw, "X Width");
    UpdateKey("NAXIS2", tint, &yw, "Y Width");
    UpdateKey("BITPIX", tint, &bp, "Number of bits per data pixel"); 
  } else {
    long axes[2];
    axes[0] = xw; axes[1] = yw;
    fits_create_img(m_file, bp, 2, axes, &m_status);
    CheckStatus("Writing image header");
  }

  CFloatType anull = nullval;
  fits_write_imgnull(m_file, CFITSFloatType, 1, m_image.GetNoPixels(),
		     m_image.GetImageBuffer(), &anull, &m_status);

  CheckStatus("Writing image");

  m_posnImage.WriteFITSHeader(*this);

  m_previousWritten = 1;
}

void CFITSFile::SetImage(const CFITSImage &copy)
{
  m_image = copy;
}

CFITSImage &CFITSFile::GetImage()
{
  return m_image;
}

void CFITSFile::SetPosn(const CFITSPosn &copy)
{
  m_posnImage = copy;
}

CFITSPosn &CFITSFile::GetPosn()
{
  return m_posnImage;
}

void CFITSFile::WriteHistory(const char *hist)
{
  fits_write_history(m_file, hist, &m_status);
  CheckStatus("Adding history to CHU");
}

void CFITSFile::UpdateChecksum()
{
  fits_write_chksum(m_file, &m_status);
  CheckStatus("Updated checksum");
}

