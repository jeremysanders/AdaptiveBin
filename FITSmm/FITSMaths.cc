// program for doing maths on FITS files
// reverse polish format
// This is _very_ old code and doesn't do this in a nice way
//  oh well.

// syntax: fits1.fits fits2.fits +
// or fits1.fits[1:100,3:200] fits2.fits[1:200:2,3:200] - (not implemented!)

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


#include <cmath>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <FITSImage.h>
#include <FITSFile.h>
#include <FITSPosn.h>

using std::cerr;
using std::cout;
using std::endl;
using std::swap;

// null value for comparisons
const double nullval = 1.320475e-37;

// bad implementation and should probably take someone else's
class CFITSStack {
public:
  CFITSStack();
  ~CFITSStack();

  void Push(const char *filename);
  void Push(const CFITSImage &image);
  void Pop(CFITSImage &retn);

  void Drop();
  void Swap();
  void Dup();
  void Add();
  void Sub();
  void Mul();
  void Div();
  void Ln();
  void Log();
  void Sqrt();
  void SignSqrt();
  void Exp();
  void Min();
  void Max();


  void Lt();
  void Gt();
  void LtRoll();
  void GtRoll();
  void TrimDown();
  void TrimUp();

  int GetNoItems();

  void GetPosn(CFITSPosn *posn);

public:
  class CBoundsException {};
  class CInvalidException {};

private:
  static const int cmaxitems = 1024;

  CFITSPosn m_posn;
  CFITSImage *m_items[cmaxitems];
  int m_noitems;
};

CFITSStack::CFITSStack()
{
  m_noitems = 0;
}

CFITSStack::~CFITSStack()
{
  for(int i = 0; i < m_noitems; i++)
    delete m_items[i];
}

void CFITSStack::GetPosn(CFITSPosn *posn)
{
  *posn = m_posn;
}

int CFITSStack::GetNoItems()
{
  return m_noitems;
}

void CFITSStack::Push(const char *filename)
{
  // copy filename part
  char fin[1024];
  int i=0;
  while( filename[i] != '[' && filename[i] != 0 ) {
    fin[i] = filename[i]; i++;
  }
  fin[i] = 0;

  if( fin[0] != '@' ) {  // normal file
    // try to open the image
    CFITSImage image;
    {
      CFITSFile file(fin, CFITSFile::existingro);
      image = file.GetImage();
      m_posn = file.GetPosn();
    }

    if( filename[i] == '[' ) {
      // implement
    }

    Push(image);
  } else { // const image
    // get number
    double num;
    if( sscanf(fin+1, "%lf", &num) != 1 ) {
      cerr << "*  Argument after @ not number - \""
	   << fin+1 << "\"\n";
      throw CInvalidException();
    }

    // get size of previous image
    if( m_noitems == 0 ) {
      cerr << "*  Cannot make constant entry first\n";
      throw CInvalidException();
    }
    CFITSImage image(m_items[m_noitems-1] -> GetXW(),
		     m_items[m_noitems-1] -> GetYW());
    for(int x = image.GetXW()-1; x >= 0; x--)
      for(int y = image.GetYW()-1; y >= 0; y--)
	image.SetPixel(x, y, num);
    Push(image);
  }

}

void CFITSStack::Push(const CFITSImage &image)
{
  m_noitems++;
  if(m_noitems >= cmaxitems) {
    cerr << "*  Cannot push - no room left on stack";
    throw CBoundsException();
  }

  m_items[m_noitems-1] = new CFITSImage( image );
}

void CFITSStack::Pop(CFITSImage &retn)
{
  if(m_noitems == 0) {
    cerr << "*  Cannot pop - no items on stack\n";
    throw CBoundsException();
  }

  retn = *m_items[m_noitems-1];
  delete m_items[m_noitems-1];
  m_noitems--;
}

void CFITSStack::Drop()
{
  if(m_noitems == 0) {
    cerr << "*  Cannot drop - no items on stack\n";
    throw CBoundsException();
  }

  CFITSImage retn;
  Pop(retn);
}

void CFITSStack::Swap()
{
  if(m_noitems < 2) {
    cerr << "*  Cannot swap, less than 2 items on stack\n";
    throw CBoundsException();
  }

  swap( m_items[m_noitems-1], m_items[m_noitems-2] );
}

void CFITSStack::Dup()
{
  if(m_noitems == 0) {
    cerr << "*  Cannot dup - no items on stack\n";
    throw CBoundsException();
  }

  m_noitems++;
  if(m_noitems >= cmaxitems) {
    cerr << "*  Cannot dup - no space left on stack\n";
    throw CBoundsException();
  }
  m_items[m_noitems-1] = new CFITSImage(*m_items[m_noitems-2]);
}

void CFITSStack::Add()
{
  CFITSImage im1, im2;
  Pop(im1);
  Pop(im2);
  im1 += im2;
  Push(im1);
}

void CFITSStack::Sub()
{
  CFITSImage im1, im2;
  Pop(im2);
  Pop(im1);
  im1 -= im2;
  Push(im1);
}

void CFITSStack::Mul()
{
  CFITSImage im1, im2;
  Pop(im1);
  Pop(im2);
  im1 *= im2;
  Push(im1);
}

void CFITSStack::Div()
{
  CFITSImage im1, im2;
  Pop(im2);
  Pop(im1);
  im1 /= im2;
  Push(im1);
}

void CFITSStack::Ln()
{
  CFITSImage image;
  Pop(image);

  for(int x = image.GetXW()-1; x >= 0; x--)
    for(int y = image.GetYW()-1; y >=0; y--) {
      double pix = image.GetPixel(x, y);
      if(pix > 0) image.SetPixel(x, y, log(pix));
      else {
	image.SetPixel(x, y, 0.0);
	cerr << "*  ln of negative value attempted\n";
      }
    }

  Push(image);
}

void CFITSStack::Log()
{
  CFITSImage image;
  Pop(image);

  unsigned noneg = 0;
  for(int x = image.GetXW()-1; x >= 0; x--)
    for(int y = image.GetYW()-1; y >=0; y--) {
      double pix = image.GetPixel(x, y);
      if(pix > 0) image.SetPixel(x, y, log10(pix));
      else {
	image.SetPixel(x, y, log10(fabs(pix)));
	noneg++;
      }
    }

  if(noneg != 0)
    {
      std::cerr << "* Warning: found " << noneg << " negative pixels\n";
    }

  Push(image);
}

void CFITSStack::Sqrt()
{
  CFITSImage image;
  Pop(image);

  unsigned noneg = 0;
  for(int x = image.GetXW()-1; x >= 0; x--)
    for(int y = image.GetYW()-1; y >=0; y--) {
      double pix = image.GetPixel(x, y);
      if(pix >= 0) image.SetPixel(x, y, sqrt(pix));
      else {
	image.SetPixel(x, y, sqrt(fabs(pix)));
	noneg++;
      }
    }

  if(noneg != 0)
    {
      std::cerr << "* Warning: found " << noneg << " negative pixels\n";
    }

  Push(image);
}

void CFITSStack::SignSqrt()
{
  CFITSImage image;
  Pop(image);

  for(int x = image.GetXW()-1; x >= 0; x--)
    for(int y = image.GetYW()-1; y >=0; y--) {
      const double pix = image.GetPixel(x, y);
      image.SetPixel( x, y, copysign( sqrt(fabs(pix)), pix ) );
    }

  Push(image);
}

void CFITSStack::Exp()
{
  CFITSImage image;
  Pop(image);

  for(int x = image.GetXW()-1; x >= 0; x--)
    for(int y = image.GetYW()-1; y >=0; y--)
      image.SetPixel(x, y, exp(image.GetPixel(x, y)) );

  Push(image);
}

void CFITSStack::Lt()
{
  CFITSImage im1, im2;
  Pop(im1);
  Pop(im2);

  int count = 0;
  for(int x = im1.GetXW()-1; x >= 0; x--)
    for(int y = im2.GetYW()-1; y >=0; y--) {
      const double p1 = im1.GetPixel(x, y);
      const double p2 = im2.GetPixel(x, y);

      if( p2 >= p1 )
	im2.SetPixel(x, y, nullval);
      else
	++count;
    }

  Push(im2);

  cout << count << " pixels satisfy lt criterion\n";
}

void CFITSStack::Min()
{
  CFITSImage im1, im2;
  Pop(im1);
  Pop(im2);

  for(int x = im1.GetXW()-1; x >= 0; x--)
    for(int y = im2.GetYW()-1; y >=0; y--) {
      const double p1 = im1.GetPixel(x, y);
      const double p2 = im2.GetPixel(x, y);

      if( p2 > p1 )
	im2.SetPixel(x, y, p1);
    }

  Push(im2);
}

void CFITSStack::Max()
{
  CFITSImage im1, im2;
  Pop(im1);
  Pop(im2);

  for(int x = im1.GetXW()-1; x >= 0; x--)
    for(int y = im2.GetYW()-1; y >=0; y--) {
      const double p1 = im1.GetPixel(x, y);
      const double p2 = im2.GetPixel(x, y);

      if( p2 < p1 )
	im2.SetPixel(x, y, p1);
    }

  Push(im2);
}

void CFITSStack::Gt()
{
  CFITSImage im1, im2;
  Pop(im1);
  Pop(im2);

  int count = 0;

  for(int x = im1.GetXW()-1; x >= 0; x--)
    for(int y = im2.GetYW()-1; y >=0; y--) {
      const double p1 = im1.GetPixel(x, y);
      const double p2 = im2.GetPixel(x, y);

      if( p2 <= p1 )
	im2.SetPixel(x, y, nullval);
      else
	++count;
    }

  Push(im2);

  cout << count << " pixels satisfy gt criterion\n";
}

void CFITSStack::LtRoll()
{
  CFITSImage im1, im2, im3;
  Pop(im1);
  Pop(im2);
  Pop(im3);

  int count = 0;
  for(int x = im1.GetXW()-1; x >= 0; x--)
    for(int y = im2.GetYW()-1; y >=0; y--) {
      const double p1 = im1.GetPixel(x, y);
      const double p2 = im2.GetPixel(x, y);

      if( p2 >= p1 )
	im3.SetPixel(x, y, nullval);
      else
	++count;
    }

  Push(im3);

  cout << count << " pixels satisfy lt criterion\n";
}

void CFITSStack::GtRoll()
{
  CFITSImage im1, im2, im3;
  Pop(im1);
  Pop(im2);
  Pop(im3);

  int count = 0;
  for(int x = im1.GetXW()-1; x >= 0; x--)
    for(int y = im2.GetYW()-1; y >=0; y--) {
      const double p1 = im1.GetPixel(x, y);
      const double p2 = im2.GetPixel(x, y);

      if( p2 <= p1 )
	im3.SetPixel(x, y, nullval);
      else
	++count;
    }

  Push(im3);

  cout << count << " pixels satisfy lt criterion\n";
}

void CFITSStack::TrimDown()
{
  CFITSImage im1, im2;
  Pop(im1);
  Pop(im2);

  int count = 0;
  for(int x = im1.GetXW()-1; x >= 0; x--)
    for(int y = im2.GetYW()-1; y >=0; y--) {
      const double p1 = im1.GetPixel(x, y);
      const double p2 = im2.GetPixel(x, y);

      if( p2 >= p1 )
	{
	  im2.SetPixel(x, y, p1);
	  ++count;
	}
    }

  Push(im2);

  cout << count << " pixels trimmed\n";
}

void CFITSStack::TrimUp()
{
  CFITSImage im1, im2;
  Pop(im1);
  Pop(im2);

  int count = 0;
  for(int x = im1.GetXW()-1; x >= 0; x--)
    for(int y = im2.GetYW()-1; y >=0; y--) {
      const double p1 = im1.GetPixel(x, y);
      const double p2 = im2.GetPixel(x, y);

      if( p2 <= p1 )
	{
	  im2.SetPixel(x, y, p1);
	  ++count;
	}
    }

  Push(im2);

  cout << count << " pixels trimmed\n";
}

//////////////////////////////////////////////////////////

class CFITSMathProgram
{
public:
  CFITSMathProgram();
  ~CFITSMathProgram();

  void ProcessArgs(int argc, char **argv);

private:
  void errormessage();
  void handlearg(const char *arg);

private:
  CFITSStack m_stack;
  
};

CFITSMathProgram::CFITSMathProgram()
{
  cout << "FITSMaths, Jeremy Sanders, 2000, 2001, 2002, 2003\n";
}

CFITSMathProgram::~CFITSMathProgram()
{
}

void CFITSMathProgram::ProcessArgs(int argc, char **argv)
{
  if(argc < 3)
    errormessage();

  for(int i=0; i<argc-2; i++)
    handlearg(argv[i+1]);

  CFITSImage outimage;
  m_stack.Pop(outimage);

  if(m_stack.GetNoItems() != 0)
    cerr << "*  Warning: Item(s) left on stack at end\n";

  CFITSPosn posn;
  m_stack.GetPosn(&posn);

  CFITSFile file(argv[argc-1], CFITSFile::create);
  file.SetImage(outimage);
  file.SetPosn(posn);
  file.WriteImageInclNull(nullval);
}

void CFITSMathProgram::errormessage()
{
  cerr <<
    "FITSMaths Version 1.14\n"
    "Usage:\n"
    " FITSMaths infile [infile|operator] ... outfile\n"
    "infile is the name of a file or a constant (format @const)\n"
    "operator can be +,-,/,*,dup,ln,log,sqrt,exp,drop,swap,lt,gt,\n"
    "                ltroll,gtroll,trimdown,trimup,signsqrt,min,max\n";
  exit(-1);
}

void CFITSMathProgram::handlearg(const char *arg)
{
  const int noops = 20;
  const char *ops[]={"+","-","/","*","dup","ln","log","sqrt",
		     "exp", "drop", "swap", "lt", "gt", "ltroll", "gtroll",
		     "trimdown", "trimup",
                     "signsqrt", "min", "max" };

  int opsel = -1;
  for(int i=0; i<noops; i++) {
    if( strcmp(ops[i], arg) == 0 )
      opsel = i;
  }

  if( opsel < 0 ) {
    try {
      m_stack.Push(arg);
    }
    catch( CFITSStack::CInvalidException x ) {
      exit(-1);
    }
    return;
  }

  try {
    switch( opsel ) {
    case 0: m_stack.Add(); break;
    case 1: m_stack.Sub(); break;
    case 2: m_stack.Div(); break;
    case 3: m_stack.Mul(); break;
    case 4: m_stack.Dup(); break;
    case 5: m_stack.Ln();  break;
    case 6: m_stack.Log(); break;
    case 7: m_stack.Sqrt(); break;
    case 8: m_stack.Exp(); break;
    case 9: m_stack.Drop(); break;
    case 10: m_stack.Swap(); break;
    case 11: m_stack.Lt(); break;
    case 12: m_stack.Gt(); break;
    case 13: m_stack.LtRoll(); break;
    case 14: m_stack.GtRoll(); break;
    case 15: m_stack.TrimDown(); break;
    case 16: m_stack.TrimUp(); break;
    case 17: m_stack.SignSqrt(); break;
    case 18: m_stack.Min(); break;
    case 19: m_stack.Max(); break;
    }
  }
  catch( CFITSStack::CBoundsException x ) {
    exit(-1);
  }

}

/////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  CFITSMathProgram prog;
  prog.ProcessArgs(argc, argv);

  return 0;
}
