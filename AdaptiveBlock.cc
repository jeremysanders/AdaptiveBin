// AdaptiveBlock - intelligent binning algorithm
// Jeremy Sanders 2000

// Idea:
//  Make image a power of two
//  Those pixels > sig/noise, keep
//  Otherwise bin 2x2, repeat, 4x4, etc...

#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <list>
#include <string>
#include <strstream>
#include <fstream>
#include "../FITSImage.h"
#include "../FITSFile.h"
#include "Coord.hh"
#include "SigCalc.hh"

namespace AdaptiveBin {

  const string progVersion = "0.2.2";
  const string progDate = __DATE__ ", " __TIME__;

  class CBlocker {
  public: // public interface
    CBlocker(CSigCalc *sigcalc, bool alternatealgoritm=false);
    ~CBlocker();
    void Block(double fracterr, int ivalue,
	       CFITSImage *outImage, CFITSImage *errImage,
	       CFITSImage *pixelImage);

  private: // private methods
    void binPass(int factor, int ivalue);          // actual algorithm
    void binPassMoveBlock(int factor, int ivalue); // alternative algorithm

  private: // private data
    CSigCalc *m_sigcalc;
    double m_fracterr;
    int m_side;
    const int m_realxw, m_realyw;

    CFITSImage *m_outImage;
    CFITSImage *m_errImage;
    CFITSImage *m_pixelImage;

    int m_pixelno;

    bool m_useAlternateAlgorithm;
  };

}

using namespace AdaptiveBin;
using std::fabs;
using std::sqrt;

//////////////////////////////////////////////////////////////

CBlocker::CBlocker(CSigCalc *sigcalc, bool alternatealgorithm) :
  m_sigcalc( sigcalc ),
  m_realxw( m_sigcalc->GetXW() ),
  m_realyw( m_sigcalc->GetYW() ),
  m_useAlternateAlgorithm(alternatealgorithm)
{
  m_fracterr = 0.0;

  // find largest power of two for edge
  int xw = 1;
  while( xw < m_realxw )
    xw *= 2;
  int yw = 1;
  while( yw < m_realyw )
    yw *= 2;

  m_side = std::max(xw, yw);

  // allocate output space
  m_outImage = new CFITSImage(m_realxw, m_realyw);
  m_errImage = new CFITSImage(m_realxw, m_realyw);
  m_pixelImage = new CFITSImage(m_realxw, m_realyw);
}

CBlocker::~CBlocker()
{
  delete m_outImage;
  delete m_errImage;
  delete m_pixelImage;
}

//////////////////////////////////////////////////////////////

// my binning algorithm
// pass over in powers of two, eliminating binned pixels
// that are below a threshold error

void CBlocker::binPass(int factor, int ivalue)
{
  cout << "    Bin pass " << factor << endl;

  int binnedside = m_side / factor;
  bool finalpass = (binnedside == 1);

  // loop over binned pixels
  for(int x = 0; x < binnedside; x++)
    for(int y = 0; y < binnedside; y++) {

      CSigCalc::cCoordList coordlist;

      // loop over displacement
      for(int xi = 0; xi < factor; xi++)
	for(int yi = 0; yi < factor; yi++) {
	  CCoord pt(x*factor+xi, y*factor+yi);

	  if( pt.m_x < m_realxw && pt.m_y < m_realyw &&
	      ! m_sigcalc->IsMasked(pt) )
	    coordlist.push_back( pt );
	}

      if( ! coordlist.empty() ) {
	double value = m_sigcalc->GetValue(coordlist, ivalue);
	double error = m_sigcalc->GetError(coordlist, 0);
	
	if( error < m_fracterr || finalpass ) { 
	  CSigCalc::cCoordList::const_iterator p = coordlist.begin();
	  double err = m_sigcalc->GetError(coordlist, ivalue);
	  while( p != coordlist.end() ) {
	    m_outImage->SetPixel( p->m_x, p->m_y, value );
	    m_errImage->SetPixel( p->m_x, p->m_y, err );
	    m_pixelImage->SetPixel( p->m_x, p->m_y, m_pixelno );
	    m_sigcalc->Mask(*p);
	    p++;
	  }
	  m_pixelno++;
	}
	
      } // coordlist != empty
      
    } // over points

}

// an alternative non-gridded algorithm
// very slow

// repeats lots of times to find the most accurate pixels

void CBlocker::binPassMoveBlock(int factor, int ivalue)
{
  cout << "    Bin pass " << factor << '\t';
  cout.flush();

  int binnedside = m_side / factor;
  bool finalpass = (binnedside == 1);

  CSigCalc::cCoordVList pixels;

  // loop over binned pixels

  while(true) {
    double minerr = 1e5;
    CSigCalc::cCoordList minerrblock;

    for(int xa = m_realxw-1; xa >= 0; xa--) {
      for(int ya = m_realyw-1; ya >= 0; ya--) {

	// get block coords
	CSigCalc::cCoordList block;
	for(int xi = factor-1; xi >= 0; xi--)
	  for(int yi = factor-1; yi >= 0; yi--) {
	    CCoord pt(xa+xi, ya+yi);

	    if(pt.m_x < m_realxw && pt.m_y < m_realyw &&
	       ! m_sigcalc->IsMasked(pt))
	      block.push_back(pt);
	  }

	// does block have minimum error
	if( ! block.empty() ) {
	  double err = m_sigcalc->GetError(block, 0);
	  if(err < minerr) {
	    minerr = err;
	    minerrblock = block;
	  }
	}
      }
    }

    if( (minerr > m_fracterr && !finalpass) || minerrblock.empty() )
      break;

    {
      CSigCalc::cCoordList::const_iterator p = minerrblock.begin();
      double val = m_sigcalc->GetValue(minerrblock, ivalue);
      double err = m_sigcalc->GetError(minerrblock, ivalue);
      while( p != minerrblock.end() ) {
	m_outImage->SetPixel( p->m_x, p->m_y, val );
	m_errImage->SetPixel( p->m_x, p->m_y, err );
	m_pixelImage->SetPixel( p->m_x, p->m_y, m_pixelno );
	m_sigcalc->Mask(*p);
	p++;
      }
      m_pixelno++;
    }

    cout << "*"; cout.flush();
  }
  cout << endl;
}

//////////////////////////////////////////////////////////////

void CBlocker::Block(double fracterr, int ivalue,
		     CFITSImage *outImage,
		     CFITSImage *errImage,
		     CFITSImage *pixelImage)
{
  m_fracterr = fracterr;
  m_outImage->SetAll(-1.0);
  m_errImage->SetAll(-1.0);
  m_pixelno = 0;

  for(int factor = 1; factor <= m_side; factor *= 2) {
    if( m_useAlternateAlgorithm )
      binPassMoveBlock(factor, ivalue);
    else
      binPass(factor, ivalue);
  }

  *outImage = *m_outImage;
  *errImage = *m_errImage;
  *pixelImage = *m_pixelImage;
}

//////////////////////////////////////////////////////////////

namespace AdaptiveBin {

  class CBlockerProgram
  {
  public:
    CBlockerProgram(int argc, char *argv[]);
    ~CBlockerProgram();
    
    void Run();
    
  private:
    void errorMessage();

    void loadImages();
    void doCalculation();
    void writeImages();

  private:
    enum cCalcType { count, colour2, colour3, colour4,
		     countmb, colour2mb, colour3mb,
		     colour4mb};
    static const int noCalcNames=8;
    
    std::list<std::string> m_fileParams;
    std::string m_outfilename, m_outerrfilename;
    cCalcType m_type;
    CSigCalc *m_sigcalc;
    double m_fracterr;

    CFITSImage m_outImage;
    CFITSImage m_outImageErr;
    CFITSImage m_pixelImage;
    CFITSPosn m_outPosn;

    bool m_useAlternateAlgorithm;
    int m_value;
  };

}

CBlockerProgram::CBlockerProgram(int argc, char *argv[])
  : m_sigcalc(0),
    m_useAlternateAlgorithm(false)
{
  if( argc < 7 )
    errorMessage();

  const char * const calcNames[noCalcNames] = {"count", "colour2",
					       "colour3", "colour4",
					       "countmb",
					       "colour2mb", "colour3mb"
					       "colour4mb"};
  const unsigned int noArgs[noCalcNames] = { 1, 2, 3, 4, 1, 2, 3, 4 };
  const int maxVal[noCalcNames] = { 0, 2, 6, 8, 0, 2, 6, 8 };

  string type = argv[1];

  m_outfilename = argv[2];
  m_outerrfilename = argv[3];

  {
    istrstream a4(argv[4]);
    if( ! (a4 >> m_fracterr) )
      errorMessage();
    istrstream a5(argv[5]);
    if( ! (a5 >> m_value) )
      errorMessage();
  }


  for(int i=6; i < argc; i++) {
    m_fileParams.push_back(argv[i]);
  }

  int it = 0;
  while( it < noCalcNames && type != calcNames[it] )
    it++;

  if( it == noCalcNames )
    errorMessage();

  if( m_value > maxVal[it] || m_value < 0 ) {
    clog << "* Val is not in range 0 - " << maxVal[it]
	 << "\n\n";
    errorMessage();
  }

  if( m_fileParams.size() != noArgs[it] )
    errorMessage();

  m_type = (cCalcType) it;
}

CBlockerProgram::~CBlockerProgram()
{
  if( m_sigcalc != 0 )
    delete m_sigcalc;
}

void CBlockerProgram::loadImages()
{
  // load in images
  std::list<std::string>::const_iterator x = m_fileParams.begin();
  switch( m_type ) {
  case count:
  case countmb:
    {
      CFITSFile file( x->c_str(), CFITSFile::existingro );
      m_outPosn = file.GetPosn();
      m_sigcalc = new CCountSig( file.GetImage() );
    }
    break;
  case colour2:
  case colour2mb:
    {
      CFITSFile filea( x->c_str(), CFITSFile::existingro );
      m_outPosn = filea.GetPosn();
      ++x;
      CFITSFile fileb( x->c_str(), CFITSFile::existingro );
      m_sigcalc = new CRatioSig( filea.GetImage(),
				 fileb.GetImage() );
    }
    break;
  case colour3:
  case colour3mb:
    {
      CFITSFile filea( x->c_str(), CFITSFile::existingro );
      m_outPosn = filea.GetPosn();
      ++x;
      CFITSFile fileb( x->c_str(), CFITSFile::existingro );
      ++x;
      CFITSFile filec( x->c_str(), CFITSFile::existingro );

      m_sigcalc = new CRatioSig3( filea.GetImage(),
				  fileb.GetImage(),
				  filec.GetImage() );
    }
    break;
  case colour4:
  case colour4mb:
    {
      CFITSFile filea( x->c_str(), CFITSFile::existingro );
      m_outPosn = filea.GetPosn();
      ++x;
      CFITSFile fileb( x->c_str(), CFITSFile::existingro );
      ++x;
      CFITSFile filec( x->c_str(), CFITSFile::existingro );
      ++x;
      CFITSFile filed( x->c_str(), CFITSFile::existingro );

      m_sigcalc = new CRatioSig4( filea.GetImage(),
				  fileb.GetImage(),
				  filec.GetImage(),
				  filed.GetImage() );
    }
    break;
  }


  if( m_type==countmb || m_type==colour2mb || m_type==colour3mb )
    m_useAlternateAlgorithm = true;
}

void CBlockerProgram::doCalculation()
{
  // do the calculations

  CBlocker blocker(m_sigcalc, m_useAlternateAlgorithm);
  blocker.Block(m_fracterr, m_value, &m_outImage, &m_outImageErr,
		&m_pixelImage);
}

void CBlockerProgram::writeImages()
{
  // write the output files
  {
    CFITSFile outfile(m_outfilename.c_str(), CFITSFile::create);
    outfile.SetImage(m_outImage);
    outfile.SetPosn(m_outPosn);
    outfile.WriteImage();
  }

  {
    CFITSFile outerrfile(m_outerrfilename.c_str(), CFITSFile::create);
    outerrfile.SetImage(m_outImageErr);
    outerrfile.SetPosn(m_outPosn);
    outerrfile.WriteImage();
  }

  {
    CFITSFile outpixelfile("pixel.fits", CFITSFile::create);
    outpixelfile.SetImage(m_pixelImage);
    outpixelfile.SetPosn(m_outPosn);
    outpixelfile.WriteImage();
  }
}

void CBlockerProgram::Run()
{
  loadImages();
  std::cout << ":  Fractional error setting = " << m_fracterr << endl;
  doCalculation();
  writeImages();
}

void CBlockerProgram::errorMessage()
{
  std::clog << "AdaptiveBlock (version " << progVersion
	    << ", " << progDate << "), Jeremy Sanders\n\n"
	    << "  Usage:  AdaptiveBlock type outf outf-err sig val"
	    << " file1 [file2] [file3] [...]\n\n"
	    << "type can be count (max val 0), colour2 (mv 2), colour3 (mv 6)\n"
	    << " colour4 (mv 8) (add mb for move-block [very slow!] method)\n"
	    << "sig is fractional error\n"
	    << "val is the value type to return (0 selects default)\n";
  exit(-1);
}

//////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

  CBlockerProgram program(argc, argv);
  program.Run();

  return 0;
}
