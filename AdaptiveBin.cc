//      Adaptive Binning Program
//      Main program module
//      Described in Sanders and Fabian (submitted)
//      Routines for adaptively binning data
//      Copyright (C) 2000, 2001 Jeremy Sanders
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

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <sstream>
#include <cassert>

#include <parammm/parammm.hh>
#include <FITSFile.h>

#include "binmodule.hh"
#include "version.hh"

using std::string;
using std::sort;
using std::vector;
using std::cout;
using std::cerr;
using std::clog;
using std::endl;
using std::ostringstream;

namespace AdaptiveBin {

  class binval
  {
  public:
    binval(const pixlist &pl, double val);
    bool operator< (const binval &cmp) const;
    const pixlist& pixels() const { return m_pl; }
    double val() const { return m_val; }
  private:
    pixlist m_pl;
    double m_val;
  };

  binval::binval(const pixlist &pl, double val)
  {
    m_pl = pl;
    m_val = val;
  }

  bool binval::operator<(const binval &cmp) const
  {
    return m_val < cmp.m_val;
  }

  typedef std::vector<binval> binval_list;

  // binner is the class which does the binning
  // it uses the binmodule to find the errors on each pixel

  class binner
  {
  public:
    binner(binmodule *bm, double threshold,
	   int subpixposn, bool contig_check);
    // bm is binning module (by which mode we're binning)
    // threshold is fractional error threshold
    // subpixposn is number of subbins/bin
    // if contig_check is true, then we only bin 'contiguous' regions

    ~binner();

    void bin(CFITSImage *out_image,
	     CFITSImage *error_image,
	     CFITSImage *binmap_image);

    void set_mask_image(const CFITSImage &mask,
			bool invert_mask = false);

  private:
    //    void binpass(int pixsize, bool finalpass);
    void apply_mask();
    void pass_bins_and_sort(int pixsize, bool finalpass);
    void check_noncontiguous_bin(const pixlist &binpixels,
				 bool finalpass,
				 binval_list *binlist);
    
  private:
    int m_latest_bin_no;                     // keep a count of the outputted bins

    double m_threshold;                      // threshold error value
    int m_subbinposn;                        // sub-bin positioning value
    bool m_contig_check;                     // only allow contig regions
    binmodule *m_binmod;                     // module to do the binning
    CFITSImage m_out_image;                  // output binned image
    CFITSImage m_error_image;                // output error image
    CFITSImage m_output_binmap_image;        // output binmap
    CFITSImage m_mask_image;                 // mask to use in binning
  };

  binner::binner(binmodule *bm, double threshold, int subpixposn,
		 bool contig_check)
    : m_threshold(threshold),
      m_subbinposn(subpixposn),
      m_contig_check(contig_check),
      m_binmod(bm),
      m_out_image(bm->xw(), bm->yw()),
      m_error_image(bm->xw(), bm->yw()),
      m_output_binmap_image(bm->xw(), bm->yw()),
      m_mask_image(bm->xw(), bm->yw())
  {
  }

  binner::~binner()
  {
  }

  void binner::set_mask_image(const CFITSImage &mask, bool invert_mask)
  {
    // check mask is the same size as the image
    assert( mask.GetXW() == m_mask_image.GetXW() &&
	    mask.GetYW() == m_mask_image.GetYW() );

    m_mask_image = mask;

    if(invert_mask)
      {
	const unsigned xw = mask.GetXW(), yw = mask.GetYW();

	for(unsigned y=0; y<yw; ++y)
	  for(unsigned x=0; x<xw; ++x)
	    {
	      double v = m_mask_image.GetPixel(x, y);
	      if( v < 1e-10 )
		v = 1;
	      else
		v = 0;
	      m_mask_image.SetPixel(x, y, v);
	    }
      }
  }

  void binner::bin(CFITSImage *out_image,
		   CFITSImage *error_image,
		   CFITSImage *binmap_image)
  {
    m_output_binmap_image.SetAll(-2.);  // -2 specifies unbinned, -1 specifies masked
    m_out_image.SetAll(-1.);
    m_error_image.SetAll(-1.);

    m_latest_bin_no = 0;
    apply_mask();

    // do passes over factor of 2
    int pass;
    for(pass=1; pass<m_binmod->xw() || pass<m_binmod->yw(); pass *= 2)
      pass_bins_and_sort(pass, false);

    // final pass
    pass_bins_and_sort(pass, true);

    // return values
    *out_image = m_out_image;
    *error_image = m_error_image;
    *binmap_image = m_output_binmap_image;
  }

  void binner::apply_mask()
  {
    for(int y = 0; y < m_binmod->yw(); ++y)
      for(int x = 0; x < m_binmod->xw(); ++x)
	{
	  if( m_mask_image.GetPixel(x, y) > 0. )
	    m_output_binmap_image.SetPixel(x, y, -1.);
	}
  }

  void binner::check_noncontiguous_bin(const pixlist &bin,
				       bool finalpass,
				       binval_list *binlist)
  {
    // this routine starts from a single pixels, and tries
    // to build a contiguous chunk of pixels from that, until it
    // can't go any further. corners touching are allowed

    // it then tries to build another until all the pixels
    // are used up

    const int nopixels = bin.size();

    vector<bool> pixel_used;
    for(int i=nopixels-1; i>=0; i--)
      pixel_used.push_back(false);

    // Surely this is far too complex - rewrite at some point for v2
    for(;;)
      {
	
	// find next available pixel
	int next_pixel = 0;
	while(next_pixel<nopixels && pixel_used[next_pixel])
	  next_pixel++;
	
	if(next_pixel == nopixels)
	  break; // no more chunks available
	
	// make an image to hold processed pixels
	const int xw = m_out_image.GetXW() + 1;
	const int yw = m_out_image.GetYW() + 1;
	bool *bin_image = new bool[xw*yw];
	for(int i=0; i<xw*yw; i++)
	  bin_image[i] = false;
	
	// try and build 'chunk'
	pixlist chunk;
	const pixel &first = bin[next_pixel];
	chunk.push_back(first);
	bin_image[ first.x() + xw*first.y() ] = true;
	
	for(;;) {  // loop, adding contig pixels, until there are no more
	  bool added_any = false;
	  
	  // loop over pixels in input bin
	  for(int pix=0; pix<nopixels; pix++) {
	    const pixel &curr = bin[pix];
	    
	    if( ! pixel_used[pix] ) {  // if the pixel is left to process
	      
	      // go over already added pixels, and see whether this pixel
	      // is contiguous
	      
	      bool is_contig = false;
	      
	      for(int y=-1; y<=1; y++)
		for(int x=-1; x<=1; x++) {
		  const int xp = curr.x() + x;
		  const int yp = curr.y() + y;
		  
		  if( xp >= 0 && yp >= 0 ) {
		    if( bin_image[xp+yp*xw] )
		      is_contig = true;
		  }
		} 
	      
	      if(is_contig) {
		added_any = true;
		pixel_used[pix] = true;
		chunk.push_back(curr);
		bin_image[curr.x()+curr.y()*xw] = true;
	      }
	      
	    } // pixel_used
	    
	  } // loop over pixels
	  
	  if(!added_any) break;
	} // repeating loop over pixels
	
	// see whether error on chunk is less than threshold
	const double error = m_binmod -> fracerror(chunk, true);
	if(error <= m_threshold || finalpass) {
	  binval p(chunk, error);
	  binlist->push_back(p);
	}

      delete[] bin_image;

    }  // try another chunk

  }

  // sorting binning version of binpass
  void binner::pass_bins_and_sort(int size, bool finalpass)
  {
    cout << "Pass " << size << endl;
    // code to allow bins to start on sub-bin boundries
    int nx, ny, ns;
    if( size < m_subbinposn ) {
      nx = m_binmod->xw() + 1;
      ny = m_binmod->yw() + 1;
      ns = 1;
    } else {
      ns = size / m_subbinposn;
      nx = m_binmod->xw() / ns + 1;
      ny = m_binmod->yw() / ns + 1;
    }

    // make an array of bins with error less than threshold
    binval_list binslist;
    // iterate over subbins
    for(int x=0; x<nx; x++)
      for(int y=0; y<ny; y++) {
	  
	// make bin with size size x size
	pixlist pixels;
	for(int sx=0; sx < size; sx++)
	  for(int sy=0; sy < size; sy++) {
	    const int tx = x*ns+sx;
	    const int ty = y*ns+sy;
	    if(tx < m_binmod->xw() && ty < m_binmod->yw() ) {

	      if( m_output_binmap_image.GetPixel(tx, ty) < -1. ) {
		pixels.push_back( pixel(tx, ty) );
	      }

	    }
	  }
	// end bin
	
	// are there any pixels in bin?
	if( pixels.size() > 0 ) {
	  // is binning error < threshold
	  const double error = m_binmod -> fracerror(pixels, true);
	  if(error <= m_threshold || finalpass) {

	    // if we need pixels to be contiguous, check for it
	    // otherwise just do it
	    if(m_contig_check) {
	      check_noncontiguous_bin(pixels, finalpass,
				      &binslist);
	    } else {
	      binval p(pixels, error);
	      binslist.push_back(p);
	    }

	  }
	}
    
      } // pixels

    // sort bins into error order
    // (not needed unless subbinning on)
    if(m_subbinposn != 1)
      sort(binslist.begin(), binslist.end());

    // set pixels which haven't already been set
    // with the values
    const int mi=binslist.size();
    for(int i=0; i<mi; i++) {
      // set output pixels to correct thing
      const pixlist &minerrpixel = binslist[i].pixels();

      bool spoilt = false;
      for(int j=minerrpixel.size()-1; j>=0 && !spoilt; j--) {
	const int x = minerrpixel[j].x(), y = minerrpixel[j].y();	
	if( ! (m_output_binmap_image.GetPixel(x, y) < 0.) )
	  spoilt = true;
      }
      if(spoilt) continue;

      const double val = m_binmod -> value(minerrpixel);
      const double outerror = m_binmod -> fracerror(minerrpixel, false);
      for(int j=minerrpixel.size()-1; j>=0; j--) {
	const int x = minerrpixel[j].x(), y = minerrpixel[j].y();
	m_output_binmap_image.SetPixel(x, y, m_latest_bin_no);
	m_out_image.SetPixel(x, y, val);
	m_error_image.SetPixel(x, y, outerror);
      }
      m_latest_bin_no ++;
    }

  } // fn

} // namespace


///////////////////////////////////////////////////////////////////////
// Program class

class prog
{
public:
  prog(int argc, char **argv);
  ~prog();
  void run();

private:
  void add_history_list(CFITSFile *file);

public:
  AdaptiveBin::binmodule *m_binmod;
  double m_threshold;    // threshold value
  string m_out_fname;    // output binned filename
  string m_err_fname;    // output error filename
  string m_binmap_fname; // output binmap filename
  string m_mask_fname;   // filename of the mask to use (optional)
  string m_value;        // quantity to bin
  int m_sub_bin;         // sub-binning value
  bool m_contig;         // only allow contiguous regions
  bool m_verbose;        // display verbose information
  bool m_invert_mask;    // invert 0 and 1 in mask

  vector<string> m_history_list;
};

prog::prog(int argc, char **argv)
  : m_binmod(0),
    m_threshold(0.1),
    m_out_fname("adbin_out.fits"),
    m_err_fname("adbin_err.fits"),
    m_binmap_fname("adbin_binmap.fits"),
    m_value("count(0)"),
    m_sub_bin(1),
    m_contig(false),
    m_verbose(false),
    m_invert_mask(false)
{
  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch("out", 'o',
				      parammm::pstring_opt(&m_out_fname),
				      "set out file (def adbin_out.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("error", 'e',
				      parammm::pstring_opt(&m_err_fname),
				      "set error out file (def adbin_err.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("binmap", 'n',
				      parammm::pstring_opt(&m_binmap_fname),
				      "set binmap out file (def adbin_binmap.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("mask", 'm',
				      parammm::pstring_opt(&m_mask_fname),
				      "set input mask filename (optional)",
				      "FILE"));
  params.add_switch( parammm::pswitch("pixel", 'p',
				      parammm::pstring_opt(&m_binmap_fname),
				      "set pixel out file (DISCOURAGED)",
				      "FILE"));
  params.add_switch( parammm::pswitch("threshold", 't',
				      parammm::pdouble_opt(&m_threshold),
				      "set threshold fraction (def 10%)",
				      "VAL"));
  params.add_switch( parammm::pswitch("value", 'v',
				      parammm::pstring_opt(&m_value),
				      "set output value (eg count(0), "
				      "ratio(1,2))", "STR") );
  params.add_switch( parammm::pswitch("subpix", 's',
				      parammm::pint_opt(&m_sub_bin),
				      "set subpixel positioning "
				      "divisior (def. 1)",
				      "INT"));
  params.add_switch( parammm::pswitch("contig", 'c',
				      parammm::pbool_noopt(&m_contig),
				      "only allow contiguous regions",
				      ""));
  params.add_switch( parammm::pswitch("invertmask", 0,
				      parammm::pbool_noopt(&m_invert_mask),
				      "invert input mask image",
				      ""));
  params.add_switch( parammm::pswitch("verbose", 0,
				      parammm::pbool_noopt(&m_verbose),
				      "display more information",
				      ""));


  params.set_autohelp("Usage: AdaptiveBin [OPTIONS] file bg=count...\n"
		      "Adaptively bins a set of images\n"
		      "Written by Jeremy Sanders, 2000, 2001.",
		      "Report bugs to <jss@ast.cam.ac.uk>");
  params.enable_autohelp();
  params.enable_autoversion(c_adbin_version,
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();

  params.interpret_and_catch();

  if(params.args().size() < 1)
    params.show_autohelp();

  try {
    // select external if specified
    const string first8(m_value, 0, 8);
    if( first8 == "external" )
      m_binmod = new AdaptiveBin::external_binmodule(params.args());
    else
      m_binmod = new AdaptiveBin::ratio_binmodule(params.args());
    m_binmod->selectvalue(m_value);
  }
  catch(AdaptiveBin::invalidargs_exception e) {
    clog << "Invalid files listed\n\n";
    params.show_autohelp();
  }
  catch(AdaptiveBin::invalidvalue_exception e) {
    clog << "Invalid output value\n\n";
    params.show_autohelp();
  }

  cout << "Using value "
       << m_binmod->get_value_descr() << endl;


  {
    // stuff to write into history in fits output files
    m_history_list.push_back(string("file created by AdaptiveBin v. ")
			     + c_adbin_version);
    
    for(int i = 0; i<int(params.args().size()); ++i) {
      ostringstream o;
      o << "arg " << i << ": " << params.args()[i] << '\0';
      m_history_list.push_back( o.str() );
    }
    
    m_history_list.push_back( string("output image: ") + m_out_fname );
    m_history_list.push_back( string("error map: ") + m_err_fname );
    m_history_list.push_back( string("bin map: ") + m_binmap_fname );

    m_history_list.push_back( string("mask: ") + m_mask_fname );
    m_history_list.push_back( string("value: ") + m_binmod->get_value_descr() );
    m_history_list.push_back( string("contig: ") + (m_contig ? "true" : "false") );
    
    {
      ostringstream o;
      o << "threshold: " << m_threshold << '\0';
      m_history_list.push_back( o.str() );
    }{
      ostringstream o;
      o << "subpix: " << m_sub_bin << '\0';
      m_history_list.push_back( o.str() );
    }
  } // end history comments

  // write history to screen if verbose option is on
  if(m_verbose) {
    cout << "\nHeader lines written to output files:\n";
    for(int i=0; i<int(m_history_list.size()); ++i)
      cout << m_history_list[i] << endl;
    cout << endl;
  }
}

prog::~prog()
{
  if(m_binmod != 0)
    delete m_binmod;
}

void prog::add_history_list(CFITSFile *file)
{
  const int no = m_history_list.size();

  for(int i=0; i<no; ++i) {
    const string line = "adbin: " + m_history_list[i];
    file -> WriteHistory(line.c_str());
  }
}

void prog::run()
{
  CFITSImage out, err, pixel;
  AdaptiveBin::binner b(m_binmod, m_threshold, m_sub_bin,
			m_contig);

  if( ! m_mask_fname.empty() ) {
    CFITSFile mask_file(m_mask_fname.c_str(), CFITSFile::existingro);
    b.set_mask_image(mask_file.GetImage(), m_invert_mask);
  }

  b.bin(&out, &err, &pixel);

  CFITSPosn posn;
  m_binmod -> getposn(&posn);

  {
    CFITSFile outf(m_out_fname.c_str(), CFITSFile::create);
    outf.SetImage(out);
    outf.SetPosn(posn);
    outf.WriteImageInclNull(-1.); // ignore masked bins
    outf.WriteHistory("adbin: file is output image");
    add_history_list( &outf );
  }{
    CFITSFile outf(m_err_fname.c_str(), CFITSFile::create);
    outf.SetImage(err);
    outf.SetPosn(posn);
    outf.WriteImageInclNull(-1.); // ignore masked bins
    outf.WriteHistory("adbin: file is error map");
    add_history_list( &outf );
  }{
    CFITSFile outf(m_binmap_fname.c_str(), CFITSFile::create);
    outf.SetImage(pixel);
    outf.SetPosn(posn);
    outf.WriteImage();
    outf.WriteHistory("adbin: file is bin map");
    add_history_list( &outf );
  }

}

int main(int argc, char *argv[])
{
  prog program(argc, argv);
  program.run();
  return 0;
}
