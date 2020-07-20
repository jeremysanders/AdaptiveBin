//      Program to make bin-map from smoothed image

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

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include <parammm/parammm.hh>
#include <FITSFile.h>

using namespace std;

class point
{
public:
  point() : m_x(0), m_y(0) {}
  point(int x, int y) { m_x = x; m_y = y; }
  
  int& x() { return m_x; }
  int& y() { return m_y; }
  const int &x() const { return m_x; }
  const int &y() const { return m_y; }

private:
  int m_x, m_y;
};

typedef vector<point> point_vector;

////////////////////////////////////////////////////////////

class contour_prog
{
public:
  contour_prog();
  ~contour_prog();

  void run();

private:
  void read_contours(const string &filename);
  void read_smoothed_image(const string &filename);

  void make_contours();
  void check_contour_contig(int cont);

  void paint_binmap();
  void write_binmap(const string &filename);

private:
  vector<double> m_contours;
  vector<point_vector*> m_contour_points;

  CFITSImage m_smoothed_image, m_out_image;
  CFITSPosn m_posn;

  int m_xw, m_yw;
};

//////////////////////////////////////////////////////////////////////

contour_prog::contour_prog()
{
}

contour_prog::~contour_prog()
{
  for(int i=0; i<int(m_contour_points.size()); ++i)
    delete m_contour_points[i];

  m_contour_points.clear();
}

void contour_prog::run()
{
  read_contours("contours.txt");
  read_smoothed_image("smoothed.fits");

  make_contours();
  paint_binmap();

  write_binmap("contour_binmap.fits");
}

//////////////////////////////////////////////////////////////////////

void contour_prog::read_contours(const string &filename)
{
  ifstream clist(filename.c_str());
  if( ! clist ) {
    cerr << "Unable to open file " << filename << endl;
    exit(1);
  }

  m_contours.push_back(-1e30);

  while( ! clist.eof() ) {
    string line;
    getline(clist, line);

    // repeat if end of file, empty line or line starting with #
    if( ! clist) continue;
    if( line.empty() ) continue;
    if( line[0] == '#' ) continue;

    istringstream is( line );

    double cont;
    is >> cont;
    if( ! is ) {
      cerr << "Invalid number in contour level file " << filename << endl;
      exit(1);
    }

    m_contours.push_back(cont);
  }

  m_contours.push_back(1e30);

  cout << "Read list of contours in " << filename << endl;

  for(int i=0; i<int(m_contours.size()); ++i)
    m_contour_points.push_back( new point_vector );
}

void contour_prog::read_smoothed_image(const string &filename)
{
  CFITSFile imagefile(filename.c_str(), CFITSFile::existingro);

  m_smoothed_image = imagefile.GetImage();
  m_posn = imagefile.GetPosn();

  m_xw = m_smoothed_image.GetXW();
  m_yw = m_smoothed_image.GetYW();

  m_out_image.Resize(m_xw, m_yw);
}

void contour_prog::make_contours()
{
  // find out which contour each pixel represents
  for(int y=0; y<m_yw; ++y)
    for(int x=0; x<m_xw; ++x) {
      const double pix = m_smoothed_image.GetPixel(x, y);

      int cont = m_contours.size()-1;
      while( m_contours[cont]>pix && cont >= 0 )
	--cont;

      ++cont;

      m_contour_points[cont] -> push_back( point(x, y) );
    }

  for(int i = m_contour_points.size()-1;  i >= 0; --i)
    check_contour_contig(i);
}

void contour_prog::check_contour_contig(int cont)
{
  vector<point_vector*> split_bins; 

  CFITSImage pixels(m_xw, m_yw);
  point_vector *bin = m_contour_points[cont];

  const int nopixels = bin->size();
  if(nopixels == 0) return;

  vector<bool> pixel_used;
  for(int i=nopixels-1; i>=0; i--)
    pixel_used.push_back(false);
  
  for(;;) {

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
    point_vector *chunk = new point_vector;
    const point &first = (*bin)[next_pixel];
    chunk->push_back(first);
    bin_image[ first.x() + xw*first.y() ] = true;
    
    for(;;) {  // loop, adding contig pixels, until there are no more
      bool added_any = false;
      
      // loop over pixels in input bin
      for(int pix=0; pix<nopixels; pix++) {
	const point &curr = (*bin)[pix];
	
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
	    chunk->push_back(curr);
	    bin_image[curr.x()+curr.y()*xw] = true;
	  }
	  
	} // pixel_used
	
      } // loop over pixels
      
      if(!added_any) break;
    } // repeating loop over pixels
    
    split_bins.push_back(chunk);
    
    delete[] bin_image;

  }  // try another chunk

  // get rid of original bin
  delete bin;

  // add new split bins on
  m_contour_points[cont] = split_bins[0];
  for(int i=1; i<int(split_bins.size()); i++)
    m_contour_points.push_back( split_bins[i] );

  split_bins.clear();
}

void contour_prog::paint_binmap()
{
  int val = 0;

  for(int i=m_contour_points.size()-1; i >= 0; --i) {
    const point_vector *bin = m_contour_points[i];

    if(bin->empty())
      continue;

    int count = 0;
    for(int j=bin->size()-1; j >= 0; --j) {
      const point &p = (*bin)[j];
      
      m_out_image.SetPixel( p.x(), p.y(), val );
      ++count;
    }
    ++val;
  }
}


void contour_prog::write_binmap(const string &filename)
{
  CFITSFile fileout(filename.c_str(), CFITSFile::create);

  fileout.SetImage(m_out_image);
  fileout.SetPosn(m_posn);

  fileout.WriteImage();
}


int main()
{
  contour_prog program;
  program.run();

  return 0;
}
