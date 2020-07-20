// a program to make a "mask image"
// which specifies which areas of an image /not/
// to look at when adaptive-binning

// --input_mask  - specify a mask to add things to
// --image       - specify an image to take dimensions from
// --output      - specify outimage

// supported shapes
// [+/~]rectangle(x1,y1,x2,y2) [+/~]circle(xc,yc,rad)
// coordinates are pixels

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <FITSFile.h>

#include <parammm/parammm.hh>

#include "version.hh"

using std::string;
using std::vector;
using std::sqrt;
using std::ostringstream;
using std::istringstream;
using std::cout;
using std::cerr;
using std::endl;

typedef vector<string> str_vec;
typedef vector<double> double_vec;

class mask
{
public:
  mask(double val);   // calls interpret_args
  virtual ~mask();

  void set_args(const str_vec &args);
  virtual string get_description() = 0;
  virtual void draw(CFITSImage *image) = 0;
  virtual void interpret_args(const double_vec &args);

protected:
  const double m_val;
};

mask::mask(double val)
  : m_val(val)
{
}

void mask::set_args(const str_vec &args)
{
  // convert args to doubles, and pass to interpret_args
  double_vec arg_list;
  for(int i=0; i<int(args.size()); ++i) {
    istringstream in(args[i]);

    double d;
    in >> d;
    if(!in) {
      cerr << args[i] << " is not a number. Aborting\n";
      exit(1);
    }

    arg_list.push_back( d );
  }

  interpret_args( arg_list );
}

mask::~mask()
{
}

void mask::interpret_args(const double_vec &args)
{
}

////////////////////////////////////////////////////////

class rectangle_mask : public mask
{
public:
  rectangle_mask(double val);
  void draw(CFITSImage *image);
  string get_description();
  void interpret_args(const double_vec &args);

private:
  int m_x1, m_y1, m_x2, m_y2;
};

rectangle_mask::rectangle_mask(double val)
  : mask(val)
{
}

void rectangle_mask::draw(CFITSImage *image)
{
  for(int y=m_y1; y<=m_y2; ++y)
    for(int x=m_x1; x<=m_x2; ++x)
      image->SetPixelBounded(x, y, m_val);
}

string rectangle_mask::get_description()
{
  ostringstream o;
  if( fabs(m_val) < 1e-8 ) o << '~'; else o << '+';

  o << "rectangle(" << m_x1 << ',' << m_y1 << ','
    << m_x2 << ',' << m_y2 << ')' << '\0';
  return o.str();
}

void rectangle_mask::interpret_args(const double_vec &args)
{
  if( args.size() != 4 ) {
    cerr << "4 arguments required for rectangle shape. Aborting\n";
    exit(1);
  }

  m_x1 = int(args[0]); m_y1 = int(args[1]);
  m_x2 = int(args[2]); m_y2 = int(args[3]);
}


class circle_mask : public mask
{
public:
  circle_mask(double val);
  void draw(CFITSImage *image);
  void interpret_args(const double_vec &args);
  string get_description();

private:
  double m_xc, m_yc;
  double m_radius;
};

circle_mask::circle_mask(double val)
  : mask(val)
{
}

string circle_mask::get_description()
{
  ostringstream o;
  if( fabs(m_val) < 1e-8 ) o << '-'; else o << '+';

  o << "circle(" << m_xc << ',' << m_yc << ','
    << m_radius << ')' << '\0';
  return o.str();
}

void circle_mask::draw(CFITSImage *image)
{
  const int xw = image->GetXW(), yw = image->GetYW();

  // crap algorithm, but it works
  for(int y=0; y<yw; ++y)
    for(int x=0; x<xw; ++x) {
      const double dx = x - m_xc, dy = y - m_yc;

      // calc distance, and see how far away it is 
      const double dist = sqrt( dx*dx + dy*dy );
      if(dist < m_radius)
	image->SetPixel(x, y, m_val);
    }

}

void circle_mask::interpret_args(const double_vec &args)
{
  if( args.size() != 3 ) {
    cerr << "3 arguments required for circle shape. Aborting\n";
    exit(1);
  }

  m_xc = args[0]; m_yc = args[1];
  m_radius = args[2];
}

////////////////////////////////////////////////////////////////

typedef vector<mask*> mask_vec;

class mask_prog
{
public:
  mask_prog(int argc, char **argv);
  ~mask_prog();

  void run();

private:
  void interpret_mask(const string &mask_str);
  void get_input_image_size();
  void draw_masks();
  void write_total_mask();

private:
  string m_output_file, m_input_file;
  mask_vec m_masks;

  int m_xw, m_yw;   // size of input image
  CFITSImage m_mask_image;
  CFITSPosn m_posn;
};

mask_prog::mask_prog(int argc, char **argv)
  : m_output_file("adbin_mask.fits")
{
  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch("out", 'o',
				      parammm::pstring_opt(&m_output_file),
				      "set out file (def adbin_mask.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("in", 'i',
				      parammm::pstring_opt(&m_input_file),
				      "set input file (required)", "FILE"));

  params.set_autohelp("Usage: MakeMask [OPTIONS] \"+/~mask(param,...)\" ...\n"
		      "Creates a mask image.\n"
		      "Options for mask are:\n"
		      " rectangle(x1,y1,x2,y2) and circle(xc, yc, diameter)\n"
		      "Written by Jeremy Sanders, 2001.",
		      "Report bugs to <jss@ast.cam.ac.uk>");
  params.enable_autohelp();
  params.enable_autoversion(c_adbin_version,
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();

  params.interpret_and_catch();

  if( m_input_file.empty() || params.args().size() < 1 )
    params.show_autohelp();

  for(int arg=0; arg<int(params.args().size()); ++arg)
    interpret_mask( params.args()[arg] );
}

mask_prog::~mask_prog()
{
  for(int i=0; i<int(m_masks.size()); ++i)
    delete m_masks[i];
  m_masks.clear();
}

void mask_prog::get_input_image_size()
{
  CFITSFile f(m_input_file.c_str(), CFITSFile::existingro);
  m_xw = f.GetImage().GetXW();
  m_yw = f.GetImage().GetYW();
  m_posn = f.GetPosn();
}

void mask_prog::draw_masks()
{
  for(int i=0; i<int(m_masks.size()); ++i)
    m_masks[i]->draw(&m_mask_image);
}

void mask_prog::write_total_mask()
{
  CFITSFile f(m_output_file.c_str(), CFITSFile::create);
  f.SetImage(m_mask_image);
  f.SetPosn(m_posn);
  f.WriteImage();
  f.WriteHistory("Mask file generated by MakeMask for AdaptiveBin");
  const string s = "Mask generated for " + m_input_file;
  f.WriteHistory(s.c_str());
  f.WriteHistory("Parameters for mask listed below");
  for(int i=0; i<int(m_masks.size()); ++i) {
    const string s = "Mask: " + m_masks[i] -> get_description();
    f.WriteHistory(s.c_str());
  }
}

void mask_prog::run()
{
  get_input_image_size();
  m_mask_image.Resize(m_xw, m_yw);
  m_mask_image.SetAll(0.);   // default, but do it again
  draw_masks();
  write_total_mask();
}

void mask_prog::interpret_mask(const string &mask_str)
{
  const int no_names = 2;
  const char *names[no_names] = {"rectangle", "circle"};

  string mask_s = mask_str;       // copy mask_str
  double val = 1.;              // value to write in output

  if( mask_s.empty() ) {
    cerr << "Blank mask specified. Aborting\n";
    exit(1);
  }

  if( mask_s[0] == '~' ) {
    val = 0;
    mask_s.erase(mask_s.begin());
  }
  if( mask_s[0] == '+' ) {
    val = 1;
    mask_s.erase(mask_s.begin());
  }

  // find mask type
  const unsigned int bracket = mask_s.find_first_of('(');
  if(bracket == string::npos || mask_s.find_first_of(')') == string::npos ) {
    cerr << "Invalid syntax for mask. Usage name(1,2,3,...).\nAborting\n";
    exit(1);
  }

  const string name = mask_s.substr(0, bracket);
  int mask_id = 0;
  while( mask_id < no_names && name != names[mask_id] )
    ++mask_id;

  if( mask_id == no_names ) {
    cerr << "Mask type not found. Valid types are:\n";
    for(int i=0; i<no_names; ++i)
      cerr << names[i] << ' ';
    cerr << "\nAborting\n";
    exit(1);
  }

  // Now, try and build list of parameters to mask
  const string params = mask_s.substr(bracket+1, string::npos);
  str_vec mask_parameters;

  const int len = params.size()-1;
  string build_str;
  for(int i=0; i<len; ++i) {
    const char c = params[i];
    if( c != ' ' && c != ',' && c != '\t' )
      build_str += c;

    if( (c==',' || i == len-1) && !(build_str.empty()) ) {
      mask_parameters.push_back( build_str );
      build_str.erase();
    }
  }

  mask *mask_ptr = 0;
  switch( mask_id ) {
  case 0: // rectangle
    mask_ptr = new rectangle_mask( val );
    break;
  case 1: // circle
    mask_ptr = new circle_mask( val );
    break;
  default:
    assert(false);
  }

  mask_ptr -> set_args(mask_parameters);

  cout << mask_ptr -> get_description() << endl;

  m_masks.push_back( mask_ptr );
}

////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  mask_prog prog(argc, argv);
  prog.run();

  return 0;
}
