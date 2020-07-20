//      Adaptive Binning Program
//      Binning module header - class to allow different methods of binning
//                              to be used
//      Described in Sanders and Fabian (submitted)
//      Routines for adaptively binning data
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

// binmodule.hh
// modules for binning files...

#include <vector>
#include <string>

#include <FITSImage.h>
#include <FITSPosn.h>

namespace AdaptiveBin
{
  class pixel
  {
  public:
    pixel(int x, int y) { m_x = x; m_y = y; }
    int x() const { return m_x; }
    int y() const { return m_y; }
  private:
    int m_x, m_y;
  };

  typedef std::vector<pixel> pixlist;
  typedef std::vector<std::string> arglist;

  class binmodule
  {
  public:
    virtual ~binmodule();

    virtual double fracerror(const pixlist &pl,
			     bool binerror) = 0;
    virtual double value(const pixlist &pl) = 0;

    virtual void getposn(CFITSPosn *out) = 0;

    virtual void selectvalue(const std::string &spec) = 0;
    virtual std::string get_value_descr() = 0;

    virtual int xw() = 0;
    virtual int yw() = 0;
  };

  class invalidargs_exception
  {};
  class invalidvalue_exception
  {};

  // normal count module
  class count_binmodule : public binmodule
  {
  public:
    count_binmodule(const arglist &al);
    count_binmodule(const std::string &fname,
		    double background);

    double fracerror(const pixlist &pl, bool binerror);
    double value(const pixlist &pl);

    void getposn(CFITSPosn *out);

    void selectvalue(const std::string &spec);
    std::string get_value_descr();

    int xw();
    int yw();

  private:
    void setf(const std::string &fname,
	      double background);

  private:
    CFITSImage m_image;
    double m_background;
    CFITSPosn m_posn;
  };

  // external module (data and error file)
  class external_binmodule : public binmodule
  {
  public:
    external_binmodule(const arglist &al);

    double fracerror(const pixlist &pl, bool binerror);
    double value(const pixlist &pl);

    void getposn(CFITSPosn *out);

    void selectvalue(const std::string &spec);
    std::string get_value_descr();

    int xw();
    int yw();

  private:
    CFITSImage m_image, m_error;
    CFITSPosn m_posn;
    bool m_absolute;  // absolute or relative error used
  };

  // bin module for the ratio of two files
  class ratio_binmodule : public binmodule
  {
  public:
    ratio_binmodule(const arglist &al);
    
    double fracerror(const pixlist &pl, bool binerror);
    double value(const pixlist &pl);

    void getposn(CFITSPosn *out);

    void selectvalue(const std::string &spec);
    std::string get_value_descr();

    int xw();
    int yw();

  private:
    enum valuet { vcount, vratio };
    valuet m_value;
    unsigned m_valparam[2];

    std::vector<count_binmodule> m_counts;
  };

}
