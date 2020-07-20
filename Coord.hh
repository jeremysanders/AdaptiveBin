// Coordinate type

#ifndef COORD_ADAPTIVE_BIN_HH
#define COORD_ADAPTIVE_BIN_HH

namespace AdaptiveBin {

  // simple x, y coordinate
  class CCoord
  {
  public:
    CCoord(int x, int y);
    CCoord();
    virtual ~CCoord();
    CCoord operator +(const CCoord &other) const;

    int m_x, m_y;
  };

  inline CCoord::CCoord(int x, int y) :
    m_x(x), m_y(y)
  {
  }

  inline CCoord::CCoord() :
    m_x(0), m_y(0)
  {
  }

  inline CCoord CCoord::operator+(const CCoord &other) const
  {
    return CCoord( m_x+other.m_x, m_y+other.m_y );
  }

  inline CCoord::~CCoord()
  {
  }

  // coordinate with value
  class CCoordV : public CCoord
  {
  public:
    CCoordV(int x, int y, double val = 0.0);
    CCoordV();

    double m_val;
  };

  inline CCoordV::CCoordV()
    : CCoord(), m_val(0.0)
  {
  }

  inline CCoordV::CCoordV(int x, int y, double val)
    : CCoord(x, y)
  {
    m_val = val;
  }

}

#endif
