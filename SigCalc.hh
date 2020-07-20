#ifndef ADAPTIVEBLOCK_SIGCALC_HH
#define ADAPTIVEBLOCK_SIGCALC_HH

namespace AdaptiveBin {

  // value specifies the possible number of results

  class CSigCalc {
  public:
    CSigCalc();
    virtual ~CSigCalc();

    typedef list<CCoord> cCoordList;
    typedef list<CCoordV> cCoordVList;

    virtual int GetNoValues() = 0;
    // get no possible values

    virtual int GetXW() = 0;
    virtual int GetYW() = 0;
    virtual double GetValue(const cCoordList &points,
			    int value = 0) = 0;
    virtual double GetError(const cCoordList &points,
			    int value = 0) = 0;

    virtual void Mask( CCoord p ) = 0;
    virtual bool IsMasked( CCoord p ) = 0;
  };


  class CCountSig : public CSigCalc {
  public:
    CCountSig(const CFITSImage &image);

    virtual int GetNoValues();
    // 0 - count

    virtual int GetXW();
    virtual int GetYW();
    virtual double GetValue(const cCoordList &points,
			    int value = 0);
    virtual double GetError(const cCoordList &points,
			    int value = 0);
    virtual void Mask( CCoord p );
    virtual bool IsMasked( CCoord p );

  private:
    const CFITSImage m_image;
    CFITSImage m_maskImage;
  };

  class CRatioSig : public CSigCalc {
  public:
    CRatioSig(const CFITSImage &imagea,
	      const CFITSImage &imageb);

    virtual int GetNoValues();
    // 0 - ratio, 1 = counta, 2 = countb

    virtual int GetXW();
    virtual int GetYW();
    virtual double GetValue(const cCoordList &points,
			    int value = 0);
    virtual double GetError(const cCoordList &points,
			    int value = 0);
    virtual void Mask( CCoord p );
    virtual bool IsMasked( CCoord p );

  protected:
    const CFITSImage m_imagea, m_imageb;

  private:
    CFITSImage m_maskImage;
  };

  class CRatioSig3 : public CRatioSig {
  public:
    CRatioSig3(const CFITSImage &imagea,
	       const CFITSImage &imageb,
	       const CFITSImage &imagec);

    virtual int GetNoValues();
    // 0 - ratio a/b (using error of a/b/c)
    // 1 - ratio a/b, 2 - ratio b/c, 3 - ratio c/a
    // 4 - count a, 5 - countb, 6 - countc

    virtual double GetValue(const cCoordList &points,
			    int value = 0);
    virtual double GetError(const cCoordList &points,
			    int value = 0);
  private:
    const CFITSImage m_imagec;
  };

  class CRatioSig4 : public CRatioSig {
  public:
    CRatioSig4(const CFITSImage &imagea,
	       const CFITSImage &imageb,
	       const CFITSImage &imagec,
	       const CFITSImage &imaged);

    virtual int GetNoValues();
    // 0 - ratio a/b (using average error of a/b/c/d)
    // 1 - ratio a/b, 2 - ratio b/c, 3 - ratio c/d, 4 - ratio d/a
    // 5 - count a, 6 - countb, 7 - countc, 8 - countd

    virtual double GetValue(const cCoordList &points,
			    int value = 0);
    virtual double GetError(const cCoordList &points,
			    int value = 0);
  private:
    const CFITSImage m_imagec, m_imaged;
  };

}

#endif
