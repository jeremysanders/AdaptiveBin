// program to take two binmaps and merge
// creating new bins from the intersections

#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <FITSFile.h>

class merge_program
{
public:
  merge_program(const string &binmap1, const string &binmap2);
  ~merge_program();

  void merge();
  void write_out(const string &filename);

private:
  CFITSImage m_binmap1, m_binmap2;
  CFITSPosn m_posn;
  int m_xw, m_yw;

  CFITSImage m_binmap_out;
};

///////////////////////////////////////////////////////////////////

merge_program::merge_program(const string &binmap1, const string &binmap2)
{
  { // read in binmap files
    CFITSFile f(binmap1.c_str(), CFITSFile::existingro);
    m_binmap1 = f.GetImage();
    m_posn = f.GetPosn();
  }{
    CFITSFile f(binmap2.c_str(), CFITSFile::existingro);
    m_binmap2 = f.GetImage();
  }

  m_xw = m_binmap1.GetXW();
  m_yw = m_binmap1.GetYW();

  m_binmap_out.Resize(m_xw, m_yw);
}

merge_program::~merge_program()
{
}

void merge_program::merge()
{
  // algorithm is compute hash value unique to intersection of
  // binmaps, and compute new bin numbers for each intersection

  // a map is used to map hash values to bin numbers

  ofstream bin_list_file1("binlist_1.txt");
  ofstream bin_list_file2("binlist_2.txt");

  map<int,int> bin_lookup;

  int next_bin = 0;

  for(int y=0; y<m_yw; ++y)
    for(int x=0; x<m_xw; ++x) {

      int outval = -1;

      const int bin1 = int(m_binmap1.GetPixel(x, y));
      const int bin2 = int(m_binmap2.GetPixel(x, y));

      assert(bin1 < 65536);

      const int hashkey = bin1 + bin2*65536;

      if( bin_lookup.find(hashkey) == bin_lookup.end() ) {
	outval = next_bin;
	bin_lookup[hashkey] = next_bin;
	++next_bin;

	bin_list_file1 << bin1 << '\t'
		       << bin_lookup[hashkey]
		       << endl;
	bin_list_file2 << bin2 << '\t'
		       << bin_lookup[hashkey]
		       << endl;
      } else {
	outval = bin_lookup[hashkey];
      }

      m_binmap_out.SetPixel(x, y, double(outval));
    }

}

void merge_program::write_out(const string &filename)
{
  CFITSFile f(filename.c_str(), CFITSFile::create);

  f.SetImage(m_binmap_out);
  f.SetPosn(m_posn);
  f.WriteImage();
}

int main(int argc, char *argv[])
{
  if(argc != 4) {
    cerr << "Usage: MergeBinMap binmap1.fits binmap2.fits out.fits\n";
    return -1;
  }

  merge_program prog(argv[1], argv[2]);
  prog.merge();
  prog.write_out(argv[3]);

  return 0;
}
