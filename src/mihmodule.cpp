#include "Python.h"
#include <boost/python.hpp>
#include <exception>
#include <stdint.h>
#include "mihasher.h"

using namespace std;
using namespace boost;
using namespace boost::python;



class MihasherPython : public mihasher {
//class MihasherPython {
  public:

  MihasherPython(B, m):mihasher(B, m) {}
  void poppulate_py(UINT8* codes, UINT32 N, int dim1codes) {
     this->populate(codes, N, dim1codes); 
  }
 
  void setK_py(int K) {
      this->setK(K);
      return; 
  }

  void query(UINT32 *results, UINT32* numres, qstat *stats, UINT8 *q, UINT64 * chunks, UINT32 * res);
 
  
};


BOOST_PYTHON_MODULE(libmihmodule)
{
  class_<MihasherPython> ("MihasherPython")
      .def("setK", &MihasherPython::setK_py);
}
