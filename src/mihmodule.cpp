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

  //MihasherPython(int B, int m ):mihasher(10, 20) {}
  MihasherPython(): mihasher() {}
  MihasherPython(int B, int m):mihasher(B, m ) {}
 
  void setK_py(int K) {
      printf("set k to be %d \n", K);
      this->setK(K);
      return; 
  }

};


BOOST_PYTHON_MODULE(libmihmodule)
{
  class_<MihasherPython> ("MihasherPython")
      .def(init<int, int>())
      .def("setK", &MihasherPython::setK_py);
}
