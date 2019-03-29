#include "masb_nodes.hpp"

namespace geoflow::nodes::mat {

  NodeRegister create_register() {
    NodeRegister R("MAT");
    R.register_node<ComputeMedialAxisNode>("ComputeMedialAxisNode");
    R.register_node<ComputeNormalsNode>("ComputeNormalsNode");
    R.register_node<testNode>("TestNode");
    R.register_node<BuildKDtree>("BuildKDtree");
    R.register_node<KDTreeNearestQurey>("KDTreeNearestQuery");
    R.register_node<KDTreeLineQurey>("KDTreeLineQuery");
    R.register_node<NumberNode>("Number");
    R.register_node<ViewPoint>("ViewPoint");

    return R;
  }

}