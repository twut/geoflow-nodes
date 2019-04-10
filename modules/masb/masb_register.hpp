#include "masb_nodes.hpp"

namespace geoflow::nodes::mat {

  NodeRegisterHandle create_register() {
    auto R = NodeRegister::create("MAT");
    R->register_node<ComputeMedialAxisNode>("ComputeMedialAxisNode");
    R->register_node<ComputeNormalsNode>("ComputeNormalsNode");
    R->register_node<testNode>("TestNode");
    R->register_node<BuildKDtree>("BuildKDtree");
    R->register_node<KDTreeNearestQurey>("KDTreeNearestQuery");
    R->register_node<KDTreeLineQurey>("KDTreeLineQuery");
    R->register_node<NumberNode>("Number");
    R->register_node<ViewPoint>("ViewPoint");
    R->register_node<MATfilter>("MATfilter");
    R->register_node<VisibiltyQurey>("VisbiltyQurey");
    R->register_node<Triangulation>("Triangulation");
    return R;
  }

}