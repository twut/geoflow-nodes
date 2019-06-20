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
    R->register_node<MATsimplification>("MATsimplification");
    R->register_node<OffsetNode>("Offset");
    R->register_node<VisibiltyQurey>("VisbiltyQurey");
    R->register_node<Triangulation>("Triangulation");    
    R->register_node<MultiKDtree>("MultiKDtree");
    R->register_node<OneQuery>("Onequery");
    R->register_node<MutiThreadsOneQuery>("MultiThreadsQuery");
    R->register_node<AMPGPUQueryTest>("GPUQuery");
    R->register_node<FromMATtoPointCloud>("FromMATtoPoints");
    R->register_node<VisiblePC>("GetVisiblePC");
    R->register_node<ParallelVector>("ParallelVector");    
    R->register_node<GetRaysResult>("GetRaysResult");
    R->register_node<VisiblePart>("GetVisiblePart");
    

    return R;
  }

}