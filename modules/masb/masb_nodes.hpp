#include <geoflow/core/geoflow.hpp>
#include"kdTree.h"
#include <compute_ma_processing.h>
#include <compute_normals_processing.h>

namespace geoflow::nodes::mat {


  class ComputeMedialAxisNode:public Node {
    public:
    masb::ma_parameters params;
    float interval = 2;
    double zero=0,pi=3.14;
    using Node::Node;
    void init() {
      add_input("points", typeid(PointCollection));
      add_input("normals", typeid(vec3f));
      add_output("ma_coords", typeid(PointCollection));
      add_output("ma_radii", typeid(vec1f));
      add_output("ma_qidx", typeid(vec1i));
      add_output("ma_is_interior", typeid(vec1i));
    }
    void gui(){
      ImGui::SliderFloat("initial_radius", &params.initial_radius, 0, 1000);
      ImGui::SliderScalar("denoise_preserve", ImGuiDataType_Double, &params.denoise_preserve, &zero, &pi);
      ImGui::SliderScalar("denoise_planar", ImGuiDataType_Double, &params.denoise_planar, &zero, &pi);
      ImGui::Checkbox("nan_for_initr", &params.nan_for_initr);
    }
    void process();
  };
  class MATfilter :public Node {
  public:
      using Node::Node;
      void init() {
          add_input("ma_coords", typeid(PointCollection));
          add_input("ma_is_interior", typeid(vec1i));
          add_output("interior_mat", typeid(PointCollection));
          add_output("exterior_mat", typeid(PointCollection));
      }
      void process();

  };

  class ComputeNormalsNode:public Node {
    public:
    masb::normals_parameters params;
    float interval = 2;
    using Node::Node;
    void init() {
      add_input("points", typeid(PointCollection));
      add_output("normals", typeid(vec3f));
    }
    void gui(){
      ImGui::SliderInt("K", &params.k, 1, 100);
    }
    void process();
  };
  class testNode :public Node {
  public:
      char filepath[256] = "";
      using Node::Node;
      void init() {
          add_input("in_radii", typeid(vec1f));
          add_output("out_radii", typeid(vec1f));

      }
      void gui() {
          ImGui::InputText("File path", filepath, IM_ARRAYSIZE(filepath));
    
      }
      
      void process();
  };
  class BuildKDtree :public Node {
  public:
      using Node::Node;
      KdTree *mp_kdTree;
      KdTree*             GetKdTree() { return mp_kdTree; }
      //number of points;
      Vector3D* mp_Points;
      long m_nPoints;

      KdTree* BuildKdTree( Vector3D* Points, long number, int nMaxBucketSize)
      {
          std::cout << "Start" << std::endl;
          //if (0 == number) return false;

          //Vector3D* pPoints = new Vector3D[m_nPoints];

          // save the radius not finished 
          
          float* m_Radius = new float[number];
 
          
          mp_kdTree = new KdTree(Points, number, nMaxBucketSize);
          delete[] Points;
          //SAFE_DELETE_ARRAY(pPoints);
          std::cout << "KDTree Built" << std::endl;


          return mp_kdTree;
      }

      void init() {
          add_input("points", typeid(PointCollection));
          add_output("KDTree", typeid(KdTree));
         

      }
      void gui() {

      }
      void process();
  };
  class KDTreeNearestQurey:public Node{
  public:
      using Node::Node;
      //KdTree* mp_kdTree;

      void init() {
          //add_input("viewpoint", TT_vec1f);
          add_input("KDTree", typeid(KdTree));
          add_input("MATpoints", typeid(PointCollection));
          add_input("Number", typeid(float));
          add_input("ViewPoint", typeid(Vector3D));
          add_output("Points", typeid(PointCollection));
          //add_output("result_points", TT_point_collection);

      }
      void process() ; 
  };
  class KDTreeLineQurey :public Node {
  public:
      using Node::Node;
      //KdTree* mp_kdTree;

      void init() {
          //add_input("viewpoint", TT_vec1f);
          add_input("KDTree", typeid(KdTree));
          add_input("MATpoints", typeid(PointCollection));
          add_input("Vector1", typeid(Vector3D));
          add_input("Vector2", typeid(Vector3D));
          add_input("Number", typeid(float));
          add_input("Distance", typeid(float));
          add_output("Points", typeid(PointCollection));
          //add_output("result_points", TT_point_collection);

      }
      void process();
  };
  class NumberNode :public Node {
  public:
      using Node::Node;

      void init() {
          add_output("result", typeid(float));
          add_param("number_value", (int)5);
      }

      void gui() {
          ImGui::InputInt("Number value", &param<int>("number_value"));
      }

      void process() {         
          output("result").set(float(param<int>("number_value")));        
      }
  };
  class ViewPoint :public Node {
  public:
      using Node::Node;

      void init() {
          add_output("result", typeid(Vector3D));

          add_param("x_value", (float)-99.0594);
          add_param("y_value", (float)-90.828);
          add_param("z_value", (float)6.59866);
      }

      void gui() {
          ImGui::InputFloat("x value", &param<float>("x_value"));
          ImGui::InputFloat("y value", &param<float>("y_value"));
          ImGui::InputFloat("z value", &param<float>("z_value"));
      }

      void process() {
          float x = param<float>("x_value");
          float y = param<float>("y_value");
          float z = param<float>("z_value");
          Vector3D viewpoint = {x,y,z};
          output("result").set(viewpoint);
          
      }
  };
  



}