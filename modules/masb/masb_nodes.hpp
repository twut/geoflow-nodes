#include <geoflow/core/geoflow.hpp>
#include"kdTree.h"
#include <compute_ma_processing.h>
#include <compute_normals_processing.h>


# define M_PI           3.14159265358979323846 


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
          add_input("ma_radii", typeid(vec1f));
          add_output("interior_mat", typeid(PointCollection));
          add_output("interior_radii", typeid(vec1f));
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
  class VisibiltyQurey :public Node 
  {
  public: 
      using Node::Node;
      void init() {
          add_input("KDTree", typeid(KdTree));
          add_input("interior_MAT", typeid(PointCollection));
          add_input("ViewPoint", typeid(Vector3D));
          add_input("interior_radii", typeid(vec1f));
          add_output("Visible_MAT", typeid(PointCollection));
          add_output("Radii_of_MAT", typeid(vec1f));
      }
      void process();
      float DistanceOfPointToLine(Vector3D a, Vector3D b, Vector3D s)
      {
          float ab = sqrt(pow((a.x - b.x), 2.0) + pow((a.y - b.y), 2.0) + pow((a.z - b.z), 2.0));
          float as = sqrt(pow((a.x - s.x), 2.0) + pow((a.y - s.y), 2.0) + pow((a.z - s.z), 2.0));
          float bs = sqrt(pow((s.x - b.x), 2.0) + pow((s.y - b.y), 2.0) + pow((s.z - b.z), 2.0));
          float cos_A = (pow(as, 2.0) + pow(ab, 2.0) - pow(bs, 2.0)) / (2 * ab*as);
          float sin_A = sqrt(1 - pow(cos_A, 2.0));
          return as * sin_A;
      }

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
          add_output("vector", typeid(Vector3D));
          add_output("point", typeid(PointCollection));

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

          PointCollection point;          
          point.push_back({ x, y, z });
                    
          output("vector").set(viewpoint);
          output("point").set(point);
          
      }
  };
  class Triangulation :public Node {
  public:
      using Node::Node;
      TriangleCollection *mp_tc = nullptr;
      vec3f *mp_normals = nullptr;
      

      void init() {
          add_input("MAT_points", typeid(PointCollection));
          add_input("radii",typeid(vec1f));
          add_output("triangle_collection", typeid(TriangleCollection));
          add_output("normals", typeid(vec3f));
      }
      void gui() {}
      void process()
      {
          auto mat_points = input("MAT_points").get<PointCollection>();
          auto radii = input("radii").get<vec1f>();
          TriangleCollection all_tc;
          vec3f all_normals;

          for (int i = 0; i < mat_points.size(); i++) {
              std::cout << "index i: " <<i<< std::endl;
              auto  tc= Triangulation::StandardSphere(mat_points[i], radii[i]);
              for (auto a : tc) {
                  all_tc.push_back(a);
              }                  
          }
          all_normals = ComputeNormals(all_tc);

          std::cout << "Triangulation Done" << std::endl;
          output("triangle_collection").set(all_tc);
          output("normals").set(all_normals);
      };
      
      geoflow::TriangleCollection StandardSphere(geoflow::arr3f center, float radius) {
          geoflow::Triangle t1,t2;
          geoflow::TriangleCollection tc;
          int Density = 10;
          std::array<float, 3> points[11][21];
          std::cout << "start triangulation" << std::endl;
          std::cout << "debug test" << std::endl;
          for (int t = 0; t <= Density; t++)
          {
              double vt = t * (M_PI / Density);
              for (int n = 0; n <= Density*2; n++) {
                  double vp = (n - Density) * (M_PI / Density);
                  points[t][n][0] = std::sin(vt)*std::cos(vp)*radius + center[0];
                  points[t][n][1] = std::sin(vt)*std::sin(vp)*radius + center[1];
                  points[t][n][2] = std::cos(vt)*radius + center[2];
                  //std::cout << "points:" << points[t][n][0] << "," << points[t][n][1] << "," << points[t][n][2] << std::endl;
              }
          }
          for (int nt = 1; nt <= Density; nt++) {
              for (int np = 0; np <= (2*Density-1); np++) {
                  tc.push_back({ points[nt][np],points[nt][np + 1],points[nt - 1][np] });
                  //(*mp_normals).push_back({ points[nt][np][0],points[nt][np][1],points[nt][np][2] });

                  tc.push_back ({ points[nt][np + 1],points[nt - 1][np + 1],points[nt - 1][np] });
                  //(*mp_normals).push_back({ points[nt][np][0],points[nt][np][1],points[nt][np][2] });
                
              }
          }
          std::cout << "size of tc:" << tc.size() << std::endl;
          return tc;          
      }
      vec3f ComputeNormals(geoflow::TriangleCollection tc) {
          vec3f normals;
          for (auto& t : tc) {
              masb::Vector a = (t[0].data());
              masb::Vector b = (t[1].data());
              masb::Vector c = (t[2].data());
              //auto v1 = b - a;
              auto n = Vrui::Geometry::cross(b - a, c - b);

              normals.push_back({ n[0],n[1],n[2] });
              normals.push_back({ n[0],n[1],n[2] });
              normals.push_back({ n[0],n[1],n[2] });
          }
          return normals;
      
      }
  };
  class TriangleNode :public Node {
  public:

      using Node::Node;
      void init() {
          add_output("triangle_collection", typeid(TriangleCollection));
          add_output("normals", typeid(vec3f));
      }
      void gui() {
      }
      
      void process() {
          typedef std::array<float, 3> point;
          point p0 = { -1.0f, -1.0f, -1.0f };
          point p1 = { 1.0f, -1.0f, -1.0f };
          point p2 = { 1.0f, 1.0f, -1.0f };
          point p3 = { -1.0f, 1.0f, -1.0f };

          point p4 = { -1.0f, -1.0f, 1.0f };
          point p5 = { 1.0f, -1.0f, 1.0f };
          point p6 = { 1.0f, 1.0f, 1.0f };
          point p7 = { -1.0f, 1.0f, 1.0f };

          TriangleCollection tc;
          tc.push_back({ p2,p1,p0 });
          tc.push_back({ p0,p3,p2 });
          tc.push_back({ p4,p5,p6 });
          tc.push_back({ p6,p7,p4 });
          tc.push_back({ p0,p1,p5 });
          tc.push_back({ p5,p4,p0 });
          tc.push_back({ p1,p2,p6 });
          tc.push_back({ p6,p5,p1 });
          tc.push_back({ p2,p3,p7 });
          tc.push_back({ p7,p6,p2 });
          tc.push_back({ p3,p0,p4 });
          tc.push_back({ p4,p7,p3 });

          vec3f normals;
          //counter-clockwise winding order
          for (auto& t : tc) {
              masb::Vector a = (t[0].data());
              masb::Vector b = (t[1].data());
              masb::Vector c = (t[2].data());
              //auto v1 = b - a;
              auto n = Vrui::Geometry::cross(b - a, c - b);

              normals.push_back({ n[0],n[1],n[2] });
              normals.push_back({ n[0],n[1],n[2] });
              normals.push_back({ n[0],n[1],n[2] });
          }
          output("triangle_collection").set(tc);
          output("normals").set(normals);
      }

  };
  



}
