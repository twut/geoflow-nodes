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
      
      using Node::Node;
      void init() {
          add_input("Vector1", typeid(Vector3D));          
          add_output("MAT_points", typeid(PointCollection));
          

      }
      void gui() {
          
    
      }
      static std::vector<Vector3D> SpherePoints(Vector3D v1, float radius) {
          std::vector<Vector3D> points;
          Vector3D point;
          int nLongitude = 100;
          int nLatitude = 2 * nLongitude;
          int p, s, i, j;
          float x, y, z, out;
          int nPitch = nLongitude + 1;
          float DEGS_TO_RAD = 3.14159f / 180.0f;

          float pitchInc = (180. / (float)nPitch) * DEGS_TO_RAD;
          float rotInc = (360. / (float)nLatitude) * DEGS_TO_RAD;
          for (p = 1; p < nPitch; p++)     // Generate all "intermediate vertices":
          {
              out = radius * sin((float)p * pitchInc);
              if (out < 0) out = -out;    // abs() command won't work with all compilers
              y = radius * cos(p * pitchInc);
              //printf("OUT = %g\n", out);    // bottom vertex
              //printf("nPitch = %d\n", nPitch);    // bottom vertex
              for (s = 0; s < nLatitude; s++)
              {
                  x = out * cos(s * rotInc);
                  z = out * sin(s * rotInc);
                  point.x = x + v1.x;
                  point.y = y + v1.y;
                  point.z = z + v1.z;
                  //outfile << x + pt.x << "," << y + pt.y << "," << z + pt.z << std::endl;
                  //numVertices++;
                  points.push_back(point);
              }
          }
          return points;
      }
      
      void process() {
          Vector3D v1 = input("Vector1").get<Vector3D>();

          PointCollection outpt;
          auto v2_list = testNode::SpherePoints(v1, 500);
          for (Vector3D a : v2_list) {
              outpt.push_back({a.x,a.y,a.z});

          }
          output("MAT_points").set(outpt);



      };
  };
  class MultiKDtree :public Node {
  public:
      using Node::Node;
     

      KdTree* NewBuildKdTree(Vector3D* Points, long number)
      {
          std::cout << "Start" << std::endl;
  
          float* m_Radius = new float[number];
          
          KdTree* mp_kdTree = new KdTree(Points, number,16);
          delete[] Points;
           return mp_kdTree;
      }
      void init() {
          add_input("points", typeid(PointCollection));
          add_output("KDTree1", typeid(KdTree));
          add_output("KDTree2", typeid(KdTree));
      }
      void gui() {

      }
      void process() ;


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
           
          float* m_Radius = new float[number];
          
          mp_kdTree = new KdTree(Points, number, nMaxBucketSize);
          delete[] Points;
          std::cout << "KDTree Built" << std::endl;
          return mp_kdTree;
      }
      std::vector<int> GetLevelPoints(std::vector<int> vec_long, std::vector<int> vec_short) {
          std::vector<int> vec_result(20);
          std::vector<int>::iterator it;
          it = std::set_difference(vec_long.begin(), vec_long.end(), vec_short.begin(), vec_short.end(), vec_result.begin());
          vec_result.resize(it - vec_result.begin());
      
          return vec_result;
      }

      std::vector<KdTree::sphere> GetLevelPoints(Vector3D max, Vector3D min, std::vector<KdTree::sphere> &points) {
          std::vector<KdTree::sphere> LevelPoints;
          std::vector<KdTree::sphere>::iterator it;
          for (auto a:points){
              if(a.pos.x<=max.x&&a.pos.x>=min.x)
                  if(a.pos.y<=max.y&&a.pos.y>=min.y)
                      if (a.pos.z <= max.z&&a.pos.z >= min.z) {

                          LevelPoints.push_back(a);
                          for (it = points.begin(); it != points.end();) {
                              if ((*it).pos == a.pos)
                                  it = points.erase(it);
                              else
                              {
                                  it++;
                              }
                          }
                      }
              
          }



          return LevelPoints;
      
      }

      void init() {
          add_input("points", typeid(PointCollection));
          add_input("radii", typeid(vec1f));
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
      static float DistanceOfPointToLine(Vector3D a, Vector3D b, Vector3D s)
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
              
              auto  tc= Triangulation::StandardSphere(mat_points[i], radii[i]);
              for (auto a : tc) {
                  all_tc.push_back(a);
              }                  
          }
          all_normals = ComputeNormals(all_tc);

          std::cout << "Triangulation Done:" <<mat_points.size()<< std::endl;
          output("triangle_collection").set(all_tc);
          output("normals").set(all_normals);
      };
      
      geoflow::TriangleCollection StandardSphere(geoflow::arr3f center, float radius) {
          geoflow::Triangle t1,t2;
          geoflow::TriangleCollection tc;
          int Density = 10;
          std::array<float, 3> points[11][21];
          //std::cout << "start triangulation" << std::endl;
          //std::cout << "debug test" << std::endl;
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
          //std::cout << "size of tc:" << tc.size() << std::endl;
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
  
  class OneQuery:public Node  {
  public:
      using Node::Node;
      void init() {
          add_input("KDTree", typeid(KdTree));
          add_input("MATpoints", typeid(PointCollection));
          add_input("radii", typeid(vec1f));
          add_input("Vector1", typeid(Vector3D));
          //add_input("Vector2", typeid(Vector3D));
          add_output("MAT_points", typeid(PointCollection));
          add_output("radii", typeid(vec1f));
          
      }
      void process();

      static float PointToPointDis(Vector3D p1, Vector3D p2) {         
          float dis = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));         
          return dis;
      }

      int inline GetIntersection(float fDst1, float fDst2, Vector3D P1, Vector3D P2, Vector3D &Hit) {
          if ((fDst1 * fDst2) >= 0.0f) return 0;
          if (fDst1 == fDst2) return 0;
          Hit = P1 + (P2 - P1) * (-fDst1 / (fDst2 - fDst1));
          return 1;
      }

      int inline InBox(Vector3D Hit, Vector3D B1, Vector3D B2, const int Axis) {
          if (Axis == 1 && Hit[2] > B1[2] && Hit[2] < B2[2]&& Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          if (Axis == 2 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[0] > B1[0] && Hit[0] < B2[0]) return 1;
          if (Axis == 3 && Hit[0] > B1[0] && Hit[0] < B2[0] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          return 0;
      }
      int CheckLineBox(Vector3D B1, Vector3D B2, Vector3D L1, Vector3D L2, Vector3D &Hit)
      {
          if (L2[0] < B1[0] && L1[0] < B1[0]) return false;
          if (L2[0] > B2[0] && L1[0] > B2[0]) return false;
          if (L2[1] < B1[1] && L1[1] < B1[1]) return false;
          if (L2[1] > B2[1] && L1[1] > B2[1]) return false;
          if (L2[2] < B1[2] && L1[2] < B1[2]) return false;
          if (L2[2] > B2[2] && L1[2] > B2[2]) return false;
          if (L1[0] > B1[0] && L1[0] < B2[0] &&
              L1[1] > B1[1] && L1[1] < B2[1] &&
              L1[2] > B1[2] && L1[2] < B2[2])
          {
              Hit = L1;
              return true;
          }
          if ((GetIntersection(L1[0] - B1[0], L2[0] - B1[0], L1, L2, Hit) && InBox(Hit, B1, B2, 1))
              || (GetIntersection(L1[1] - B1[1], L2[1] - B1[1], L1, L2, Hit) && InBox(Hit, B1, B2, 2))
              || (GetIntersection(L1[2] - B1[2], L2[2] - B1[2], L1, L2, Hit) && InBox(Hit, B1, B2, 3))
              || (GetIntersection(L1[0] - B2[0], L2[0] - B2[0], L1, L2, Hit) && InBox(Hit, B1, B2, 1))
              || (GetIntersection(L1[1] - B2[1], L2[1] - B2[1], L1, L2, Hit) && InBox(Hit, B1, B2, 2))
              || (GetIntersection(L1[2] - B2[2], L2[2] - B2[2], L1, L2, Hit) && InBox(Hit, B1, B2, 3)))
              return true;

          return false;
      }

       static std::vector<Vector3D> SpherePoints(Vector3D v1, float radius) {
          std::vector<Vector3D> points;
          Vector3D point;
          int nLongitude = 100;
          int nLatitude = 2 * nLongitude;
          int p, s, i, j;
          float x, y, z, out;
          int nPitch = nLongitude + 1;
          float DEGS_TO_RAD = 3.14159f / 180.0f;

          float pitchInc = (180. / (float)nPitch) * DEGS_TO_RAD;
          float rotInc = (360. / (float)nLatitude) * DEGS_TO_RAD;
          for (p = 1; p < nPitch; p++)     // Generate all "intermediate vertices":
          {
              out = radius * sin((float)p * pitchInc);
              if (out < 0) out = -out;    // abs() command won't work with all compilers
              y = radius * cos(p * pitchInc);
              //printf("OUT = %g\n", out);    // bottom vertex
              //printf("nPitch = %d\n", nPitch);    // bottom vertex
              for (s = 0; s < nLatitude; s++)
              {
                  x = out * cos(s * rotInc);
                  z = out * sin(s * rotInc);
                  point.x = x + v1.x;
                  point.y = y + v1.y;
                  point.z = z + v1.z;
                  //outfile << x + pt.x << "," << y + pt.y << "," << z + pt.z << std::endl;
                  //numVertices++;
                  points.push_back(point);
              }
          }
          return points;
      }

  };


}
