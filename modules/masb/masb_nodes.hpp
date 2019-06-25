#include <geoflow/core/geoflow.hpp>
#include"kdTree.h"
#include "Vector3DNew.h"
#include <compute_ma_processing.h>
#include <compute_normals_processing.h>
#include <thread>
#include <mutex>
#include <amp.h>
#include<amp_math.h>

# define M_PI           3.14159265358979323846 


namespace geoflow::nodes::mat {
   

  static std::mutex mtx;
  static std::mutex mtx2;
  class FromMATtoPointCloud :public Node {
  public:
      using Node::Node;
      void init() 
      { 
          add_input("points", typeid(PointCollection));
          //add_input("Visible_MAT", typeid(PointCollection));

          //add_input("ma_qidx", typeid(vec1i));
          add_input("vis_idx", typeid(vec1i));
          

          add_output("Vis_points", typeid(PointCollection));
      }
      void gui() 
      {
      }
      void process() 
      {
          // ---------------input -------------------//
          auto input_points = input("points").get<PointCollection>();
          auto vis_idx = input("vis_idx").get<vec1i>();
          //---------------output ---------------//
          PointCollection output_points;
          
          for (int idx : vis_idx)
          {
              //std::cout << idx << std::endl;
              output_points.push_back(input_points[idx]);
          }

          output("Vis_points").set(output_points);         
      }

  };

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
      ImGui::SliderFloat("initial_radius", &params.initial_radius, 0, 1);
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
          //add_input("offset", typeid(float));
          //add_input("min_z", typeid(float));          
          
          add_output("exterior_mat", typeid(PointCollection));
          add_output("exterior_radii", typeid(vec1f));
          add_output("exterior_idx", typeid(vec1i));

          add_output("interior_mat", typeid(PointCollection));
          add_output("interior_radii", typeid(vec1f));
          add_output("interior_idx", typeid(vec1i));
      }
      void process();

  };
  class MATsimplification :public Node {
  public:
      using Node::Node;
      void init() {
          add_input("interior_mat", typeid(PointCollection));
          add_input("interior_radii", typeid(vec1f));
          add_input("interior_idx", typeid(vec1i));
          add_input("threshold", typeid(float));

          add_output("interior_mat", typeid(PointCollection));
          add_output("interior_radii", typeid(vec1f));
          add_output("interior_idx", typeid(std::vector<vec1i>));
      }
      void process();
      static bool ifSimplify(Vector3D v1, float radius,int index, float threshold, std::vector<KdTree::sphere> &MAT_to_KD)
      {
          //std::cout << "Test " << std::endl;
          KdTree::sphere tempSP;
          
          for (int i = 0; i < MAT_to_KD.size(); i++)
          {
              if (abs(MAT_to_KD[i].pos.x - v1.x) <= threshold && abs(MAT_to_KD[i].pos.y - v1.y) <= threshold && abs(MAT_to_KD[i].pos.z - v1.z) <= threshold && abs(MAT_to_KD[i].radius - radius) <= threshold)
              {
                  MAT_to_KD[i].index.push_back(index);
                  //break;
                  return 1;

              }
          }                                
          tempSP.pos = v1;
          tempSP.radius = radius;
          tempSP.index.push_back(index);

          MAT_to_KD.push_back(tempSP);
          //break;
          return 0;
              
          
         
      }
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
          int nLongitude = 1000;
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

      std::vector<KdTree::sphere> GetLevelPoints(Vector3D max, Vector3D min, std::vector<KdTree::sphere> *points) {
          std::vector<KdTree::sphere> LevelPoints;       
          std::vector<KdTree::sphere> diff;
          std::vector<KdTree::sphere>::iterator it;
          for (auto a:*points){
              if(a.pos.x<=max.x&&a.pos.x>=min.x)
                  if(a.pos.y<=max.y&&a.pos.y>=min.y)
                      if (a.pos.z <= max.z&&a.pos.z >= min.z) {
                          LevelPoints.push_back(a);
                          for (it = (*points).begin(); it != (*points).end();) {
                              if ((*it).pos == a.pos)
                                  it = (*points).erase(it);
                              else
                              {
                                  ++it;
                              }
                          }
                      }
              
          }          
          /*it = std::set_difference((*points).begin(), (*points).end(), LevelPoints.begin(), LevelPoints.end(), diff.begin());
          diff.resize(it - diff.begin());
          *points = diff;*/
          return LevelPoints;      
      }
      static bool BoxIntersectsSphere(Vector3D Bmin, Vector3D Bmax, Vector3D C, float r) {
          float r2 = r * r;
          float dmin = 0;
          for (int i = 0; i < 3; i++) {
              if (C[i] < Bmin[i]) dmin += ((C[i] - Bmin[i])*(C[i] - Bmin[i]));
              else if (C[i] > Bmax[i]) dmin += ((C[i] - Bmax[i])*(C[i] - Bmax[i]));
          }
          return dmin <= r2;
      }

      std::vector<KdTree::sphere> NewGetLevelPoints(Vector3D max, Vector3D min, std::vector<KdTree::sphere> *points) {
          std::vector<KdTree::sphere> LevelPoints;
          std::vector<KdTree::sphere> diff;
          std::vector<KdTree::sphere>::iterator it;
          for (auto a : *points) {              
              if (BoxIntersectsSphere(min,max,a.pos,a.radius))
              {
                  LevelPoints.push_back(a);
                  for (it = (*points).begin(); it != (*points).end();) {
                      if ((*it).pos == a.pos)
                          it = (*points).erase(it);
                      else
                      {
                          ++it;
                      }
                  }
              }

          }
          
          return LevelPoints;
      }

      void init() {
          add_input("points", typeid(PointCollection));
          add_input("radii", typeid(vec1f));
          
          //add_input("indice", typeid(vec1i));
          add_input("indice", typeid(std::vector <vec1i>));
          add_output("KDTree", typeid(KdTree));
         

      }
      void gui() {

      }

      void process();
  };
  class VisiblePC :public Node 
  {
  public:
      using Node::Node;
      void init()
      {
          //add_input("KDTree", typeid(KdTree));
          add_input("interior_MAT", typeid(PointCollection));
          add_input("interior_radii", typeid(vec1f));
          add_input("original_pc", typeid(PointCollection));
          add_input("viewPoint", typeid(Vector3D));
          
          add_output("visible_pc", typeid(PointCollection));
      }
      void process();      

      static void GetVisblePT(std::vector<arr3f> pc, PointCollection interior_MAT,vec1f radii,Vector3D viewpoint,PointCollection &visible_pc)
      {
          for (int j = 0; j < pc.size(); j++)
          {
              bool visflag = true;
              Vector3D v2(pc[j][0], pc[j][1], pc[j][2]);
              for (int i = 0; i < interior_MAT.size(); i++)
              {
                  if (i != j) {
                      Vector3D centre(interior_MAT[i][0], interior_MAT[i][1], interior_MAT[i][2]);

                      float dis = DistancePointToSegment(viewpoint, v2, centre);
                      if (dis < radii[i])
                      {
                          visflag = false;
                          break;
                      }
                  }
                  else
                      continue;
              }

              if (visflag == false)continue;
              mtx2.lock();
              visible_pc.push_back({ pc[j][0], pc[j][1], pc[j][2] });
              mtx2.unlock();
          }

      }

      static float PointToPointDis(Vector3D p1, Vector3D p2) 
      {
          float dis = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
          return dis;
      }

      static inline float dotProduct(const Vector3D &a, const Vector3D &b)
      {
          return(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
      }

      static float DistancePointToSegment(Vector3D v1, Vector3D v2, Vector3D centre )
      {
          //v1 viewpoint v2 pt in pc;
     

          Vector3D v_12 = v2 - v1;
          Vector3D v_1c = centre - v1;
          float f = dotProduct(v_12 ,v_1c);
          if (f < 0)
              return PointToPointDis(v1, centre);

          float d = dotProduct(v_12, v_12);
          if( f>=d)
              return PointToPointDis(v2, centre);

          f = f / d;
          Vector3D D = v1 + f * v_12;
          return PointToPointDis(centre, D);

      }

      //v1 viewpoint v2 target point//
      static bool GetOneLineResult(Vector3D v1, Vector3D v2, KdTree *kd)
      {         
          bool visflag = true;
          Vector3D hit;
          int count = 0;
          for (int i = 0; i < (*kd).m_maxpoint.size(); i++) {
              bool a = CheckLineBox((*kd).m_minpoint[i], (*kd).m_maxpoint[i], v1, v2, hit);
              if (a == 1) 
              {
                  for (auto pt : (*kd).m_levelpoints[(*kd).m_maxpoint.size() - i - 1]) {
                      count++;
                      float dis = DistanceOfPointToLine(v1, v2, pt.pos);
                      if (dis < pt.radius) 
                      {
                          visflag = false;
                          break;
                      }
                  }
                  if (visflag == false) 
                  {
                      break;
                  }
              }              
          }
          return visflag;                         

      };
      // overload //
      static float DistanceOfPointToLine(Vector3D a, Vector3D b, Vector3D s)
      {
          float ab = sqrt(pow((a.x - b.x), 2.0f) + pow((a.y - b.y), 2.0f) + pow((a.z - b.z), 2.0f));
          float as = sqrt(pow((a.x - s.x), 2.0f) + pow((a.y - s.y), 2.0f) + pow((a.z - s.z), 2.0f));
          float bs = sqrt(pow((s.x - b.x), 2.0f) + pow((s.y - b.y), 2.0f) + pow((s.z - b.z), 2.0f));
          float cos_A = (pow(as, 2.0f) + pow(ab, 2.0f) - pow(bs, 2.0f)) / (2.0f * ab*as);
          float sin_A = sqrt(1.0f - pow(cos_A, 2.0f));
          return as * sin_A;
      }
      static int inline GetIntersection(float fDst1, float fDst2, Vector3D P1, Vector3D P2, Vector3D &Hit)
      {
          if ((fDst1 * fDst2) >= 0.0f) return 0;
          if (fDst1 == fDst2) return 0;
          Hit = P1 + (P2 - P1) * (-fDst1 / (fDst2 - fDst1));
          return 1;
      }

      static int inline InBox(Vector3D Hit, Vector3D B1, Vector3D B2, const int Axis)
      {
          if (Axis == 1 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          if (Axis == 2 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[0] > B1[0] && Hit[0] < B2[0]) return 1;
          if (Axis == 3 && Hit[0] > B1[0] && Hit[0] < B2[0] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          return 0;
      }
      //minbounding maxbounding viewpoint target point hit
      static int CheckLineBox(Vector3D B1, Vector3D B2, Vector3D L1, Vector3D L2, Vector3D &Hit)
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
      
      

  };
  class VisibiltyQurey :public Node //old version cylinder check
  {
  public: 
      using Node::Node;
      void init() {
          add_input("KDTree", typeid(KdTree));
          add_input("interior_MAT", typeid(PointCollection));
          add_input("ViewPoint", typeid(Vector3D));
          add_input("interior_radii", typeid(vec1f));
          add_input("original_pc", typeid(PointCollection));
          add_output("Visible_MAT", typeid(PointCollection));
          add_output("Radii_of_MAT", typeid(vec1f));
          //add_output("indices", typeid(vec1i));
      }
      void process();
      static float DistanceOfPointToLine(Vector3D a, Vector3D b, Vector3D s)
      {
          //concurrency::precise_math::
          float ab = sqrt(pow((a.x - b.x), 2.0f) + pow((a.y - b.y), 2.0f) + pow((a.z - b.z), 2.0f));
          float as = sqrt(pow((a.x - s.x), 2.0f) +pow((a.y - s.y), 2.0f) + pow((a.z - s.z), 2.0f));
          float bs = sqrt(pow((s.x - b.x), 2.0f) + pow((s.y - b.y), 2.0f) +pow((s.z - b.z), 2.0f));
          float cos_A = (pow(as, 2.0f) + pow(ab, 2.0f) - pow(bs, 2.0f)) / (2.0f * ab*as);
          float sin_A = sqrt(1.0f - pow(cos_A, 2.0f));
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
  class OffsetNode:public Node {
  public:
      using Node::Node;
      void init() {
          add_output("Numfloat", typeid(float));
          add_param("float_value", (float)0.5);
      }
      void gui() {
          ImGui::InputFloat("Num", &param<float>("float_value"));
      }
      void process() {
          float num = param<float>("float_value");
          output("Numfloat").set(num);
      }

  };
  class ParallelVector :public Node 
  {
  public:
      using Node::Node;
      void init(){
          add_input("vector1", typeid(Vector3D));
          add_input("vector2", typeid(Vector3D));
          add_output("Headvectors", typeid(std::vector<Vector3D>));
          add_output("Endvectors", typeid(std::vector<Vector3D>));

          
          add_param("Density", (float)500);
          add_param("Radius", (float)500);
      }
      void gui(){
          /*ImGui::InputFloat("x value", &param<float>("x_value"));
          ImGui::InputFloat("y value", &param<float>("y_value"));
          ImGui::InputFloat("z value", &param<float>("z_value"));*/
          ImGui::InputFloat("Density", &param<float>("Density"));
          ImGui::InputFloat("Radius", &param<float>("Radius"));
          

      }
      void process();
      static inline Vector3D crossProduct(const Vector3D &a, const Vector3D &b)
      {
          Vector3D result;

          result[0] = a[1] * b[2] - a[2] * b[1];
          result[1] = a[2] * b[0] - a[0] * b[2];
          result[2] = a[0] * b[1] - a[1] * b[0];

          return(result);
      }
      
      
  };
  class ViewPoint :public Node {
  public:
      using Node::Node;

      void init() {
          add_output("vector", typeid(Vector3D));
          add_output("point", typeid(PointCollection));

          add_param("x_value", (float)5.0);
          add_param("y_value", (float)0.0);
          add_param("z_value", (float)5.0);
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
                                              
          ///////////////////
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

      
      geoflow::TriangleCollection StandardSphere(geoflow::arr3f center, float radius) 
      {

          //// 10
          //// [11][21]
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
  class GetRaysResult : public Node 
  {
  public:
      using Node::Node;
      void init() 
      {
          add_input("Headvectors", typeid(std::vector<Vector3D>));
          add_input("Endvectors", typeid(std::vector<Vector3D>));
          add_input("KDTree", typeid(KdTree));
          add_input("MATpoints", typeid(PointCollection));
          add_input("radii", typeid(vec1f));


          add_output("MAT_points", typeid(PointCollection));
          add_output("radii", typeid(vec1f));
          add_output("indice", typeid(vec1i));
          
      }
      void process();
      //--------------overload-------------------//
      static float PointToPointDis(Vector3D p1, Vector3D p2) {
          float dis = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
          return dis;
      }

      static int inline GetIntersection(float fDst1, float fDst2, Vector3D P1, Vector3D P2, Vector3D &Hit)
      {
          if ((fDst1 * fDst2) >= 0.0f) return 0;
          if (fDst1 == fDst2) return 0;
          Hit = P1 + (P2 - P1) * (-fDst1 / (fDst2 - fDst1));
          return 1;
      }

      static int inline InBox(Vector3D Hit, Vector3D B1, Vector3D B2, const int Axis)
      {
          if (Axis == 1 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          if (Axis == 2 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[0] > B1[0] && Hit[0] < B2[0]) return 1;
          if (Axis == 3 && Hit[0] > B1[0] && Hit[0] < B2[0] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          return 0;
      }
      static int CheckLineBox(Vector3D B1, Vector3D B2, Vector3D L1, Vector3D L2, Vector3D &Hit)
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
      
  };
 
  class AMPGPUQueryTest :public Node {
  public:
      using Node::Node;
      
      struct  sphere 
      {
          Vector3D pos;
          float r;
      };
      void init() {
          add_input("KDTree", typeid(KdTree));
          add_input("MATpoints", typeid(PointCollection));
          add_input("radii", typeid(vec1f));
          //add_input("indice", typeid(vec1i));
          //add_input("indice", typeid(std::vector<vec1i>));
          add_input("Vector1", typeid(Vector3D));
          add_input("interval", typeid(float));

          add_output("MAT_points", typeid(PointCollection));
          add_output("radii", typeid(vec1f));
          add_output("indice", typeid(vec1i));
      }

      void process();

      static std::vector<Vector3D> SpherePoints(Vector3D v1, float radius, float interval) {
          std::vector<Vector3D> points;
          Vector3D point;
          int nLongitude = (int)interval;
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
      
      static float PointToPointDis(Vector3DNew p1, Vector3DNew p2) restrict(amp,cpu)
      {
          float dis = concurrency::precise_math::sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
          return dis;
      }

      static int inline GetIntersection(float fDst1, float fDst2, Vector3DNew P1, Vector3DNew P2, Vector3DNew &Hit)  restrict(cpu, amp)
      {
          
          if ((fDst1 * fDst2) >= 0.0f) return 0;
          if (fDst1 == fDst2) return 0;
          Hit.x = P1.x + (P2.x - P1.x) * (-fDst1 / (fDst2 - fDst1));
          Hit.y = P1.y + (P2.y - P1.y) * (-fDst1 / (fDst2 - fDst1));
          Hit.z = P1.z + (P2.z - P1.z) * (-fDst1 / (fDst2 - fDst1));
          return 1;
      }

      static int inline InBox (Vector3DNew Hit, Vector3DNew B1, Vector3DNew B2, const int Axis)  restrict(cpu, amp)
      {
          if (Axis == 1 && Hit.z > B1.z && Hit.z < B2.z && Hit.y > B1.y && Hit.y < B2.y) return 1;
          if (Axis == 2 && Hit.z > B1.z && Hit.z < B2.z && Hit.x > B1.x && Hit.x < B2.x) return 1;
          if (Axis == 3 && Hit.x > B1.x && Hit.x < B2.x && Hit.y > B1.y && Hit.y < B2.y) return 1;
          return 0;
      }
      //minbounding maxbounding viewpoint target point hit
      static int CheckLineBox(Vector3DNew B1, Vector3DNew B2, Vector3DNew L1, Vector3DNew L2, Vector3DNew &Hit)  restrict(cpu, amp)
      {
          if (L2.x < B1.x && L1.x < B1.x) return false;
          if (L2.x > B2.x && L1.x > B2.x) return false;
          if (L2.y < B1.y && L1.y < B1.y) return false;
          if (L2.y > B2.y && L1.y > B2.y) return false;
          if (L2.z < B1.z && L1.z < B1.z) return false;
          if (L2.z > B2.z && L1.z > B2.z) return false;
          if (L1.x > B1.x && L1.x < B2.x &&
              L1.y > B1.y && L1.y < B2.y &&
              L1.z > B1.z && L1.z < B2.z)
          {
              Hit = L1;
              return 1;
          }
          if (GetIntersection(L1.x - B1.x, L2.x - B1.x, L1, L2, Hit)&&(InBox(Hit, B1, B2, 1))
              || (GetIntersection(L1.y - B1.y, L2.y - B1.y, L1, L2, Hit) && InBox(Hit, B1, B2, 2))
              || (GetIntersection(L1.z - B1.z, L2.z - B1.z, L1, L2, Hit) && InBox(Hit, B1, B2, 3))
              || (GetIntersection(L1.x - B2.x, L2.x - B2.x, L1, L2, Hit) && InBox(Hit, B1, B2, 1))
              || (GetIntersection(L1.y - B2.y, L2.y - B2.y, L1, L2, Hit) && InBox(Hit, B1, B2, 2))
              || (GetIntersection(L1.z - B2.z, L2.z - B2.z, L1, L2, Hit) && InBox(Hit, B1, B2, 3)))
              return true;   

          return false;
      }
      static float DistanceOfPointToLine(Vector3DNew a, Vector3DNew b, Vector3DNew s)restrict(amp,cpu)
      {
          //concurrency::precise_math::
          float ab = concurrency::precise_math::sqrt(concurrency::precise_math::pow((a.x - b.x), 2.0f) + concurrency::precise_math::pow((a.y - b.y), 2.0f) + concurrency::precise_math::pow((a.z - b.z), 2.0f));
          float as = concurrency::precise_math::sqrt(concurrency::precise_math::pow((a.x - s.x), 2.0f) + concurrency::precise_math::pow((a.y - s.y), 2.0f) + concurrency::precise_math::pow((a.z - s.z), 2.0f));
          float bs = concurrency::precise_math::sqrt(concurrency::precise_math::pow((s.x - b.x), 2.0f) + concurrency::precise_math::pow((s.y - b.y), 2.0f) + concurrency::precise_math::pow((s.z - b.z), 2.0f));
          float cos_A = (concurrency::precise_math::pow(as, 2.0f) + concurrency::precise_math::pow(ab, 2.0f) - concurrency::precise_math::pow(bs, 2.0f)) / (2.0f * ab*as);
          float sin_A = concurrency::precise_math::sqrt(1.0f - concurrency::precise_math::pow(cos_A, 2.0f));
          return as * sin_A;
      }
      
      

      static void GPUQuery(std::vector<Vector3DNew> v2, Vector3DNew v1,int size, KdTree *kd,  std::vector<int> &IfIntersect)
      {
          int levelsize = (*kd).m_maxpoint.size();

          std::vector<Vector3DNew> new_m_maxpoint;
          std::vector<Vector3DNew> new_m_minpoint;
          for (auto a : (*kd).m_maxpoint) 
          {
              Vector3DNew new_a(a);
              new_m_maxpoint.push_back(new_a);
          }

          for (auto a2 : (*kd).m_minpoint)
          {
              Vector3DNew new_a2(a2);
              new_m_minpoint.push_back(new_a2);
          }

          //v1 is the viewpoint
          concurrency::array_view<Vector3DNew, 1> v1_av(size);
          for (int i = 0; i < size; i++) 
          {
              v1_av[i] = v1;
          }
          
          
          // bv is the vector of input target
          concurrency::array_view<Vector3DNew,1> v2_av(size, v2);
          
          concurrency::array_view<Vector3DNew,1> pointlist(levelsize);
          concurrency::array_view<float,1> radiilist(levelsize);

          concurrency::array_view<Vector3DNew,1> kd_maxpoint(levelsize, new_m_maxpoint);
          concurrency::array_view<Vector3DNew,1> kd_minpoint(levelsize, new_m_minpoint);
          concurrency::array_view<Vector3DNew, 1> hit_av(size);


          //concurrency::array_view<KdTree::sphere, 1> kd_levelpoints(levelsize,50);
          
          //concurrency::array_view<int> kd_levelpoint_size(levelsize);
          /*for (int i = 0; i < levelsize; i++) 
              for (int j = 0; j < 50; j++) 
              {
                  kd_levelpoints[i *
              }*/
          /*for(int i =0;i<levelsize;i++)
          {
              for (int j = 0; j < (*kd).m_levelpoints[i].size();j++) 
              {
                  kd_levelpoints[i * 50 + j] = (*kd).m_levelpoints[i][j];
              }
          }*/

          

          //for (int i = 0; i < levelsize; i++) 
          //{
          //    kd_levelpoint_size[i] = (*kd).m_levelpoints[i].size();
          //}


          //result bool if intersect
          concurrency::array_view<int, 1> result_av(size);
          for (int i = 0; i < size; i++) 
          {
              result_av[i] = 0;
          }
          

          
          concurrency::parallel_for_each(v2_av.extent,
              [v1_av, v2_av, kd_maxpoint, kd_minpoint,levelsize, hit_av, result_av](concurrency::index<1> idx) restrict(amp)
          {    
              int flag = 0;          
              /*for (int i = 0; i < levelsize; i++) 
              {
                  
                  flag = AMPGPUQueryTest::CheckLineBox(kd_minpoint[i],kd_maxpoint[i],v1_av[idx],v2_av[idx],hit_av[idx]);
                  if (flag == 1) 
                  {   
                      result_av[idx] = 1;     
                      break;
                  }     
              }*/
              flag = AMPGPUQueryTest::CheckLineBox(kd_minpoint[0], kd_maxpoint[0], v1_av[idx], v2_av[idx], hit_av[idx]);
              if (flag == 1)
              {
                  result_av[idx] = 1;
              }

          });
          result_av.synchronize();
          

          for (int i = 0; i < size; i++) 
          {
              IfIntersect.push_back(result_av[i]);
          }
          
          
      }
      static void AMPGPUQueryTest::GetQueryResult(std::vector<Vector3D> vecList, Vector3D v1, KdTree *kd, std::vector<KdTree::sphere> &visble_sph) 
      {
          for (Vector3D v2 : vecList)
          {
              KdTree::sphere result;
              std::vector<Vector3D> pointlist;
              std::vector<float> radiilist;
              Vector3D hit;
              int count = 0;
              for (int i = 0; i < (*kd).m_maxpoint.size(); i++) {
                  bool a = AMPGPUQueryTest::CheckLineBox((*kd).m_minpoint[i], (*kd).m_maxpoint[i], v1, v2, hit);
                  if (a == 1) {
                      for (auto pt : (*kd).m_levelpoints[(*kd).m_maxpoint.size() - i - 1]) {
                          count++;
                          float dis = VisibiltyQurey::DistanceOfPointToLine(v1, v2, pt.pos);
                          if (dis <= pt.radius) {
                              pointlist.push_back(pt.pos);
                              radiilist.push_back(pt.radius);
                          }
                      }
                  }
              }
              if (pointlist.size() > 0) {
                  //std::cout << "----------------this direction has :" << pointlist.size() << "  intersected" << std::endl;
                  //float minDis = MutiThreadsOneQuery::PointToPointDis(v1, pointlist[0]);
                  float minDis = AMPGPUQueryTest::PointToPointDis(v1, pointlist[0]) - radiilist[0];

                  int flag = 0;
                  for (int i = 0; i < pointlist.size(); i++)
                  {
                      float temp = AMPGPUQueryTest::PointToPointDis(v1, pointlist[i]) - radiilist[i];
                      //float temp = MutiThreadsOneQuery::PointToPointDis(v1, pointlist[i]);
                      if (temp < minDis) {
                          minDis = temp;
                          flag = i;
                      }
                  }
                  result.pos = { pointlist[flag].x,pointlist[flag].y,pointlist[flag].z };
                  result.radius = radiilist[flag];
                  //std::cout << "result :" << result.pos.x<<"," << result.pos.y<<","<< result.pos.z<< std::endl;
                  mtx.lock();
                  visble_sph.push_back(result);
                  mtx.unlock();
              }
          }

      };

      static KdTree::sphere GetOneLineResult(Vector3D v1, Vector3D v2,KdTree *kd)
      {                       
          KdTree::sphere result;
          std::vector<Vector3D> pointlist;
          std::vector<float> radiilist;
          //std::vector<int> indice;
          std::vector<std::vector<int>> indice;
          Vector3DNew hit;
          int count = 0;
          //--------------------------allpoints-----------------------//          
          bool a = AMPGPUQueryTest::CheckLineBox((*kd).m_minpoint[0], (*kd).m_maxpoint[0], v1, v2, hit);
          if (a == 1) {

              for (auto pt : (*kd).m_allballs)
              {
                  count++;
                  float dis = VisibiltyQurey::DistanceOfPointToLine(v1, v2, pt.pos);
                  if (dis <= pt.radius) {
                      pointlist.push_back(pt.pos);
                      radiilist.push_back(pt.radius);
                      indice.push_back(pt.index);
                  }
              }
          }          
          //------------------------levelpoints--------------------------//
          /*for (int i = 0; i < (*kd).m_maxpoint.size(); i++) {
              bool a = AMPGPUQueryTest::CheckLineBox((*kd).m_minpoint[i], (*kd).m_maxpoint[i], v1, v2, hit);
              if (a == 1) {

                  
                  for (auto pt : (*kd).m_levelpoints[(*kd).m_maxpoint.size() - i - 1])                   
                  {
                      count++;
                      float dis = VisibiltyQurey::DistanceOfPointToLine(v1, v2, pt.pos);
                      if (dis <= pt.radius) {
                          pointlist.push_back(pt.pos);
                          radiilist.push_back(pt.radius);
                          indice.push_back(pt.index);
                      }
                  }
              }
          }*/
          //------------------------------------------------//
          if (pointlist.size() > 0) {
              //std::cout << "----------------this direction has :" << pointlist.size() << "  intersected" << std::endl;
              //float minDis = MutiThreadsOneQuery::PointToPointDis(v1, pointlist[0]);
              float minDis = AMPGPUQueryTest::PointToPointDis(v1, pointlist[0]) - radiilist[0];

              int flag = 0;
              for (int i = 0; i < pointlist.size(); i++)
              {
                  float temp = AMPGPUQueryTest::PointToPointDis(v1, pointlist[i]) - radiilist[i];
                  //float temp = MutiThreadsOneQuery::PointToPointDis(v1, pointlist[i]);
                  if (temp < minDis) {
                      minDis = temp;
                      flag = i;
                  }
              }
              result.pos = { pointlist[flag].x,pointlist[flag].y,pointlist[flag].z };
              result.radius = radiilist[flag];
              result.index = indice[flag];             
              
          }
          return result;

      };
      
      //--------------------   overload -----------------------//
      static float PointToPointDis(Vector3D p1, Vector3D p2) {
          float dis = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
          return dis;
      }

      static int inline GetIntersection(float fDst1, float fDst2, Vector3D P1, Vector3D P2, Vector3D &Hit)
      {
          if ((fDst1 * fDst2) >= 0.0f) return 0;
          if (fDst1 == fDst2) return 0;
          Hit = P1 + (P2 - P1) * (-fDst1 / (fDst2 - fDst1));
          return 1;
      }

      static int inline InBox(Vector3D Hit, Vector3D B1, Vector3D B2, const int Axis)
      {
          if (Axis == 1 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          if (Axis == 2 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[0] > B1[0] && Hit[0] < B2[0]) return 1;
          if (Axis == 3 && Hit[0] > B1[0] && Hit[0] < B2[0] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          return 0;
      }
      static int CheckLineBox(Vector3D B1, Vector3D B2, Vector3D L1, Vector3D L2, Vector3D &Hit)
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
  };
  
  class MutiThreadsOneQuery :public Node {
  public:
      using Node::Node;
      
      
      struct  sphere
      {
          Vector3D pos;
          float r;
      };
      //std::vector<sphere> m_visble_sph;
      
      void init() {
          add_input("KDTree", typeid(KdTree));
          add_input("MATpoints", typeid(PointCollection));
          add_input("radii", typeid(vec1f));
          add_input("Vector1", typeid(Vector3D));

          add_output("MAT_points", typeid(PointCollection));
          add_output("radii", typeid(vec1f));
      }
      static std::vector<std::vector<Vector3D>> CutVecList(std::vector<Vector3D> vecList, int CutNum) {
          std::vector<std::vector<Vector3D>> OutVecList;
          float interval = 1.0 / CutNum;
          float interval_size = interval * vecList.size();;
          
          for (int i = 0; i < CutNum; i++) {
              std::vector<Vector3D> TemVecList;
              
              std::for_each(begin(vecList)+(i*interval_size), begin(vecList) + ((i+1)*interval_size), [&TemVecList](Vector3D x) {
                  TemVecList.push_back(x);
              }); 

              OutVecList.push_back(TemVecList);       
          
          }
          return OutVecList;
      }
      
      static float PointToPointDis(Vector3D p1, Vector3D p2) {
          float dis = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
          return dis;
      }

      static int inline GetIntersection(float fDst1, float fDst2, Vector3D P1, Vector3D P2, Vector3D &Hit) 
      {
          if ((fDst1 * fDst2) >= 0.0f) return 0;
          if (fDst1 == fDst2) return 0;
          Hit = P1 + (P2 - P1) * (-fDst1 / (fDst2 - fDst1));
          return 1;
      }

      static int inline InBox(Vector3D Hit, Vector3D B1, Vector3D B2, const int Axis)  
      {
          if (Axis == 1 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          if (Axis == 2 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[0] > B1[0] && Hit[0] < B2[0]) return 1;
          if (Axis == 3 && Hit[0] > B1[0] && Hit[0] < B2[0] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          return 0;
      }
      static int CheckLineBox(Vector3D B1, Vector3D B2, Vector3D L1, Vector3D L2, Vector3D &Hit) 
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
     
      static void GetQueryResult(std::vector<Vector3D> vecList, Vector3D v1, KdTree *kd, std::vector<MutiThreadsOneQuery::sphere> &visble_sph) {
          
          
          for (Vector3D v2 : vecList) 
          {
              MutiThreadsOneQuery::sphere result;
              std::vector<Vector3D> pointlist;
              std::vector<float> radiilist;
              Vector3D hit;
              int count = 0;
              for (int i = 0; i < (*kd).m_maxpoint.size(); i++) {
                  bool a = MutiThreadsOneQuery::CheckLineBox((*kd).m_minpoint[i], (*kd).m_maxpoint[i], v1, v2, hit);
                  if (a == 1) {                      
                      for (auto pt : (*kd).m_levelpoints[(*kd).m_maxpoint.size() - i - 1]) {
                          count++;                          
                          float dis = VisibiltyQurey::DistanceOfPointToLine(v1, v2, pt.pos);                       
                          if (dis <= pt.radius) {                          
                              pointlist.push_back(pt.pos);
                              radiilist.push_back(pt.radius);                              
                          }
                      }                      
                  }
              }
              if (pointlist.size() > 0) {
                  //std::cout << "----------------this direction has :" << pointlist.size() << "  intersected" << std::endl;
                  //float minDis = MutiThreadsOneQuery::PointToPointDis(v1, pointlist[0]);
                  float minDis = MutiThreadsOneQuery::PointToPointDis(v1, pointlist[0]) -radiilist[0];
                  
                  int flag = 0;
                  for (int i = 0; i < pointlist.size(); i++)
                  {
                      float temp = MutiThreadsOneQuery::PointToPointDis(v1, pointlist[i]) - radiilist[i];
                      //float temp = MutiThreadsOneQuery::PointToPointDis(v1, pointlist[i]);
                      if (temp < minDis) {
                          minDis = temp;
                          flag = i;
                      }
                  }
                  result.pos = { pointlist[flag].x,pointlist[flag].y,pointlist[flag].z };
                  result.r = radiilist[flag];
                  //std::cout << "result :" << result.pos.x<<"," << result.pos.y<<","<< result.pos.z<< std::endl;
                  mtx.lock();
                  visble_sph.push_back(result);
                  mtx.unlock();
              }
            }
          
      };
      void process();
  };
  class VisiblePart :public Node 
  {
  public:
      using Node::Node;
      void init() 
      {
          add_input("viewpoint", typeid(Vector3D));
          add_input("targetPC", typeid(PointCollection));
          add_input("KDTree", typeid(KdTree));
          add_input("MATpoints", typeid(PointCollection));
          add_input("radii", typeid(vec1f));

          add_output("visible_parts", typeid(PointCollection));          
      }
      void process();
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

      static int inline GetIntersection(float fDst1, float fDst2, Vector3D P1, Vector3D P2, Vector3D &Hit) {
          if ((fDst1 * fDst2) >= 0.0f) return 0;
          if (fDst1 == fDst2) return 0;
          Hit = P1 + (P2 - P1) * (-fDst1 / (fDst2 - fDst1));
          return 1;
      }

      static int inline InBox(Vector3D Hit, Vector3D B1, Vector3D B2, const int Axis) {
          if (Axis == 1 && Hit[2] > B1[2] && Hit[2] < B2[2]&& Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          if (Axis == 2 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[0] > B1[0] && Hit[0] < B2[0]) return 1;
          if (Axis == 3 && Hit[0] > B1[0] && Hit[0] < B2[0] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
          return 0;
      }
      static int CheckLineBox(Vector3D B1, Vector3D B2, Vector3D L1, Vector3D L2, Vector3D &Hit)
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
