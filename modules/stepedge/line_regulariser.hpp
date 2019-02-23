#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include<geoflow/core/geoflow.hpp>

namespace linereg {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_2 Point_2;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 Point_3;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Vector_2 Vector_2;

  typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
  typedef CGAL::Polygon_2<EK> Polygon_2;

  class LineRegulariser {
    static constexpr double pi = 3.14159265358979323846;

    struct ValueCluster {
      Vector_2 ref_vec;
      Vector_2 ref_point;
      std::vector<size_t> idx;
    };

    typedef std::tuple<double, Point_2, double, size_t, size_t, double> linetype; 
      // new angle, midpoint, distance in angle cluster, priority, segment_id, sqlength
    typedef std::vector<EK::Segment_2> vec_ek_seg;

    // geoflow::SegmentCollection& input_segments;
    public:
    std::vector<linetype> lines;
    // vec_ek_seg input_reg_exact;
    double angle_threshold, dist_threshold;

    std::unordered_map<size_t, vec_ek_seg> segments;

    LineRegulariser() {};

    void add_segments(size_t priority, geoflow::SegmentCollection& segs) {
      size_t i=0;
      segments[priority] = vec_ek_seg();
      for(auto& edge : segs) {
        auto source = EK::Point_2(edge[0][0], edge[0][1]);
        auto target = EK::Point_2(edge[1][0], edge[1][1]);
        auto v = target-source;
        auto p_ = source + v/2;
        auto p = Point_2(CGAL::to_double(p_.x()),CGAL::to_double(p_.y()));
        auto l = CGAL::to_double(v.squared_length());
        auto angle = std::atan2(CGAL::to_double(v.x()), CGAL::to_double(v.y()));
        if (angle < 0) angle += pi;
        lines.push_back(std::make_tuple(angle,p,0,priority,i++,l));
        segments[priority].push_back(EK::Segment_2(source, target));
      }
    }

    vec_ek_seg& get_segments(const size_t& priority) {
      return segments[priority];
    }

    // template<typename element> void get_max(const std::vector<size_t>& idx, const size_t element) {
      
    //   for(auto& i : idx) {
    //     std::get<element>(lines[i]);
        
    //   }
    // }

    void cluster() {
      // cluster by angle
      std::vector<size_t> edge_idx(lines.size());
      for (size_t i=0; i<lines.size(); ++i) {
        edge_idx[i]=i;
      }
      std::sort(edge_idx.begin(), edge_idx.end(), [&lines=lines](size_t a, size_t b) {
        return std::get<0>(lines[a]) < std::get<0>(lines[b]);   
      });

      std::vector<ValueCluster> angle_clusters(1);
      auto angle_sum = std::get<0>(lines[edge_idx[0]]);
      angle_clusters.back().idx.push_back(edge_idx[0]);
      for(size_t i=1; i < edge_idx.size(); ++i) {
        auto& edge_id = edge_idx[i];
        auto& line = lines[edge_id];
        auto& angle = std::get<0>(line);
        auto mean_angle = angle_sum/angle_clusters.back().idx.size();
        if(std::fmod((angle - mean_angle), pi) < angle_threshold) {
          angle_clusters.back().idx.push_back(edge_id);
          angle_sum += angle;
        } else {
          angle_clusters.resize(angle_clusters.size()+1);
          angle_clusters.back().idx.push_back(edge_id);
          angle_sum = angle;
        }
      }

      // determine angle for each cluster
      for(auto& cluster : angle_clusters) {
        // find maximum priority in this cluster
        size_t max_pr=0;
        for(auto& i : cluster.idx)
          max_pr = std::max(max_pr, std::get<3>(lines[i]));

        // collect all ids that have the max priority
        std::vector<size_t> max_idx;
        for(auto& i : cluster.idx)
          if (std::get<3>(lines[i]) == max_pr)
            max_idx.push_back(i);

        // computed average angle weighted by segment length
        double sum_len=0;
        double sum_all=0;
        std::vector<double> weighted_angles;
        for(auto& i : max_idx) {
          auto& len = std::get<5>(lines[i]);
          auto& angle = std::get<0>(lines[i]);
          sum_all += len * angle;
          sum_len += len;
        }
        double angle = sum_all/sum_len;
        Vector_2 n(-1.0, std::tan(angle));
        cluster.ref_vec = n/std::sqrt(n.squared_length()); // normalize
        
        // or median angle:
        // size_t median_id = cluster.idx[cluster.idx.size()/2];
        // cluster.value = std::get<0>(lines[median_id]);
      }

      // cluster parallel lines by distance
      std::vector<ValueCluster> dist_clusters;
      for(auto& cluster : angle_clusters) {
        auto n = cluster.ref_vec;
        // compute distances along n wrt to first line in cluster
        auto p = std::get<1>(lines[cluster.idx[0]]);
        for(auto& i : cluster.idx) {
          auto q = std::get<1>(lines[i]);
          auto v = p-q;
          std::get<2>(lines[i]) = v*n;
          // distances.push_back(v*n);
        }
        // sort by distance, ascending
        auto dist_idx = cluster.idx;
        std::sort(dist_idx.begin(), dist_idx.end(), [&lines=lines](size_t a, size_t b){
            return std::get<2>(lines[a]) < std::get<2>(lines[b]);
        });
        // cluster nearby lines using separation threshold
        double dist_sum = std::get<2>(lines[dist_idx[0]]);
        dist_clusters.resize(dist_clusters.size()+1);
        dist_clusters.back().ref_vec = n;
        dist_clusters.back().idx.push_back(dist_idx[0]);
        for(size_t i=1; i < dist_idx.size(); ++i) {
          auto& edge_id = dist_idx[i];
          auto& line = lines[edge_id];
          double mean_dist = dist_sum / dist_clusters.back().idx.size();
          double& dist = std::get<2>(line);
          if (std::abs(mean_dist-dist) < dist_threshold) {
            dist_clusters.back().idx.push_back(edge_id);
            dist_sum += dist;
          } else {
            dist_clusters.resize(dist_clusters.size()+1);
            dist_clusters.back().ref_vec = n;
            dist_clusters.back().idx.push_back(edge_id);
            dist_sum = dist;
          }
        }
      }

      // find average direction and center point for each cluster
      for(auto& cluster : dist_clusters) {
        // find max priority
        size_t max_pr=0;
        for(auto& i : cluster.idx)
          max_pr = std::max(max_pr, std::get<3>(lines[i]));

        // collect all ids that have the max priority
        std::vector<size_t> max_idx;
        for(auto& i : cluster.idx)
          if (std::get<3>(lines[i]) == max_pr)
            max_idx.push_back(i);

        // computed average angle weighted by segment length
        double sum_len=0;
        Vector_2 sum_all(0,0);
        std::vector<double> weighted_angles;
        for(auto& i : max_idx) {
          auto& len = std::get<5>(lines[i]);
          auto& q = std::get<1>(lines[i]);
          sum_all += len*Vector_2(q.x(), q.y());
          sum_len += len;
        }
        // cluster.distance = sum/cluster.idx.size();
        cluster.ref_point = sum_all/sum_len;
        cluster.ref_vec = Vector_2(cluster.ref_vec.y(), -cluster.ref_vec.x());
      }

      // project input segments on the cluster lines
      for(auto& cluster : dist_clusters) {
        auto ref_v = EK::Vector_2(cluster.ref_vec.x(), cluster.ref_vec.y());
        auto ref_p = EK::Point_2(cluster.ref_point.x(), cluster.ref_point.y());

        EK::Line_2 ref_line(ref_p, ref_v);

        for(auto& i : cluster.idx) {
          auto pri = std::get<3>(lines[i]);
          auto j = std::get<4>(lines[i]);
          auto& edge = segments[pri][j];
          auto s_new = ref_line.projection(edge.source());
          auto t_new = ref_line.projection(edge.target());
          segments[pri][j] = EK::Segment_2(s_new, t_new);
          // input_segments[i][0] = 
          //   {float(CGAL::to_double(s_new.x())), float(CGAL::to_double(s_new.y())), 0};
          // input_segments[i][1] = 
          //   {float(CGAL::to_double(t_new.x())), float(CGAL::to_double(t_new.y())), 0};
        }
      }
    };

  };

  void chain(const EK::Segment_2& a, const EK::Segment_2& b, Polygon_2& ring, const float& snap_threshold) {
    auto l_a = a.supporting_line();
    auto l_b = b.supporting_line();
    EK::Segment_2 s(a.target(), b.source());
    auto result = CGAL::intersection(l_a, l_b);
    if (result) {
      if (auto p = boost::get<EK::Point_2>(&*result)) {
        if (CGAL::squared_distance(*p, s) < snap_threshold) {
          ring.push_back(*p);
        } else {
          ring.push_back(a.target());
          ring.push_back(b.source());
        }
      }
    // } else if (auto l = boost::get<K::Line_2>(&*result)) {

    // }
    } else { // there is no intersection
      ring.push_back(a.target());
      ring.push_back(b.source());
    }
  }

  // void chain(Segment& a, Segment& b, LinearRing& ring, const float& snap_threshold) {
  Polygon_2 chain_ring(const std::vector<size_t>& idx, const std::vector<EK::Segment_2>& segments, const float& snap_threshold) {
    Polygon_2 ring, fixed_ring;

    if (idx.size()>2) {
      for (size_t i=idx[1]; i<idx[0]+idx.size(); ++i) {
        chain(segments[i-1], segments[i], ring, snap_threshold);
      }
      chain(segments[idx[idx.size()-1]], segments[idx[0]], ring, snap_threshold);

      // get rid of segments with zero length
      // check again the size, to ignore degenerate case of input ring that consists of 3 co-linear segments (would get chained to eg 0 vertices)
      if (ring.size()>2) {
        auto circ = ring.vertices_circulator();
        auto curr = circ;
        do {
          auto d = CGAL::squared_distance(*circ, *(circ+1));
          if (d > 1E-6) {
            fixed_ring.push_back(*circ);
          }
          ++circ;
        } while (curr != circ);
      }
    }

    return fixed_ring;
  }
}