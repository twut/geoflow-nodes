#include "masb_nodes.hpp"
#include "iostream"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <time.h>
#include<amp.h>
#include "region_grower.hpp"
#include "region_grower_testers.hpp"



namespace geoflow::nodes::mat {
    

    void ComputeMedialAxisNode::process() {
        std::cout << "computing MAT " << std::endl;
        auto point_collection = input("points").get<PointCollection>();
        auto normals_vec3f = input("normals").get<vec3f>();

        masb::ma_parameters params;
        params.initial_radius = param<float>("initial_radius");
        params.denoise_preserve = param<float>("denoise_preserve");
        params.denoise_planar = param<float>("denoise_planar");
        params.nan_for_initr = param<bool>("nan_for_initr");

        //std::cout << "where is the error " << std::endl;
        // prepare data structures and transfer data
        masb::ma_data madata;
        madata.m = point_collection.size();

        masb::PointList coords;
        coords.reserve(madata.m);
        for (auto& p : point_collection) {
            coords.push_back(masb::Point(p.data()));
        }

        std::cout << "where is the error " << std::endl;

        masb::VectorList normals;
        normals.reserve(madata.m);
        for (auto& n : normals_vec3f) {
            normals.push_back(masb::Vector(n.data()));
        }
        masb::PointList ma_coords_(madata.m * 2);
        std::vector<int> ma_qidx_(madata.m * 2);

        madata.coords = &coords;
        madata.normals = &normals;
        madata.ma_coords = &ma_coords_;
        madata.ma_qidx = ma_qidx_.data();

        // compute mat points
        masb::compute_masb_points(params, madata);

        // retrieve mat points
        vec1i ma_qidx;
        ma_qidx.reserve(madata.m * 2);
        for (size_t i = 0; i < madata.m * 2; ++i) {
            ma_qidx.push_back(madata.ma_qidx[i]);
        }

        PointCollection ma_coords;
        ma_coords.reserve(madata.m * 2);
        for (auto& c : ma_coords_) {
            ma_coords.push_back({ c[0], c[1], c[2] });
        }

        // Compute medial geometry
        vec1f ma_radii(madata.m * 2);
        vec1f ma_sepangle(madata.m * 2);
        vec3f ma_spoke_f1(madata.m * 2);
        vec3f ma_spoke_f2(madata.m * 2);
        vec3f ma_bisector(madata.m * 2);
        vec3f ma_spokecross(madata.m * 2);
        for (size_t i = 0; i < madata.m * 2; ++i) {
            auto i_ = i % madata.m;
            auto& c = ma_coords_[i];
            // feature points
            auto& f1 = coords[i_];
            auto& f2 = coords[ma_qidx[i]];
            // radius
            ma_radii[i] = Vrui::Geometry::dist(f1, c);
            // spoke vectors
            auto s1 = f1 - c;
            auto s2 = f2 - c;
            ma_spoke_f1[i] = { s1[0], s1[1], s1[2] };
            ma_spoke_f2[i] = { s2[0], s2[1], s2[2] };
            // bisector
            s1.normalize();
            s2.normalize();
            auto b = (s1 + s2).normalize();
            ma_bisector[i] = { b[0], b[1], b[2] };
            // separation angle
            ma_sepangle[i] = std::acos(s1*s2);
            // cross product of spoke vectors
            auto scross = Vrui::Geometry::cross(s1, s2).normalize();
            ma_spokecross[i] = { scross[0], scross[1], scross[2] };
        }
        vec1i ma_is_interior(madata.m * 2, 0);
        std::fill_n(ma_is_interior.begin(), madata.m, 1);

        output("ma_coords").set(ma_coords);
        output("ma_qidx").set(ma_qidx);
        output("ma_radii").set(ma_radii);
        output("ma_is_interior").set(ma_is_interior);
        output("ma_sepangle").set(ma_sepangle);
        output("ma_bisector").set(ma_bisector);
        output("ma_spoke_f1").set(ma_spoke_f1);
        output("ma_spoke_f2").set(ma_spoke_f2);
        output("ma_spokecross").set(ma_spokecross);

        std::cout << "computing MAT done" << std::endl;
    }


    void NegNormalDetector::process() 
    {
        //---------input -------------//
        auto pc = input("originalPC").get<PointCollection>();
        auto normals_vec3f = input("normals").get<vec3f>();
        auto offset = input("offset").get<float>();

        //---------output--------------//
        PointCollection neg_pc;
        vec3f normal_fixed;

        // -----------process -------------//
        for (int i = 0; i < normals_vec3f.size(); i++) 
        {
            if (normals_vec3f[i][2] < offset )
            {
                neg_pc.push_back({ pc[i][0],pc[i][1],pc[i][2] });
                normal_fixed.push_back({ -normals_vec3f[i][0] ,-normals_vec3f[i][1] ,-normals_vec3f[i][2] });
            }
            else
            {
                normal_fixed.push_back({ normals_vec3f[i][0] ,normals_vec3f[i][1] ,normals_vec3f[i][2] });
            }
        }        
        output("pc").set(neg_pc);
        output("normals_fixed").set(normal_fixed);
    }
    

    void MATfilter::process() {
        std::cout << "Filter running" << std::endl;
        //auto pc = input("originalPC").get<PointCollection>();
        auto matpoints = input("ma_coords").get<PointCollection>();
        auto interior_index = input("ma_is_interior").get<vec1i>();
        auto ma_radii = input("ma_radii").get<vec1f>();
        //auto normals_vec3f = input("normals").get<vec3f>();

        //float offset = input("offset").get<float>();
        //float min_z = input("min_z").get<float>();

        
        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\ma_is_interior.txt";
        std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);
        for (auto a : interior_index)
        {
            outfile << a << std::endl;
        }
        outfile.close();

        
        //-----output-----------------//
        PointCollection interior_mat;
        PointCollection exterior_mat;
        vec1f  interior_radii;
        vec1i interior_idx;
        vec1f  exterior_radii;
        vec1i exterior_idx;

        //----------process-----------//

       /* std::string filepath2 = "c:\\users\\tengw\\documents\\git\\Results\\ex_indices_out.txt";
        std::ofstream outfile2(filepath2, std::fstream::out | std::fstream::trunc);

        std::string filepath1 = "c:\\users\\tengw\\documents\\git\\Results\\in_filtered_MAT.txt";
        std::ofstream outfile1(filepath1, std::fstream::out | std::fstream::trunc);*/
        

        for (int i = 0; i < matpoints.size(); i++) 
        {
            if (interior_index[i] == 1 )            
            {
                interior_mat.push_back(matpoints[i]);
                interior_radii.push_back(ma_radii[i]);
                interior_idx.push_back(i);
                //outfile1 << matpoints[i][0] << "," << matpoints[i][1] << "," << matpoints[i][2] << "," << ma_radii[i] << "," << i << std::endl;
                
            }
            else
            {
                exterior_mat.push_back(matpoints[i]);
                exterior_radii.push_back(ma_radii[i]);
                exterior_idx.push_back(i-0.5*matpoints.size());
                //outfile2 << i - 0.5*matpoints.size() << std::endl;
            }            
        }
        
        for (int i = 0; i < interior_mat.size(); i++) 
        {
            //Vector3D in_pt(interior_mat[i][0], interior_mat[i][1], interior_mat[i][2]);
            //Vector3D pc_point(pc[i][0], pc[i][1], pc[i][2]);

            // normal points upwards
            /*float flag = if_interiorMAT(in_pt, pc_point);
            
            if (normals_vec3f[i][2] < 0)
            {
                if (flag < 0 || flag == 0) continue;
                if (flag > 0)
                {
                    auto temp_mat = interior_mat[i];
                    auto temp_radii = interior_radii[i];
                    auto temp_index = interior_idx[i];

                    interior_mat[i] = exterior_mat[i];
                    interior_radii[i] = exterior_radii[i];
                    interior_idx[i] = exterior_idx[i];

                    exterior_mat[i] = temp_mat;
                    exterior_radii[i] = temp_radii;
                    exterior_idx[i] = temp_index;
                }
            }
            if (normals_vec3f[i][2] > 0) 
            {
                if (flag > 0 || flag == 0) continue;
                if(flag<0)
                {
                    auto temp_mat = interior_mat[i];
                    auto temp_radii = interior_radii[i];
                    auto temp_index = interior_idx[i];

                    interior_mat[i] = exterior_mat[i];
                    interior_radii[i] = exterior_radii[i];
                    interior_idx[i] = exterior_idx[i];

                    exterior_mat[i] = temp_mat;
                    exterior_radii[i] = temp_radii;
                    exterior_idx[i] = temp_index;
                }

            }*/

        }

        //outfile1.close();
        //outfile2.close();

        std::cout << "Number of input points:" << matpoints.size() << std::endl;
        std::cout << "Number of interior MAT points :" << interior_mat.size() << std::endl;
        std::cout << "Number of exterior MAT points :" << exterior_mat.size() << std::endl;
        output("interior_mat").set(interior_mat);        
        output("interior_radii").set(interior_radii);
        output("interior_idx").set(interior_idx);

        output("exterior_mat").set(exterior_mat);
        output("exterior_radii").set(exterior_radii);
        output("exterior_idx").set(exterior_idx);
    }
    void RegionGrowMedialAxisNode::process() 
    {
        auto ma_coords = input("ma_coords").get<PointCollection>();
        auto ma_bisector = input("ma_bisector").get<vec3f>();
        auto ma_sepangle = input("ma_sepangle").get<vec1f>();
        auto ma_radii = input("ma_radii").get<vec1f>();

        regiongrower::RegionGrower<MaData, Region> R;
        R.min_segment_count = param<int>("min_count");

        MaData D(ma_coords, ma_bisector, ma_sepangle, ma_radii, param<int>("k"));

        switch (param<int>("method"))
        {
        case 0: {
            AngleOfVectorsTester T_bisector_angle(param<float>("bisector_angle"));
            R.grow_regions(D, T_bisector_angle); break;
        } case 1: {
            DiffOfAnglesTester T_separation_angle(param<float>("separation_angle"));
            R.grow_regions(D, T_separation_angle); break;
        } case 2: {
            BallOverlapTester T_ball_overlap(param<float>("ball_overlap"));
            R.grow_regions(D, T_ball_overlap); break;
        } case 3: {
            CountTester T_shape_count(param<int>("shape_count"));
            R.grow_regions(D, T_shape_count); break;
        } default: break;
        };

        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\segment_id.txt";
        std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);

        vec1i segment_ids;
        for (auto& region_id : R.region_ids) {
            segment_ids.push_back(int(region_id));
        }

        std::set<int> id_values;
        for (int i = 0; i < ma_coords.size(); i++) 
        {
            id_values.insert(segment_ids[i]);
            outfile << ma_coords[i][0] << "," << ma_coords[i][1] << "," << ma_coords[i][2] << "," << segment_ids[i] << std::endl;
        }
        outfile.close();
        output("segment_ids").set(segment_ids);
        std::cout << "total classes:" << id_values.size() << std::endl;
        for (set<int>::iterator it = id_values.begin(); it != id_values.end(); it++)
        {
            std::cout << *it << " occurs " << std::endl;
        }
    }

    void ShowClusterMAT::process()
    {
        //----------input-----------------//
        int classifiaction = param<int>("classification");
        auto segment_ids = input("segment_ids").get<vec1i>();
        auto MAT_points = input("ma_coords").get<PointCollection>();
        // ------------output------------//
        PointCollection selected_mat;
        // ----------process -----------//
        std::set<int> s;
        for (int i = 0; i < segment_ids.size(); i++) 
        {
            s.insert(segment_ids[i]);
        }       
        
        for (int i = 0; i < segment_ids.size(); i++) 
        {
            if (segment_ids[i] == classifiaction) 
            {
                selected_mat.push_back({ MAT_points[i][0],MAT_points[i][1],MAT_points[i][2] });
            }
        }

        output("mat").set(selected_mat);
    }
    void GetClusterSheets::process() 
    {
        std::cout << "Trying to get all the sheets" << std::endl;
        //-------------input---------------//
        auto segment_ids = input("segment_ids").get<vec1i>();
        auto MAT_points = input("ma_coords").get<PointCollection>();
        auto bisectors = input("ma_bisector").get<vec3f>();
        auto radii = input("ma_radii").get<vec1f>();
        //---------output-----------//
        std::vector<PointCollection> vec_sheets;
        std::vector<float> vec_bisectors;
        std::vector<vec1f> vec_radii;
        std::vector<vec1i> vec_index;
        //-------process------//
        std::set<int> id_values;
        
        //std::cout << "size of ALL MAT points:" << MAT_points.size() << std::endl;
        //std::cout << "size of segment_id:" << segment_ids.size() << std::endl;
        // segment_id is the cluster id of all MAT, size = size of ALL MAT

        for (int i = 0; i < segment_ids.size(); i++)
        {
            id_values.insert(segment_ids[i]);            
        }

        int num = id_values.size(); // how many different clusters generated

        for (set<int>::iterator it = id_values.begin(); it != id_values.end(); it++)
        {
            PointCollection one_mat_sheet;
            float one_sheet_bisector =0;
            vec1f one_radii;
            vec1i one_index;

            for (int i = 0; i < segment_ids.size(); i++) 
            {       
                
                if (*it == segment_ids[i]) 
                {
                    
                    one_mat_sheet.push_back({ MAT_points[i][0],MAT_points[i][1],MAT_points[i][2] });
                    one_sheet_bisector += bisectors[i][2];
                    one_index.push_back(i);
                    one_radii.push_back(radii[i]);
                    
                }
                //one_sheet_bisector = one_sheet_bisector / count;
            }
            vec_sheets.push_back(one_mat_sheet);
            vec_index.push_back(one_index);
            vec_radii.push_back(one_radii);
            std::cout << "one sheet bisector value:" << one_sheet_bisector << std::endl;
            vec_bisectors.push_back(one_sheet_bisector);
        }
                     
        output("vec_sheets").set(vec_sheets);
        output("vec_bisectors").set(vec_bisectors);
        output("vec_radii").set(vec_radii);
        output("vec_index").set(vec_index);

        std::cout << "sheets done" << std::endl;
        std::cout << "bisector group size:" << vec_bisectors.size() << std::endl;
        std::cout << "in total:" << vec_sheets.size() << " obtained" << std::endl;
    }
    void MATSeparation::process() 
    {
        std::cout << "MAT separation starts" << std::endl;
        // ---------  input ------------//
        auto vec_sheets = input("vec_sheets").get<std::vector<PointCollection>>();
        auto vec_bisectors = input("vec_bisectors").get<std::vector<float>>();
        auto vec_radii = input("vec_radii").get<std::vector<vec1f>>();
        auto vec_index = input("vec_index").get<std::vector<vec1i>>();

        //auto points = input("points").get<PointCollection>();
        int offset = param<float>("offset");
        
        // ------------output-----------//
        PointCollection interior_MAT;
        PointCollection exterior_MAT;
        PointCollection unclassified_MAT;

        
        vec1f ex_radii;
        vec1f in_radii;
        vec1i in_index;
        vec1i ex_index;


        //-----------process-------------//
        int num = 0;
        for (auto sheet : vec_sheets) 
        {
            for (auto pt : sheet)
                num++;
        }
        std::cout << "input mat point size:" << num << std::endl;

        for (int i=1;i<vec_bisectors.size();i++) 
        {
            if (vec_bisectors[i] < 0) 
            {
                //std::cout << "ex sheet size:" << vec_sheets[i].size() << std::endl;
                for (auto pt : vec_sheets[i]) 
                {
                    
                    exterior_MAT.push_back(pt);
                    
                }
                for (auto r : vec_radii[i]) 
                {
                    ex_radii.push_back(r);
                }
                for (auto id : vec_index[i]) 
                {
                    ex_index.push_back(id);
                }
                
            }
            //if (vec_bisectors[i] >= 0) 
            else
            {
                //std::cout << "in sheet size:" << vec_sheets[i].size() << std::endl;
                for (auto pt : vec_sheets[i])
                {
                    interior_MAT.push_back(pt);
                }
                for (auto r : vec_radii[i])
                {
                    in_radii.push_back(r);
                }
                for (auto id : vec_index[i])
                {
                    in_index.push_back(id);
                }
            }
        }
                       



        output("interior_mat").set(interior_MAT);
        output("exterior_mat").set(exterior_MAT);
        output("unclassified_mat").set(unclassified_MAT);
        output("in_radii").set(in_radii);
        output("ex_radii").set(ex_radii);
        output("in_index").set(in_index);
        output("ex_index").set(ex_index);

    }

    void MATsimplification::process() {
        std::cout << "simplification is running" << std::endl;
        //-------------input----------------//
        auto point_collection = input("interior_mat").get<PointCollection>();
        auto radii = input("interior_radii").get<vec1f>();
        auto indice = input("interior_idx").get<vec1i>();
        float threshold = input("threshold").get<float>();
        //---------output----------------//
        PointCollection sim_mat;        
        vec1f  sim_radii;
        std::vector<std::vector<int>> sim_idx; 


        


        if (threshold == 0) 
        {
            sim_mat = point_collection;
            sim_radii = radii;
            for (auto idx : indice) 
            {
                std::vector<int> idx_vec;
                idx_vec.push_back(idx);
                sim_idx.push_back(idx_vec);
            }
        }
        else 
        {

            std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\Sim_MAT_out.txt";
            std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);

            std::vector<KdTree::sphere> MAT_sph;
            KdTree::sphere sp1;
            Vector3D p1 = { point_collection[0][0],point_collection[0][1],point_collection[0][2] };
            sp1.pos = p1;
            sp1.radius = radii[0];
            sp1.index.push_back(indice[0]);

            MAT_sph.push_back(sp1);


            for (int i = 0; i < point_collection.size(); i++)
            {
                Vector3D v1 = { point_collection[i][0],point_collection[i][1],point_collection[i][2] };
                bool ifsim = MATsimplification::ifSimplify(v1, radii[i], indice[i], threshold, MAT_sph);
            }
            
            std::cout << "MAT_to_kd done" << std::endl;
            
            for (int i = 0; i < MAT_sph.size(); i++)
            {
                arr3f point = { MAT_sph[i].pos.x,MAT_sph[i].pos.y,MAT_sph[i].pos.z };
                sim_mat.push_back(point);
                sim_radii.push_back(MAT_sph[i].radius);
                sim_idx.push_back(MAT_sph[i].index);
                outfile << MAT_sph[i].pos.x << "," << MAT_sph[i].pos.y << "," << MAT_sph[i].pos.z << "," << MAT_sph[i].radius << "," << MAT_sph[i].index.size() << std::endl;

            }
            outfile.close();
        }
        

                                     
        output("interior_mat").set(sim_mat);
        output("interior_radii").set(sim_radii);
        output("interior_idx").set(sim_idx);
        std::cout << "Simplified MAT size:" << sim_mat.size() << std::endl;
        
    }

    void ComputeNormalsNode::process() {
        //--------------input----------------//
        auto point_collection = input("points").get<PointCollection>();

        


        masb::ma_data madata;
        madata.m = point_collection.size();
        masb::PointList coords;
        coords.reserve(madata.m);
        for (auto& p : point_collection) {
            coords.push_back(masb::Point(p.data()));
        }
        masb::VectorList normals(madata.m);
        madata.coords = &coords;
        madata.normals = &normals;

        masb::compute_normals(params, madata);

        //std::string filepath = "c:\\users\\tengw\\documents\\git\\Reorien\\las_normal.txt";
        //std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);

        vec3f normals_vec3f;
        normals_vec3f.reserve(madata.m);
        for (auto& n : *madata.normals) {
            normals_vec3f.push_back({ n[0], n[1], n[2] });
            //outfile << n[0] << "," << n[1] << "," << n[2] << std::endl;
        }
        //outfile.close();
        output("normals").set(normals_vec3f);
    }

    
    void MultiKDtree::process() {
       /* auto point_collection = input("points").get<PointCollection>();
        masb::ma_data madata;
        madata.m = point_collection.size();
        auto number = point_collection.size();
        vec3f points1 ;
        vec3f points2 ;
        Vector3D* mp_Points = new Vector3D[number];

        KdTree* kd = NewBuildKdTree(mp_Points, number);
        Vector3D center = (*kd).centerofBoundingBox();
        for (int i = 0; i < number; i++) {
            if (mp_Points[i][0] < center.x) {
                points1.push_back({mp_Points[i][0],mp_Points[i][1],mp_Points[i][2] });
            }
            else {
                points2.push_back({ mp_Points[i][0],mp_Points[i][1],mp_Points[i][2] });
            }
        }
        int count1 = points1.size();
        int count2 = points2.size();
        Vector3D* m_Points1 = new Vector3D[count1];
        Vector3D* m_Points2 = new Vector3D[count2];
        for (int i = 0; i < count1;i++) {
            m_Points1[i].x = points1[i][0];
            m_Points1[i].y = points1[i][1];
            m_Points1[i].z = points1[i][2];
        }
        KdTree* kd1 = NewBuildKdTree(m_Points1, count1);

        for (int j = 0; j < count2; j++) {
            m_Points2[j].x = points2[j][0];
            m_Points2[j].y = points2[j][1];
            m_Points2[j].z = points2[j][2];
        }
        KdTree* kd2 = NewBuildKdTree(m_Points2, count2);
        
        output("KDTree1").set(kd1);
        output("KDTree2").set(kd2);*/
    }

    void BuildKDtree::process()
    {
        clock_t starttime, endtime;
        starttime = clock();

        auto point_collection = input("points").get<PointCollection>();
        auto radii = input("radii").get<vec1f>();
        //////////////////////////////
        //auto indice = input("indice").get<vec1i>();
        auto indice = input("indice").get<std::vector<vec1i>>();
        masb::ma_data madata;
        madata.m = point_collection.size();
        m_nPoints = point_collection.size();
        std::cout << "MAT point size:" << madata.m << std::endl;
        KdTree::sphere* mp_Points = new KdTree::sphere[m_nPoints];
        Vector3D* Points = new Vector3D[m_nPoints];


        Vector3D center;

        vec3f points;
        points.reserve(madata.m);
        for (auto& p : point_collection) {
            points.push_back({ p[0], p[1], p[2] });
        }
        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\MAT_points_out.txt";
        std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);

        std::vector<KdTree::sphere> allPointsVec;


        for (int i = 0; i < m_nPoints; i++) {




            mp_Points[i].pos.x = points[i][0];
            mp_Points[i].pos.y = points[i][1];
            mp_Points[i].pos.z = points[i][2];
            mp_Points[i].radius = radii[i];
            mp_Points[i].index = indice[i];
            //mp_Points[i].index.push_back(indice[i]);
            

            Points[i].x = points[i][0];
            Points[i].y = points[i][1];
            Points[i].z = points[i][2];




            outfile<< mp_Points[i].pos.x << "," << mp_Points[i].pos.y << "," << mp_Points[i].pos.z <<","<< mp_Points[i].radius << std::endl;
            allPointsVec.push_back(mp_Points[i]);
        }        
        outfile.close();
        
        KdTree* kd = BuildKDtree::BuildKdTree(Points, m_nPoints,20);
        center = (*kd).centerofBoundingBox();
        (*kd).m_allballs = allPointsVec;
       
        //std::cout << "center point of bounding box:" << center.x << "," << center.y << "," << center.z << std::endl;
        std::cout << "Total number of points in kd-tree"<< allPointsVec.size() << std::endl;
        std::cout<< "number of bounding box:" << (*kd).m_maxpoint.size() << std::endl;
        //std::cout << "Vector size of kdtreepoints:" << (*kd).m_boxKdTreePoint.size() << std::endl;
        std::cout << "size of level:" << (*kd).m_currentlevel.size() << std::endl;
        
        std::cout << "Level points are saving" << std::endl;

        
        int new_count = 0;
        for (int i =0; i<(*kd).m_maxpoint.size() ; i++)
        {

            //std::vector<KdTree::sphere> levelpoints = BuildKDtree::NewGetLevelPoints((*kd).m_maxpoint[i], (*kd).m_minpoint[i], &allPointsVec);
            std::vector<KdTree::sphere> levelpoints = BuildKDtree::GetLevelPoints((*kd).m_maxpoint[i], (*kd).m_minpoint[i], &allPointsVec);


            new_count +=levelpoints.size() ;
            //std::cout <<"Number of points in each level:"<< levelpoints.size() << std::endl;
            (*kd).m_levelpoints.push_back(levelpoints);
        }
        std::cout << "total level points:" << new_count << std::endl;
        std::cout << "KDTree output done" << std::endl;
        


        /*std::cout << "Total level points:" << new_count << std::endl;
        std::cout << "At the end the size of allPointsVec are:" << allPointsVec.size() << std::endl;
        std::cout << "The size of m_levelpoints:" << (*kd).m_levelpoints.size() << std::endl;
        std::cout << "The size of boxs:" << (*kd).m_maxpoint.size() << std::endl;

        std::cout << "checking points in each box" << std::endl;
        int flag = 0;
        for (int i = 0; i < (*kd).m_maxpoint.size(); i++) {
            for (auto point : (*kd).m_levelpoints[(*kd).m_maxpoint.size() - i-1]) {
                if(point.pos.x<= (*kd).m_maxpoint[i].x && point.pos.x >= (*kd).m_minpoint[i].x)
                    if(point.pos.y <= (*kd).m_maxpoint[i].y && point.pos.y >= (*kd).m_minpoint[i].y)
                        if (point.pos.z <= (*kd).m_maxpoint[i].z && point.pos.z >= (*kd).m_minpoint[i].z)
                        {
                            flag++;
                        }
            
            }
        }
        std::cout << "Correct points:" << flag<<std::endl;
        std::cout << "Wrong points:" << new_count -flag << std::endl;*/
        
 
        output("KDTree").set(kd);
        endtime = clock();
        std::cout << "KD-Tree running time:" << endtime - starttime << std::endl;
    }
    
    void AMPGPUQueryTest::process() {

        clock_t starttime, endtime;
        starttime = clock();

        //----------------input------------------------//
        std::cout << "GPU query start" << std::endl;
        auto kd = input("KDTree").get<KdTree*>();
        auto point_collection = input("MATpoints").get<PointCollection>();
        auto radii = input("radii").get<vec1f>();
        //auto indice = input("indice").get<vec1i>();
        //auto indice = input("indice").get<std::vector<vec1i>>();

        auto interval = input("interval").get<float>();

        // -------------output----------------------//
        PointCollection visible_mat;
        vec1f visible_radii;
        vec1i visible_indice;

        // -------------------------------------//
        masb::ma_data madata;
        madata.m = point_collection.size();
        int number = point_collection.size();
        std::cout << "MAT point size:" << madata.m << std::endl;


        Vector3D* mp_Points = new Vector3D[number];




        vec3f points;
        points.reserve(madata.m);
        for (auto& p : point_collection) {
            points.push_back({ p[0], p[1], p[2] });
        }
        for (int i = 0; i < number; i++) {
            mp_Points[i].x = points[i][0];
            mp_Points[i].y = points[i][1];
            mp_Points[i].z = points[i][2];
        }

        std::vector<KdTree::sphere> visble_sph;

        Vector3D v1 = input("Vector1").get<Vector3D>();
        Vector3DNew v1_new(v1);

        std::vector<Vector3D> v2_list = AMPGPUQueryTest::SpherePoints(v1, 500, interval);
        int linesize = v2_list.size();

        std::vector<Vector3DNew> v2_new;
        for (auto a : v2_list) 
        {
            Vector3DNew new_a(a);
            v2_new.push_back(new_a);
        }
        std::vector<int> result;
                
        //run the query function here check if ray intersects with bounding boxes.
        std::cout << "Ray-boxes check starts" << std::endl;

        AMPGPUQueryTest::GPUQuery(v2_new, v1_new, linesize, kd, result);
        
        std::cout << "Ray-boxes check done" << std::endl;

        int GPU_box_check_count = 0;

        std::vector<Vector3D> v2_intersected;
        // processing bar                      //

        //float progress = 0.0;
        //while (progress < 1.0) {
        //    int barWidth = 70;

        //    std::cout << "[";
        //    int pos = barWidth * progress;
        //    for (int i = 0; i < barWidth; ++i) {
        //        if (i < pos) std::cout << "=";
        //        else if (i == pos) std::cout << ">";
        //        else std::cout << " ";
        //    }
        //    std::cout << "] " << int(progress * 100.0) << " %\r";
        //    std::cout.flush();

        //    progress += 0.16; // for demonstration only
        //}
        //std::cout << std::endl;

        for (int i = 0; i < result.size();i++) {
            if (result[i] == 1) 
            { 

                //v2_intersected.push_back(v2_list[i]);
                
                GPU_box_check_count++;
                
                auto sph1 = AMPGPUQueryTest::GetOneLineResult(v1,v2_list[i],kd);
                visble_sph.push_back(sph1);
            }
            
        }
        std::cout << "Number of intersected directions:" << GPU_box_check_count << std::endl;
        //---------------Mutiple Threads-----------------------//        
       /* auto CutVecList = MutiThreadsOneQuery::CutVecList(v2_intersected, 4);
        
        std::thread t[4];


        int Threads_Count = 0;
        for (std::vector<Vector3D> vec_cut : CutVecList) {
            t[Threads_Count] = std::thread(AMPGPUQueryTest::GetQueryResult, vec_cut, v1, kd, std::ref(visble_sph));
            t[Threads_Count].join();
            Threads_Count++;
        }*/
                              
        



        //std::cout << "Visble MAT size:" << visble_sph.size() << std::endl;

        

        /*for (auto item : visble_sph) {
            visible_mat.push_back({ item.pos.x,item.pos.y,item.pos.z });
            visible_radii.push_back(item.radius);
            visible_indice.push_back(item.index);
        }*/
        ////////////////////
        for (auto item : visble_sph) {
            visible_mat.push_back({ item.pos.x,item.pos.y,item.pos.z });
            visible_radii.push_back(item.radius);
            for(auto idx: item.index)
                visible_indice.push_back(idx);
        }


        endtime = clock();

        output("MAT_points").set(visible_mat);
        output("radii").set(visible_radii);
        output("indice").set(visible_indice);
        std::cout << "GPUNode running time:" << endtime - starttime << std::endl;

    }

    void MutiThreadsOneQuery::process() {
        //----------------input------------------------//
        std::cout << "Mutiple Threads Query start" << std::endl;
        auto kd = input("KDTree").get<KdTree*>();
        auto point_collection = input("MATpoints").get<PointCollection>();
        auto radii = input("radii").get<vec1f>();
        // -------------output----------------------//
        PointCollection visible_mat;
        vec1f visible_radii;

        // -------------------------------------//
        masb::ma_data madata;
        madata.m = point_collection.size();
        int number = point_collection.size();
        std::cout << "MAT point size:" << madata.m << std::endl;


        Vector3D* mp_Points = new Vector3D[number];




        vec3f points;
        points.reserve(madata.m);
        for (auto& p : point_collection) {
            points.push_back({ p[0], p[1], p[2] });
        }
        for (int i = 0; i < number; i++) {
            mp_Points[i].x = points[i][0];
            mp_Points[i].y = points[i][1];
            mp_Points[i].z = points[i][2];
        }

        std::vector<sphere> visble_sph;

        Vector3D v1 = input("Vector1").get<Vector3D>();

        auto v2_list = OneQuery::SpherePoints(v1, 500);

        clock_t starttime, endtime;
        starttime = clock();
        auto CutVecList = MutiThreadsOneQuery::CutVecList(v2_list, 4);


        std::thread t[4];


        int Threads_Count = 0;
        for (std::vector<Vector3D> vec_cut : CutVecList) {
            t[Threads_Count]=std::thread(MutiThreadsOneQuery::GetQueryResult, vec_cut, v1, kd, std::ref(visble_sph));
            t[Threads_Count].join();
            Threads_Count++;
        }


        


        ////////////////////////////////
        /*int vetsize = v2_list.size();
        std::vector<Vector3D> threadvec1;
        std::vector<Vector3D> threadvec2;
        std::vector<Vector3D> threadvec3;
        std::vector<Vector3D> threadvec4;


        std::for_each(begin(v2_list), begin(v2_list) + 0.25*vetsize, [&threadvec1](Vector3D x) {
            threadvec1.push_back(x);
        });
        std::for_each(begin(v2_list) + 0.25*vetsize, begin(v2_list) + 0.5*vetsize, [&threadvec2](Vector3D y) {
            threadvec2.push_back(y);
        });
        std::for_each(begin(v2_list) + 0.5*vetsize, begin(v2_list) + 0.75*vetsize, [&threadvec3](Vector3D z) {
            threadvec3.push_back(z);
        });
        std::for_each(begin(v2_list) + 0.75*vetsize, end(v2_list) , [&threadvec4](Vector3D v) {
            threadvec4.push_back(v);
            
        });


        std::cout << "Vector divide done" << std::endl;


        
        std::thread t1(MutiThreadsOneQuery::GetQueryResult, threadvec1, v1,kd,std::ref(visble_sph) );
        t1.join();
        std::thread t2(MutiThreadsOneQuery::GetQueryResult, threadvec2, v1, kd, std::ref(visble_sph));
        t2.join();
        std::thread t3(MutiThreadsOneQuery::GetQueryResult, threadvec3, v1, kd, std::ref(visble_sph));
        t3.join();
        std::thread t4(MutiThreadsOneQuery::GetQueryResult, threadvec4, v1, kd, std::ref(visble_sph));
        t4.join();
        */

        //////////////////////////////
        //2739

        std::cout << "Multi Threads query done" << std::endl;
        std::cout << "Visble MAT size:" << visble_sph.size() << std::endl;


        for (auto item : visble_sph) {
            visible_mat.push_back({ item.pos.x,item.pos.y,item.pos.z });
            visible_radii.push_back(item.r);
        }
        endtime = clock();
        
        output("MAT_points").set(visible_mat);
        output("radii").set(visible_radii);
        std::cout << "Multiple threads running time:" << endtime - starttime << std::endl;

    }
    void GetRadialRayResults::process() 
    {
        std::cout << "get radial rays result starts" << std::endl;
        clock_t starttime, endtime;
        starttime = clock();
        // ---------------input -------------------//
        auto kd = input("KDTree").get<KdTree*>();
        auto point_collection = input("MATpoints").get<PointCollection>();
        auto radii = input("radii").get<vec1f>();
        auto viewpoint = input("viewpoint").get<Vector3D>();
        auto radial_vectors = input("radial_rays").get<std::vector<Vector3D>>();
        std::cout << "number of rays:" << radial_vectors.size()<< std::endl;
        // ----------------- output -----------------//
        std::vector<KdTree::sphere> visble_sph;

        PointCollection visible_mat;
        vec1f visible_radii;
        vec1i visible_indice;

        //----------process------------//
        std::vector<bool> ifinter;

        for (int j = 0; j < radial_vectors.size(); j++)
        {
            Vector3D hit;
            bool inter = GetRaysResult::CheckLineBox((*kd).m_min, (*kd).m_max, viewpoint, radial_vectors[j], hit);
            ifinter.push_back(inter);
        }
        int count_intersect = 0;
        for (int j = 0; j < radial_vectors.size(); j++)
        {
            if (ifinter[j] == 1)
            {
                count_intersect++;
                auto sph1 = AMPGPUQueryTest::GetOneLineResult(viewpoint, radial_vectors[j], kd);
                visble_sph.push_back(sph1);
            }
        }
        std::cout << "!!!!! number of intersection:" << count_intersect << std::endl;

        for (auto item : visble_sph) {
            visible_mat.push_back({ item.pos.x,item.pos.y,item.pos.z });
            visible_radii.push_back(item.radius);
            for (auto idx : item.index)
                visible_indice.push_back(idx);
        }

        std::cout <<"vis_id size"<< visible_indice.size() << std::endl;

        endtime = clock();

        output("MAT_points").set(visible_mat);
        output("radii").set(visible_radii);
        output("indice").set(visible_indice);
        std::cout << "Get rays result running time:" << endtime - starttime << std::endl;
        
    }
    void GetRaysResult::process()
    {
        std::cout << "get rays result starts" << std::endl;
        clock_t starttime, endtime;
        starttime = clock();
        // --------------input---------------------//
        auto kd = input("KDTree").get<KdTree*>();
        auto point_collection = input("MATpoints").get<PointCollection>();
        auto radii = input("radii").get<vec1f>();
        auto headvectors = input("Headvectors").get<std::vector<Vector3D>>();
        auto endvectors = input("Endvectors").get<std::vector<Vector3D>>();
        std::cout << "the number of rays: " << endvectors.size() << std::endl;
        //-------output-------------------//
        std::vector<KdTree::sphere> visble_sph;

        PointCollection visible_mat;
        vec1f visible_radii;
        vec1i visible_indice;

        //----------process------------//
        // check ray intersect with box first here AMP query //
        std::vector<bool> ifinter;

        for (int j = 0; j < headvectors.size(); j++) 
        {
            Vector3D hit;
            bool inter = CheckLineBox((*kd).m_min,(*kd).m_max,headvectors[j], endvectors[j],hit);
            ifinter.push_back(inter);
        }

        int count_intersect = 0;
        for (int j = 0; j < headvectors.size(); j++) 
        {
            if (ifinter[j] == 1) 
            {
                count_intersect++;
                auto sph1 = AMPGPUQueryTest::GetOneLineResult(headvectors[j], endvectors[j], kd);
                visble_sph.push_back(sph1);
            }
        }        
        std::cout << "!!!!! number of intersection:" << count_intersect << std::endl;
        for (auto item : visble_sph) {
            visible_mat.push_back({ item.pos.x,item.pos.y,item.pos.z });
            visible_radii.push_back(item.radius);
            for (auto idx : item.index)
                visible_indice.push_back(idx);
        }



        endtime = clock();

        output("MAT_points").set(visible_mat);
        output("radii").set(visible_radii);
        output("indice").set(visible_indice);
        std::cout << "Get rays result running time:" << endtime - starttime << std::endl;
    }

    void OneQuery::process() {
        std::cout << "OneQuery start" << std::endl;

        //----------------input------------------------//

        auto kd = input("KDTree").get<KdTree*>();
        auto point_collection = input("MATpoints").get<PointCollection>();
        auto radii = input("radii").get<vec1f>();


        // -------------output----------------------//
        PointCollection visible_mat;
        vec1f visible_radii;

        // -------------------------------------//

        masb::ma_data madata;
        madata.m = point_collection.size();
        int number = point_collection.size();
        std::cout << "MAT point size:" << madata.m << std::endl;


        Vector3D* mp_Points = new Vector3D[number];




        vec3f points;
        points.reserve(madata.m);
        for (auto& p : point_collection) {
            points.push_back({ p[0], p[1], p[2] });
        }
        for (int i = 0; i < number; i++) {
            mp_Points[i].x = points[i][0];
            mp_Points[i].y = points[i][1];
            mp_Points[i].z = points[i][2];
        }

        Vector3D v1 = input("Vector1").get<Vector3D>();

        auto v2_list = OneQuery::SpherePoints(v1, 500);
        //Vector3D v2 = input("Vector2").get<Vector3D>();

        

        for (Vector3D v2 : v2_list) {
            std::vector<Vector3D> pointlist;
            std::vector<float> radiilist;
            Vector3D hit;
            int count = 0;
            for (int i = 0; i < (*kd).m_maxpoint.size(); i++) {
                bool a = OneQuery::CheckLineBox((*kd).m_minpoint[i], (*kd).m_maxpoint[i], v1, v2, hit);
                //std::cout << "a��" <<a<< std::endl;
                if (a == 1) {
                    //std::cout << "intersect bounding box level:" << i << std::endl;
                    //std::cout << "points inside" << std::endl;
                    for (auto pt : (*kd).m_levelpoints[(*kd).m_maxpoint.size() - i - 1]) {
                        count++;
                        //std::cout << pt_index << std::endl;

                        //std::cout << mp_Points[pt_index].x << "," << mp_Points[pt_index].y << "," << mp_Points[pt_index].z << std::endl;
                        float dis = VisibiltyQurey::DistanceOfPointToLine(v1, v2, pt.pos);
                        //std::cout << "distance to line:" << dis << std::endl;
                        
                        if (dis <= pt.radius) {
                            //std::cout << "Radius:" << pt.radius<<std::endl;
                            pointlist.push_back(pt.pos);
                            radiilist.push_back(pt.radius);
                            //visible_mat.push_back({ mp_Points[pt_index].x,mp_Points[pt_index].y,mp_Points[pt_index].z });
                            //visible_radii.push_back(radii[pt_index]);
                        }
                    }
                    //std::cout << "total points inside:" << count << std::endl;
                }
            }
            

            if (pointlist.size() > 0) {
                std::cout << "----------------this direction has :" << pointlist.size()<<"   intersected"<< std::endl;
                float minDis = OneQuery::PointToPointDis(v1, pointlist[0]);
                //float minDis = OneQuery::PointToPointDis(v1, pointlist[0]) - radiilist[0];
                
                int flag = 0;
                for (int i = 0; i < pointlist.size(); i++)
                {
                    //float temp = OneQuery::PointToPointDis(v1, pointlist[i]) - radiilist[i];
                    float temp = OneQuery::PointToPointDis(v1, pointlist[i]);
                    if (temp < minDis) {
                        minDis = temp;
                        flag = i;
                    }
                }
                visible_mat.push_back({ pointlist[flag].x,pointlist[flag].y,pointlist[flag].z });
                visible_radii.push_back(radiilist[flag]);
            }
        }

        output("MAT_points").set(visible_mat);
        output("radii").set(visible_radii);

    }

    void KDTreeNearestQurey::process() {
        //auto viewpoint = input("viewpoint").get<PointCollection>();
        std::cout << "Nearest points qurey starting" << std::endl;
        int number = int(input("Number").get<float>());

        PointCollection resultPoints;
        resultPoints.reserve(number);

        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\Query_points_output.txt";
        std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);

        if (input("KDTree").has_data())
        {
            auto kdtree = input("KDTree").get<KdTree*>();
            std::cout << "kdtree output successfully" << std::endl;
            auto MATpoints = input("MATpoints").get<PointCollection>();


            //Vector3D v1 = { -99.0594,-90.828,6.59866 };

            //Vector3D v2 = { 0,0,0 };
            //(*kdtree).queryLineIntersection(v1, v2, 100, 1, 1);
            //for (int i = 0; i < number; i++)
            //{
            //    int index = (*kdtree).getNeighbourPositionIndex(i);
            //    std::cout << "index   " << index << std::endl;
            //    std::cout << "x   " << MATpoints[index][0] << std::endl;
            //    std::cout << "y   " << MATpoints[index][1] << std::endl;
            //    std::cout << "z   " << MATpoints[index][2] << std::endl;
            //    resultPoints.push_back({ MATpoints[index][0], MATpoints[index][1] ,MATpoints[index][2] });
            //}
            //output("Points").set(resultPoints);

            //Vector3D v1 = { -99.0594,-90.828,6.59866 };
            Vector3D v1 = input("ViewPoint").get<Vector3D>();
            (*kdtree).setNOfNeighbours(number);
            (*kdtree).queryPosition(v1);
            for (int i = 0; i < number; i++) {
                int index = (*kdtree).getNeighbourPositionIndex(i);                
                outfile << MATpoints[index][0] << "," << MATpoints[index][1] << "," << MATpoints[index][2] << std::endl;
                resultPoints.push_back({ MATpoints[index][0], MATpoints[index][1] ,MATpoints[index][2] });
            }
            output("Points").set(resultPoints);

        }
        else std::cout << "KDTree is empty " << std::endl;
        outfile.close();  
        std::cout << "KDtree nearest query results done." << std::endl;
    }
    void KDTreeLineQurey::process() {

        std::cout << "Line qurey starting" << std::endl;
        //---------------Input Data--------------
        int number = int(input("Number").get<float>());
        float distance = input("Distance").get<float>();

        Vector3D v1 = input("Vector1").get<Vector3D>();
        Vector3D v2 = input("Vector2").get<Vector3D>();
        auto kdtree = input("KDTree").get<KdTree*>();
        auto MATpoints = input("MATpoints").get<PointCollection>();
        //---------------------Output Data-----------
        PointCollection resultPoints;
        resultPoints.reserve(number);
        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\LineQuery_points_output.txt";
        std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);
        (*kdtree).queryLineIntersection(v1, v2, distance, 1, 1);
        for (int i = 0; i < number; i++)
        {
            int index = (*kdtree).getNeighbourPositionIndex(i);
            std::cout << "index   " << index << std::endl;
            std::cout << "x   " << MATpoints[index][0] << std::endl;
            std::cout << "y   " << MATpoints[index][1] << std::endl;
            std::cout << "z   " << MATpoints[index][2] << std::endl;
            resultPoints.push_back({ MATpoints[index][0], MATpoints[index][1] ,MATpoints[index][2] });
            outfile << MATpoints[index][0] << "," << MATpoints[index][1] << "," << MATpoints[index][2] << std::endl;
        }
        output("Points").set(resultPoints);
        outfile.close();
    }
    void VisiblePCbyRTree::process()
    {
        std::cout << "Visible PC query by RTree starts " << std::endl;
        clock_t starttime, endtime, endtimeRTree;
        starttime = clock();
        ////---------------------input----------------------------//
        Vector3D viewpoint = input("viewPoint").get<Vector3D>();
        auto interior_MAT = input("interior_MAT").get<PointCollection>();
        auto radii = input("interior_radii").get<vec1f>();
        auto pc = input("original_pc").get<PointCollection>();
        auto in_idx = input("in_index").get<vec1i>();
        std::cout << "size of interior MAT" << interior_MAT.size() << std::endl;
        std::cout << "size of in_idx" << in_idx.size() << std::endl;

        std::cout << "input pc size:" << pc.size() << std::endl;
        ////-------------output -----------------------//
        PointCollection visible_pc;

        //std::vector<int> vec_result_id;
        std::set<int>  vec_result_id2;

        ////-----------------process-------------------//        
        namespace bg = boost::geometry;
        namespace bgi = boost::geometry::index;
        typedef bg::model::point<float, 3, bg::cs::cartesian> point;
        typedef bg::model::box<point> box;
        typedef bg::model::segment<point> seg;
        typedef std::pair<box, unsigned> value;

        // create the rtree using default constructor
        bgi::rtree< value, bgi::quadratic<16> > rtree;

        for (int i = 0; i < interior_MAT.size(); i++)
        {
            box b(point(interior_MAT[i][0] - radii[i], interior_MAT[i][1] - radii[i], interior_MAT[i][2] - radii[i]), point(interior_MAT[i][0] + radii[i], interior_MAT[i][1] + radii[i], interior_MAT[i][2] + radii[i]));
            rtree.insert(std::make_pair(b, i));
        }
        std::cout << "RTree node inserted done!" << std::endl;
        endtimeRTree = clock();
        std::cout << "Running time of RTree construction: " << endtimeRTree - starttime << std::endl;
        for (int i = 0; i < interior_MAT.size(); i++) 
        {
            std::vector<value> result_s;
            seg query_seg(point(viewpoint[0], viewpoint[1], viewpoint[2]), point(interior_MAT[i][0], interior_MAT[i][1], interior_MAT[i][2]));
            rtree.query(bgi::intersects(query_seg), std::back_inserter(result_s));

            if (result_s.size() == 0)
                std::cout << "no box intersection!" << "ID" << i << std::endl;      
            if (result_s.size() == 1)
                std::cout << "only intersect with target pt itself!" << i << std::endl;

            //int result_id = result_s[0].second;
            int result_id = i;
            float min_dis = PointToPointDis(viewpoint, Vector3D(interior_MAT[i][0], interior_MAT[i][1], interior_MAT[i][2])) - radii[i];
           
            for (auto item : result_s) 
            {                                
                int id = item.second;
                if (DistancePointToSegment(viewpoint, Vector3D(interior_MAT[i][0], interior_MAT[i][1], interior_MAT[i][2]), Vector3D(interior_MAT[id][0], interior_MAT[id][1], interior_MAT[id][2]))< radii[id])
                {
                    float current_dis = PointToPointDis(viewpoint, Vector3D(interior_MAT[id][0], interior_MAT[id][1], interior_MAT[id][2])) -radii[id];
                    if (current_dis < min_dis)
                    {
                        result_id = id;
                        min_dis = current_dis;
                    }
                }
            }
            //vec_result_id.push_back(result_id);     
            vec_result_id2.insert(in_idx[result_id]);
        }

        //std::cout << "Size of vector result id: " << vec_result_id.size() << std::endl;

        for (auto id : vec_result_id2) 
        {
            if (id > pc.size()) 
            {
                id = id - pc.size();
                //std::cout << "index out of range!" << std::endl;
            }
            visible_pc.push_back(pc[id]);
        }

        //std::cout << "Size of set: " << vec_result_id2.size() << std::endl;
        std::cout << "Size of visible pc: " << visible_pc.size() << std::endl;
        endtime = clock();
        std::cout << "running time:" << endtime - starttime <<"ms"<< std::endl;
        output("visible_pc").set(visible_pc);
        

    }
    void VisiblePC::process() 
    {
        std::cout << "Visible PC query starts" << std::endl;
        std::cout << "Using KD-tree:" << if_checkbox<<std::endl ;

        clock_t starttime, endtime;
        starttime = clock();
        //-------------input---------------//
        Vector3D viewpoint = input("viewPoint").get<Vector3D>();
        auto kd = input("KDTree").get<KdTree*>();
        auto interior_MAT = input("interior_MAT").get<PointCollection>();
        auto radii = input("interior_radii").get<vec1f>();


        auto pc = input("original_pc").get<PointCollection>();
        std::cout << "input pc size:" << pc.size() << std::endl;
        //------------output----------------//       
        PointCollection visible_pc;
        // ---------process ------------------//


        /*int vetsize = pc.size();
        std::vector<arr3f> threadvec1;
        std::vector<arr3f> threadvec2;
        std::vector<arr3f> threadvec3;
        std::vector<arr3f> threadvec4;


        std::for_each(begin(pc), begin(pc) + 0.25*vetsize, [&threadvec1](arr3f x) {
            threadvec1.push_back(x);
        });
        std::for_each(begin(pc) + 0.25*vetsize, begin(pc) + 0.5*vetsize, [&threadvec2](arr3f y) {
            threadvec2.push_back(y);
        });
        std::for_each(begin(pc) + 0.5*vetsize, begin(pc) + 0.75*vetsize, [&threadvec3](arr3f z) {
            threadvec3.push_back(z);
        });
        std::for_each(begin(pc) + 0.75*vetsize, end(pc) , [&threadvec4](arr3f v) {
            threadvec4.push_back(v);
        });

        std::thread th[4];
        
        th[0] = std::thread(VisiblePC::GetVisblePT, threadvec1, interior_MAT, radii, viewpoint, std::ref(visible_pc));
        th[0].join();
        th[1] = std::thread(VisiblePC::GetVisblePT, threadvec2, interior_MAT, radii, viewpoint, std::ref(visible_pc));
        th[1].join();
        th[2] = std::thread(VisiblePC::GetVisblePT, threadvec3, interior_MAT, radii, viewpoint, std::ref(visible_pc));
        th[2].join();
        th[3] = std::thread(VisiblePC::GetVisblePT, threadvec4, interior_MAT, radii, viewpoint, std::ref(visible_pc));
        th[3].join();
            */
            
        if (if_checkbox == true) 
        {
            GetVisblePT(pc, kd, viewpoint, visible_pc);
        }
        else 
        {
            GetVisblePTWithoutKD(pc, viewpoint, interior_MAT,radii,visible_pc);
        }       
        endtime = clock();

        //outfile.close();
        //----------------set result ----------------//
        std::cout << "visible pc done" << std::endl;
        std::cout << "visible points size:" << visible_pc.size() << std::endl;
        std::cout << "running time:" << endtime - starttime << std::endl;
        
        output("visible_pc").set(visible_pc);

    }
    void ParallelVector::process()
    {
        std::cout << "sight vectors starts" << std::endl;
        //-----------input ---------------------//        
        Vector3D v1 = input("vector1").get<Vector3D>();
        Vector3D v2 = input("vector2").get<Vector3D>();

        float density = param<float>("Density");
        float radius = param<float>("Radius");

        
        //------------output-----------------//
        std::vector<Vector3D> Headvectors;
        std::vector<Vector3D> Endvectors;

        //-------------process-----------------//
        Vector3D normal = (v2 - v1);
        Vector3D perpen_normal(0, normal[2], -normal[1]);
        auto u = perpen_normal.normalize();
        std::cout << "Normalization:" << perpen_normal[0] << "," << perpen_normal[1] << "," << perpen_normal[2] << std::endl;
        Vector3D v = crossProduct(perpen_normal, normal);
        auto length = v.normalize();
        //std::cout << "v:" << v[0] << "," << v[1] << "," << v[2] << std::endl;
        // change radius and N to dynamic later
        //float radius = 200;
        //float density = 100;

        float delta = radius / density;
        float epsilon = delta * 0.5f;
        //std::cout << "Test" << std::endl;

        for (float y = -radius; y < radius + epsilon; y += delta)
        {
            for (float x = -radius; x < radius + epsilon; x += delta)
            {
                Headvectors.push_back(v1 + x * perpen_normal + y * v); // v1 is the point on the plane
                Endvectors.push_back(v2 + x * perpen_normal + y * v);
            }
        }
        //std::cout << "Test done" << std::endl;
        
        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\Vectors_head_end.txt";
        std::ofstream outfileVec(filepath, std::fstream::out | std::fstream::trunc);

        for (int i = 0; i < Headvectors.size(); i++)
        {
            outfileVec << Headvectors[i][0] << "," << Headvectors[i][1] << "," << Headvectors[i][2] << ";" << Endvectors[i][0] << "," << Endvectors[i][1] << "," << Endvectors[i][2] << std::endl;
        }
        outfileVec.close();

        output("Headvectors").set(Headvectors);
        output("Endvectors").set(Endvectors);
        std::cout << "number of parallel rays:" << Headvectors.size() << std::endl;
        std::cout << "sight vectors done" << std::endl;
    }
    void VisiblePart::process()     
    {
        std::cout << "Trying to get the visible part" << std::endl;
        //-------------input ---------------- //
        Vector3D viewpoint = input("viewpoint").get<Vector3D>();
        auto targetPC = input("targetPC").get<PointCollection>();
        auto kd = input("KDTree").get<KdTree*>();
        auto interior_mat = input("MATpoints").get<PointCollection>();
        auto radii = input("radii").get<vec1f>();
        // ------------ output ---------------//
        PointCollection vis_PC;
        // ------------process ---------------//
        for (auto pc : targetPC) 
        {
            Vector3D v2(pc[0], pc[1], pc[2]);
            Vector3D hit;

            //------------Chech each point visibility ---------//
            bool inter = GetRaysResult::CheckLineBox((*kd).m_min, (*kd).m_max, viewpoint, v2, hit);
            if (inter==0) 
            {
                vis_PC.push_back({v2[0],v2[1],v2[2]});
                continue;
            }
            else 
            {
                bool ifblock = 0;
                for (auto ball : (*kd).m_allballs)                {
                    
                    float dis = VisibiltyQurey::DistanceOfPointToLine(viewpoint, v2, ball.pos);
                    if (dis <= ball.radius) {
                        ifblock = 1;
                        break;
                    }                    
                }
                if (ifblock == 1) 
                {
                    continue;
                }
                vis_PC.push_back({ v2[0],v2[1],v2[2]});
            }
        }




        output("visible_parts").set(vis_PC);
        std::cout << "Visible part is done" << std::endl;
    }

    void VisibiltyQurey::process() {
        std::cout << "Visibility Qurey starts" << std::endl;
        //----------------input----------------//
        Vector3D viewpoint = input("ViewPoint").get<Vector3D>();
        auto kdtree = input("KDTree").get<KdTree*>();
        auto interior_mat = input("interior_MAT").get<PointCollection>();
        auto interior_radii = input("interior_radii").get<vec1f>();
        auto original_pc = input("original_pc").get<PointCollection>();
        //-----------output-------------------//
        PointCollection visible_mat;

        vec1f visible_radii;
        vec1i visible_index;
        // ------------process -------------//

        //std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\Visible_MAT_out.txt";
        //std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);
        

        for (int i = 0; i < interior_mat.size(); i++) {
            bool visibility = true;
            Vector3D target = { interior_mat[i][0], interior_mat[i][1],interior_mat[i][2] };
            (*kdtree).queryLineIntersection(viewpoint, target, 10, 1, 1);
            int number = (*kdtree).getNOfFoundNeighbours();
            for (int j = 0; j < number; j++)
            {
                int index = (*kdtree).getNeighbourPositionIndex(j);

                Vector3D s = { interior_mat[index][0],interior_mat[index][1],interior_mat[index][2] };
                float dis = VisibiltyQurey::DistanceOfPointToLine(viewpoint, target, s);                
                if (dis <= interior_radii[index])
                {
                    visibility = false;
                    break;
                }
            }
            if (visibility == false)
            {
                continue;
            }

            visible_mat.push_back({ interior_mat[i][0], interior_mat[i][1],interior_mat[i][2] });
            visible_radii.push_back(interior_radii[i]);
           
            //outfile << interior_mat[i][0] << "," << interior_mat[i][1] << "," << interior_mat[i][2] << "," << interior_radii[i] << std::endl;
        }
        //outfile.close();
        output("Visible_MAT").set(visible_mat);
        output("Radii_of_MAT").set(visible_radii);
        //output("indices").set(visible_index);
        std::cout << "Visiblity Qurey Done" << std::endl;
    }
 
    void WritePC2File::process() 
    {
        auto pc = input("points").get<PointCollection>();
        auto filepath = param<std::string>("filepath");
        //std::string filepath = "c:\\users\\tengw\\documents\\git\\Reorien\\PC_points.txt";
        std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);
        for (auto point : pc)
        {
            outfile << point[0] << "," << point[1] << "," << point[2] << std::endl;
        }
        outfile.close();
    }
    void ReadNormal::process() 
    {
        //----------input-------//
        auto filepath = param<std::string>("filepath");
        //----------output-----------//
        vec3f normal;
        // -----------proccess--------//
        std::ifstream infile(filepath);
        std::string line;
        if (infile)
        {
            while (getline(infile, line))
            {
                std::vector<std::string> result;
                ReadNormal::split(line, ',', result);
                float x = ReadNormal::str2num(result[0]);
                float y = ReadNormal::str2num(result[1]);
                float z = ReadNormal::str2num(result[2]);
                normal.push_back({ x,y,z });
            }
        }
        infile.close();
        std::cout <<"read normal size:" <<normal.size() << std::endl;
        output("normal").set(normal);

    }

    void TreeRemover::process() 
    {
        std::cout << "Removing trees..." << std::endl;
        //--------------input-------------//
        auto pc = input("points").get<PointCollection>(); 
        auto classification = input("classification").get<std::vector<int>>();
        int id = param<int>("number_value");

        //-------------output-------------//
        PointCollection outpc;
        //--------process-----------------//

        std::set<int> s;
        for (int i=0;i<pc.size();i++) 
        {
            s.insert(classification[i]);
            if (classification[i] != id) outpc.push_back(pc[i]);        
        } 
        for (set<int>::iterator it = s.begin(); it != s.end(); it++)
        {
            std::cout << *it << " occurs " << std::endl;
        }        

        output("points").set(outpc);
        std::cout << "Tree removing done." << std::endl;
    }
}