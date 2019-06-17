#include "masb_nodes.hpp"
#include "iostream"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <time.h>
#include<amp.h>



namespace geoflow::nodes::mat {
    

    void ComputeMedialAxisNode::process() {
        auto point_collection = input("points").get<PointCollection>();
        auto normals_vec3f = input("normals").get<vec3f>();

        /*float min_z= point_collection[0][2];
        for (int i = 0; i < point_collection.size(); i++) 
        {
            float temp = point_collection[i][2];
            if (temp < min_z) min_z = temp;
        }*/


        masb::ma_data madata;
        madata.m = point_collection.size();

        masb::PointList coords;
        coords.reserve(madata.m);
        for (auto& p : point_collection) {
            coords.push_back(masb::Point(p.data()));
        }
        
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

        masb::compute_masb_points(params, madata);

        vec1i ma_qidx;
        ma_qidx.reserve(madata.m * 2);
        for (size_t i = 0; i < madata.m * 2; ++i) {
            ma_qidx.push_back(madata.ma_qidx[i]);
        }

        PointCollection ma_coords;
        ma_coords.reserve(madata.m * 2);
        for (auto& c : *madata.ma_coords) {
            ma_coords.push_back({ c[0], c[1], c[2] });
        }



        //----- radii -----------------//
        vec1f ma_radii;
        ma_radii.reserve(madata.m * 2);
        for (size_t i = 0; i < madata.m * 2; ++i) {
            double r = Vrui::Geometry::dist(coords[i%madata.m], ma_coords_[i]);            
            ma_radii.push_back(r);
        }


        vec1i ma_is_interior(madata.m * 2, 0);
        std::fill_n(ma_is_interior.begin(), madata.m, 1);

        output("ma_coords").set(ma_coords);
        output("ma_radii").set(ma_radii);
        output("ma_qidx").set(ma_qidx);
        output("ma_is_interior").set(ma_is_interior);
        //output("min_z").set(min_z);
    }

    void MATfilter::process() {
        std::cout << "Filter running" << std::endl;
        auto matpoints = input("ma_coords").get<PointCollection>();
        auto interior_index = input("ma_is_interior").get<vec1i>();
        auto ma_radii = input("ma_radii").get<vec1f>();
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

        std::string filepath2 = "c:\\users\\tengw\\documents\\git\\Results\\ex_indices_out.txt";
        std::ofstream outfile2(filepath2, std::fstream::out | std::fstream::trunc);

        std::string filepath1 = "c:\\users\\tengw\\documents\\git\\Results\\in_filtered_MAT.txt";
        std::ofstream outfile1(filepath1, std::fstream::out | std::fstream::trunc);
        

        for (int i = 0; i < matpoints.size(); i++) 
        {
            if (interior_index[i] == 1 )            
            {
                interior_mat.push_back(matpoints[i]);
                interior_radii.push_back(ma_radii[i]);
                interior_idx.push_back(i);
                outfile1 << matpoints[i][0] << "," << matpoints[i][1] << "," << matpoints[i][2] << "," << ma_radii[i] << "," << i << std::endl;
                //outfile2 << ma_radii[i] << std::endl;
            }
            if(interior_index[i]==0)
            {
                exterior_mat.push_back(matpoints[i]);
                exterior_radii.push_back(ma_radii[i]);
                exterior_idx.push_back(i-0.5*matpoints.size());
                outfile2 << i - 0.5*matpoints.size() << std::endl;
            }
            
        }
        outfile1.close();
        outfile2.close();

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
            std::cout << "Simplified MAT size:" << MAT_sph.size() << std::endl;
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
    }

    void ComputeNormalsNode::process() {
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

        vec3f normals_vec3f;
        normals_vec3f.reserve(madata.m);
        for (auto& n : *madata.normals) {
            normals_vec3f.push_back({ n[0], n[1], n[2] });
        }

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
       
        //std::cout << "center point of bounding box:" << center.x << "," << center.y << "," << center.z << std::endl;
        std::cout << "Total number of points in kd-tree"<< allPointsVec.size() << std::endl;
        std::cout<< "number of bounding box:" << (*kd).m_maxpoint.size() << std::endl;
        //std::cout << "Vector size of kdtreepoints:" << (*kd).m_boxKdTreePoint.size() << std::endl;
        std::cout << "size of level:" << (*kd).m_currentlevel.size() << std::endl;
        
        std::cout << "Level points are saving" << std::endl;

        
        int new_count = 0;
        for (int i = (*kd).m_maxpoint.size()-1 ; i >=0; i--)
        {

            std::vector<KdTree::sphere> levelpoints = BuildKDtree::NewGetLevelPoints((*kd).m_maxpoint[i], (*kd).m_minpoint[i], &allPointsVec);
            new_count +=levelpoints.size() ;
            //std::cout <<"Number of points in each level:"<< levelpoints.size() << std::endl;
            (*kd).m_levelpoints.push_back(levelpoints);
        }
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
        //-------output-------------------//
        std::vector<KdTree::sphere> visble_sph;

        PointCollection visible_mat;
        vec1f visible_radii;
        vec1i visible_indice;

        //----------process------------//
        // check ray intersect with box first here AMP query //
        for (int j = 0; j < headvectors.size(); j++) 
        {
            auto sph1 = AMPGPUQueryTest::GetOneLineResult(headvectors[j],endvectors[j], kd);
            visble_sph.push_back(sph1);
        }


        
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
                //std::cout << "a£º" <<a<< std::endl;
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
    void VisiblePC::process() 
    {
        std::cout << "Visible PC query starts" << std::endl;

        clock_t starttime, endtime;
        starttime = clock();
        //-------------input---------------//
        Vector3D viewpoint = input("viewPoint").get<Vector3D>();
        //auto kdtree = input("KDTree").get<KdTree*>();
        auto interior_MAT = input("interior_MAT").get<PointCollection>();
        auto radii = input("interior_radii").get<vec1f>();


        auto pc = input("original_pc").get<PointCollection>();
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
        th[3].join();*/

        

        
        


        for (int j=0;j<pc.size();j++)
        {
            bool visflag = true;
            Vector3D v2(pc[j][0], pc[j][1], pc[j][2]);
            for (int i = 0; i < interior_MAT.size(); i++)
            {
                Vector3D centre(interior_MAT[i][0], interior_MAT[i][1], interior_MAT[i][2]);
                
                float dis = DistancePointToSegment(viewpoint, v2, centre);                                
                if (dis < radii[i])
                {
                    visflag = false;
                    break;
                }                                
            }
            
            if (visflag == false)continue;
            
            visible_pc.push_back({ pc[j][0], pc[j][1], pc[j][2] });            
            
        }

        endtime = clock();

        //outfile.close();
        //----------------set result ----------------//
        std::cout << "visible pc done" << std::endl;
        std::cout << "visible points size:" << visible_pc.size() << std::endl;
        std::cout << "running time:" << endtime - starttime << std::endl;
        
        output("visible_pc").set(visible_pc);

    }
    void SightVector::process() 
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
        std::cout << "v:" << v[0] << "," << v[1] << "," << v[2] << std::endl;
        // change radius and N to dynamic later
        //float radius = 200;
        //float density = 100;

        float delta = radius / density;
        float epsilon = delta * 0.5f;
        std::cout << "Test" << std::endl;

        for (float y = -radius; y < radius + epsilon; y += delta)
        {
            for (float x = -radius; x < radius + epsilon; x += delta)
            {
                Headvectors.push_back(v1 + x * perpen_normal + y * v); // v1 is the point on the plane
                Endvectors.push_back(v2 + x * perpen_normal + y * v);
            }
        }
        std::cout << "Test done" << std::endl;
        
        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\Vectors_head_end.txt";
        std::ofstream outfileVec(filepath, std::fstream::out | std::fstream::trunc);

        for (int i = 0; i < Headvectors.size(); i++)
        {
            outfileVec << Headvectors[i][0] << "," << Headvectors[i][1] << "," << Headvectors[i][2] << ";" << Endvectors[i][0] << "," << Endvectors[i][1] << "," << Endvectors[i][2] << std::endl;
        }
        outfileVec.close();

        output("Headvectors").set(Headvectors);
        output("Endvectors").set(Endvectors);
        std::cout << "sight vectors done" << std::endl;
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
 
    


}