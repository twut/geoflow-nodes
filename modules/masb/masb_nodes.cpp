#include "masb_nodes.hpp"
#include "iostream"
#include <fstream>
#include <sstream>
#include <algorithm>



namespace geoflow::nodes::mat {

    void ComputeMedialAxisNode::process() {
        auto point_collection = input("points").get<PointCollection>();
        auto normals_vec3f = input("normals").get<vec3f>();

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
        vec1f ma_radii(madata.m * 2);
        for (size_t i = 0; i < madata.m * 2; ++i) {
            ma_radii.push_back(Vrui::Geometry::dist(coords[i%madata.m], ma_coords_[i]));
        }


        vec1i ma_is_interior(madata.m * 2, 0);
        std::fill_n(ma_is_interior.begin(), madata.m, 1);

        output("ma_coords").set(ma_coords);
        output("ma_radii").set(ma_radii);
        output("ma_qidx").set(ma_qidx);
        output("ma_is_interior").set(ma_is_interior);
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

    void testNode::process() 
    {
        //filepath = "C:\Users\tengw\Documents\Git\test.las";
        std::cout << "Test 1 running" << std::endl;
        auto radii = input("in_radii").get<vec1f>();
        //std::cout <<"Filepath: "<< filepath << std::endl;
        char filename[256] = "C:\\Users\\tengw\\Documents\\Git\\radii.txt";
        //std::cout <<"Writing to file:" <<filename << std::endl;
        for (auto a = radii.begin(); a != radii.end(); ++a) 
        {
            std::cout << "radius:" << *a << std::endl;
        }
             

    }

    void BuildKDtree::process()
    {
        auto point_collection = input("points").get<PointCollection>();
        masb::ma_data madata;
        madata.m = point_collection.size();
        m_nPoints = point_collection.size();
        std::cout << "MAT point size:" << madata.m << std::endl;
        Vector3D* mp_Points = new Vector3D[m_nPoints];

        vec3f points;
        points.reserve(madata.m);
        for (auto& p : point_collection) {
            points.push_back({ p[0], p[1], p[2] });
        }
        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\MAT_points_out.txt";
        std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);
        for (int i = 0; i < m_nPoints; i++) {
            mp_Points[i].x = points[i][0];
            mp_Points[i].y = points[i][1];
            mp_Points[i].z = points[i][2];
            outfile<< mp_Points[i].x << "," << mp_Points[i].y << "," << mp_Points[i].z << std::endl;
        }
        outfile.close();
        /*std::cout << "Test X: " << mp_Points[18655].x << std::endl;
        std::cout << "Test X: " << mp_Points[18655].y << std::endl;
        std::cout << "Test X: " << mp_Points[18656].z << std::endl;*/
        KdTree* kd = BuildKdTree(mp_Points, m_nPoints,16);                    
        output("KDTree").set(kd);
        std::cout << "KDTree output done"<<std::endl;
                          
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

}