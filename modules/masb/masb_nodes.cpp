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
    }

    void MATfilter::process() {
        std::cout << "Filter running" << std::endl;
        auto matpoints = input("ma_coords").get<PointCollection>();
        auto interior_index = input("ma_is_interior").get<vec1i>();
        auto ma_radii = input("ma_radii").get<vec1f>();
        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\radii_out.txt";
        std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);
        for (float a : ma_radii) 
        {
            outfile << a << std::endl;
        }
        outfile.close();

        int size = matpoints.size();
        //-----output-----------------//
        PointCollection interior_mat;
        vec1f  interior_radii;

        std::string filepath2 = "c:\\users\\tengw\\documents\\git\\Results\\filtered_radii_out.txt";
        std::ofstream outfile2(filepath2, std::fstream::out | std::fstream::trunc);

        for (int i = 0; i < size; i++) {
            if (interior_index[i] == 1 && ma_radii[i] < 199)
            {
                interior_mat.push_back(matpoints[i]);
                interior_radii.push_back(ma_radii[i]);
                outfile2 << ma_radii[i] << std::endl;
            }
        }

        outfile2.close();
        std::cout << "Number of Mat points filtered:" << interior_mat.size() << std::endl;
        output("interior_mat").set(interior_mat);
        output("interior_radii").set(interior_radii);
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
        auto point_collection = input("points").get<PointCollection>();
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
        output("KDTree2").set(kd2);
    }

    void BuildKDtree::process()
    {
        auto point_collection = input("points").get<PointCollection>();
        auto radii = input("radii").get<vec1f>();
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

            Points[i].x = points[i][0];
            Points[i].y = points[i][1];
            Points[i].z = points[i][2];




            outfile<< mp_Points[i].pos.x << "," << mp_Points[i].pos.y << "," << mp_Points[i].pos.z << std::endl;
            allPointsVec.push_back(mp_Points[i]);
        }
        outfile.close();
        
        KdTree* kd = BuildKdTree(Points, m_nPoints,20);       
        //center = (*kd).centerofBoundingBox();
       
        //std::cout << "center point of bounding box:" << center.x << "," << center.y << "," << center.z << std::endl;
        std::cout<< "number of bounding box:" << (*kd).m_maxpoint.size() << std::endl;
        //std::cout << "Vector size of kdtreepoints:" << (*kd).m_boxKdTreePoint.size() << std::endl;
        std::cout << "size of level:" << (*kd).m_currentlevel.size() << std::endl;
        std::cout << "KDTree output done"<<std::endl;

        
        int new_count = 0;
        for (int i = (*kd).m_maxpoint.size()-1 ; i >=0; i--)
        {

            auto levelpoints = BuildKDtree::GetLevelPoints((*kd).m_maxpoint[i], (*kd).m_minpoint[i], allPointsVec); 
            new_count +=levelpoints.size() ;
            //std::cout <<"Number of points in each level:"<< levelpoints.size() << std::endl;
            (*kd).m_levelpoints.push_back(levelpoints);
        }


        std::cout << "Total level points:" << new_count << std::endl;
        std::cout << "The size of allPointsVec:" << allPointsVec.size() << std::endl;
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
        //(*kd).m_levelPointsIndex.resize((*kd).m_currentlevel.size());
        //int testflag = 0;
        //int count = 0;
        //for (int i = 0; i < (*kd).m_levelPointsIndex.size(); i++) {
        //    auto PointIndex = BuildKDtree::GetLevelPoints((*kd).m_boxKdTreePoint[i+1], (*kd).m_boxKdTreePoint[i]);
        //    (*kd).m_levelPointsIndex[i]= PointIndex;
        //    // whether points are inside bounding box//
        //    for (int ptindex : PointIndex) {
        //        if((mp_Points[ptindex].x<=(*kd).m_maxpoint[i].x )&&(mp_Points[ptindex].x >= (*kd).m_minpoint[i].x))
        //            if((mp_Points[ptindex].y <= (*kd).m_maxpoint[i].y) && (mp_Points[ptindex].y >= (*kd).m_minpoint[i].y))
        //                if ((mp_Points[ptindex].z <= (*kd).m_maxpoint[i].z) && (mp_Points[ptindex].z >= (*kd).m_minpoint[i].z))
        //                {
        //                    testflag++;
        //                }
        //    }
        //    //
        //    count += PointIndex.size();
        //}







        //std::vector<int> all_index;
        //for (int i = 0; i < m_nPoints; i++) {
        //    all_index.push_back(i);
        //}
        //std::vector<int> data_index;
        //for (int j = 0; j < (*kd).m_levelPointsIndex.size(); j++) {
        //    for (int a : (*kd).m_levelPointsIndex[j]) {
        //        data_index.push_back(a);
        //    }
        //} 
        //std::sort(data_index.begin(), data_index.end());
        //for (int i = 0;i< m_nPoints-1; i++) {
        //    if (data_index[i] != all_index[i]) {
        //        std::cout << "missing point index:" << i << std::endl;
        //        
        //        //no 47 directly to 48; 46 () 48
        //        break;
        //    }

        //}




        /*std::cout << all_index[46] << std::endl;
        std::cout << all_index[47] << std::endl;
        std::cout << data_index[46] << std::endl;
        std::cout << data_index[47] << std::endl;
        std::cout << points[47][0] << "," << points[47][1] << "," << points[47][2]<< "," << std::endl;*/
        
        

        //std::cout << "missing point index:" << diff_index[0] << std::endl;



        //std::cout << "total points in all level of kdtree:" << count << std::endl;
        //std::cout << "test flag:" << testflag << std::endl;
        
    
        /*std::cout << "boundingbox" << std::endl;
        for (int i = 0; i < (*kd).m_maxpoint.size(); i++) {
            std::cout << "box index:" << i << std::endl;
            std::cout << (*kd).m_maxpoint[i][0] << "," << (*kd).m_maxpoint[i][1] << "," << (*kd).m_maxpoint[i][2] << std::endl;
            std::cout << (*kd).m_minpoint[i][0] << "," << (*kd).m_minpoint[i][1] << "," << (*kd).m_minpoint[i][2] << std::endl;
        }*/
        
        
        
        //std::cout << "points in the boundingbox"<<std::endl;    
        //

        ///*int size = (*kd).m_boxKdTreePoint.size();
        //std::cout << "size:" << size << std::endl;
        //auto a0 = (*kd).m_boxKdTreePoint[1];
        //std::cout << "size:" << a0.size() << std::endl;
        //for (auto b0 : a0) {
        //    std::cout << b0  << std::endl;
        //    std::cout << mp_Points[b0].x << "," << mp_Points[b0].y << "," << mp_Points[b0].z << std::endl;
        //}*/
 
        output("KDTree").set(kd);
        
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

    void VisibiltyQurey::process() {
        std::cout << "Visibility Qurey start" << std::endl;
        //----------------input----------------//
        Vector3D viewpoint = input("ViewPoint").get<Vector3D>();
        auto kdtree = input("KDTree").get<KdTree*>();
        auto interior_mat = input("interior_MAT").get<PointCollection>();
        auto interior_radii = input("interior_radii").get<vec1f>();
        //-----------output-------------------//
        PointCollection visible_mat;

        vec1f visible_radii;

        std::string filepath = "c:\\users\\tengw\\documents\\git\\Results\\Visible_MAT_out.txt";
        std::ofstream outfile(filepath, std::fstream::out | std::fstream::trunc);
        

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
            outfile << interior_mat[i][0] << "," << interior_mat[i][1] << "," << interior_mat[i][2] << "," << interior_radii[i] << std::endl;
        }
        outfile.close();
        output("Visible_MAT").set(visible_mat);
        output("Radii_of_MAT").set(visible_radii);
        std::cout << "Visiblity Qurey Done" << std::endl;
    }
 
    


}