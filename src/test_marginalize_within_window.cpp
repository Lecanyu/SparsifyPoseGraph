//
// Created by lecanyu on 18-11-2.
//

#include <iostream>
#include <unordered_set>
#include "graph_wrapper_g2o.h"
#include "flameprofiler.h"

struct MyNode2d{
    int id_;
    g2o::SE2 pose_;
    MyNode2d():id_(-1), pose_(g2o::SE2()){}
    MyNode2d(int id, const g2o::SE2& pose):id_(id), pose_(pose){}
};

struct MyEdge2d{
    int from_, to_;
    g2o::SE2 odometry_;   // relative transformation
    Eigen::MatrixXd info_;
    MyEdge2d():from_(-1),to_(-1),odometry_(g2o::SE2()),info_(Eigen::MatrixXd::Identity(3, 3)){}
    MyEdge2d(int from, int to, const g2o::SE2& odometry, const Eigen::MatrixXd& info)
        : from_(from),to_(to),odometry_(odometry),info_(info){}
};

struct MyEdge2dHash{
    int operator()(const MyEdge2d& edge) const{
        return edge.from_*10000 + edge.to_;
    }
};

struct MyEdge2dEq{
    bool operator()(const MyEdge2d& e1, const MyEdge2d& e2) const {
        return e1.from_==e2.from_ && e1.to_==e2.to_;
    }
};

int main(int argc, char** argv) {
    if(argc != 6)
    {
        std::cout<<"./test save_path window_len loop_closure_times use_glc use_tree2sparsify\n";
        return 0;
    }
    std::string save_path = argv[1];
    int window_len = atoi(argv[2]);
    int loop_closure_times = atoi(argv[3]);
    int use_glc = atoi(argv[4]);
    int use_tree = atoi(argv[5]);

//    std::string save_path = "/home/lecanyu/Desktop/sparsifier_test";
//    int window_len = 15;
//    int loop_closure_times = 8;
//    int use_glc = 0;

    // NFR or GLC
    std::string opt_type;
    bool glc;
    if(use_glc==0) {
        glc = false;
        opt_type = "NFR";
    }
    else {
        glc = true;
        opt_type = "GLC";
    }

    // use tree to sparsify
    std::string sparse_type;
    if(use_tree)
        sparse_type = "tree";
    else
        sparse_type = "dense";

    std::vector<int> nodeIds;
    GraphWrapperG2O slam(false, glc);
    /**
     * Canyu Le
     * You can use below code to load a graph.
     * */
//    slam.getOptimizer()->load("/home/lecanyu/Desktop/slam_status.g2o");

    // initial nodes and edges
    std::vector<MyNode2d> nodes;
    std::unordered_set<MyEdge2d, MyEdge2dHash, MyEdge2dEq> edges;

    g2o::SE2 odo(1., 0., 0.);
    Eigen::MatrixXd noise = 0.01 * Eigen::MatrixXd::Identity(3, 3);
    noise(2, 2) = 100;

    for (int i = 0; i < window_len; ++i) {
        MyNode2d node(i, g2o::SE2());
        nodes.push_back(node);
        nodeIds.push_back(i);
    }
    for (int i = 1; i < window_len; ++i) {
        MyEdge2d edge(i - 1, i, odo, noise);
        edges.insert(edge);
    }

    for (auto node : nodes)
        slam.addVertex(node.id_, node.pose_);
    for (auto edge : edges)
        slam.addEdge(edge.from_, edge.to_, edge.odometry_, edge.info_);

    slam.optimize();

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, window_len - 2);

    std::stringstream stat_filename_s, save_filename_s;
    stat_filename_s<< save_path <<"/stats_"<<opt_type <<"_"<< sparse_type <<"_win_"<<window_len<<"_loop_"<<loop_closure_times<<".txt";
    save_filename_s<<"time_"<<opt_type <<"_"<< sparse_type <<"_win_"<<window_len<<"_loop_"<<loop_closure_times<<".json";

    std::ofstream file2(stat_filename_s.str());

    for (int i = 0; i < 200000; ++i)
    {
        // add new node and edge
        int last_node_id = nodeIds.back();
        int new_node_id = last_node_id + 1;
        nodeIds.push_back(new_node_id);

        slam.addVertex(new_node_id, g2o::SE2(i + window_len, 0, 0));
        slam.addEdge(last_node_id, new_node_id, odo, noise);


        // add random loop closure edge between current node and node in the sliding window
        std::unordered_set<MyEdge2d, MyEdge2dHash, MyEdge2dEq> loop_edge_set;
        for (int j = 0; j < loop_closure_times; ++j)
        {
            // a uniform random node in the sliding window
            int previous_node_index = dis(gen);
            int previous_node_id = nodeIds[previous_node_index];
            int current_node_id = nodeIds.back();

            MyEdge2d e(previous_node_id, current_node_id, g2o::SE2(), Eigen::MatrixXd::Identity(3, 3));
            if(loop_edge_set.find(e)==loop_edge_set.end())
            {
                loop_edge_set.insert(e);
                slam.addEdge(previous_node_id, current_node_id, g2o::SE2(nodeIds.size() - previous_node_index - 1, 0., 0.), noise);
            }
        }


        {
            PZonePath("marg_and_opt_time", "default", save_path.c_str(), save_filename_s.str().c_str());
            // marginalize to keep fixed length window
            {
                PZonePath("marg_time", "default", save_path.c_str(), save_filename_s.str().c_str());
                std::vector<int> marg;
                marg.push_back(nodeIds.front());
                // most accurate sparsification, but a little bit slow.
                SparsityOptions opts;
                opts.linPoint = SparsityOptions::Global;
                if (use_tree)
                    opts.topology = SparsityOptions::Tree;
                else
                    opts.topology = SparsityOptions::Dense;
                slam.marginalizeNoOptimize(marg, opts);

                nodeIds.erase(nodeIds.begin());
            }

            {
                PZonePath("opt_time", "default", save_path.c_str(), save_filename_s.str().c_str());
                slam.optimizeFromId(nodeIds.front());
            }

            {
                PZonePath("current_pose_cov_calc_time", "default", save_path.c_str(), save_filename_s.str().c_str());
                Eigen::MatrixXd cov = slam.covariance();
            }
        }

//        std::string stats = slam.getStats();
//        file2 << stats << "\n";
//        std::cout<<stats<<"\n";
//        slam.debugPrint(std::cout);
//        std::cout<<std::endl;
//        ProfilerWriteFileHeadToDisk();
        ProfilerWriteFileHeadToDiskWithPath(save_path, save_filename_s.str());


//        std::vector<double> result;
//        if(slam.getOptimizer()->vertex(new_node_id)->getEstimateData(result)){
//            std::cout<<"vertex "<<new_node_id<<": ";
//            for(auto r:result)
//                std::cout<<r<<", ";
//            std::cout<<"\n";
//        }
    }

    /**
     * Canyu Le
     * You can use below codes to save a graph.
     *
     * If you have questions about the format (i.e. the GLC factor), go to "glc_edge.cpp" -> write(std::ostream& os) to figure out the meaning.
     * I'd like to add a simple explanation here for glc_edge format:
     * GLC_EDGE 10 16 || GLC_REPARAM_SE2_ISAM 3 6 10 0 0 6 0 0 0 6.44942e-14 7.37069e-14 -5.38607e-12 -0.0781786 -2.09314e-05 -0 -2.5976e-16 -3.96149e-14 0.0782534 -5.39122e-12 -1.42468e-15 0 5.67138e-15 -1.09704e-14 -1.87717e-15 -0.00209596 7.82839 1 0 0 1 0 1
     * 10 16: connected two vertices.
     * 3 6: dimension of matrix _W
     * 10 0 0: the pose of vertex 10 .
     * 6 0 0: the relative transformation from vertex 10 to 16.
     * 0 6.44942e-14 ... 7.82839: the 3x6 elements in matrix _W.
     * 1 0 0 1 0 1: information matrix (upper triangular matrix).
     * */
//    std::ofstream f("/home/lecanyu/Desktop/slam_status.g2o");
//    slam.write(f);

    return 0;
}