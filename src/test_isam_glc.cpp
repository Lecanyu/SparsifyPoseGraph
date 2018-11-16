/**
 * @file glc.cpp
 * @brief 2D Generic Linear Constraints Example with poses and a landmark
 * @author Nicholas Carlevaris-Bianco
 * @author Michael Kaess
 * @version $Id: glc.cpp 8578 2013-07-01 00:28:49Z kaess $
 *
 * Copyright (C) 2009-2013 Massachusetts Institute of Technology.
 * Michael Kaess, Hordur Johannsson, David Rosen,
 * Nicholas Carlevaris-Bianco and John. J. Leonard
 *
 * This file is part of iSAM.
 *
 * iSAM is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * iSAM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with iSAM.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <ctime>
#include <isam/Slam.h>
#include <isam/slam2d.h>
#include <isam/glc.h>
#include <unordered_set>
#include <set>

#include <fstream>

#include "flameprofiler.h"

using namespace std;
using namespace isam;
using namespace Eigen;

MatrixXd get_info(Slam *slam)
{
    const SparseSystem &R = slam->get_R();
    const int *a_to_r = R.a_to_r(); // column ordering
    MatrixXd Rsub(R.num_rows(), R.num_cols());
    // get colums of R for each of the inds
    for (int i = 0; i < R.num_rows(); i++)
    {
        for (int j = 0; j < R.num_cols(); j++)
        {
            Rsub(i, j) = R(a_to_r[i], a_to_r[j]);
        }
    }

    return Rsub.transpose() * Rsub;
}

void print_all(Slam &slam)
{
    list<Node *> ids = slam.get_nodes();
    for (list<Node *>::iterator iter = ids.begin(); iter != ids.end(); iter++)
    {
        Pose2d_Node *pose = dynamic_cast<Pose2d_Node *>(*iter);
        cout << pose->value() << endl;
    }
}

double snrv(void)
{
    srand(time(0));
    double U = rand() / double(RAND_MAX);
    double V = rand() / double(RAND_MAX);
    // Box-Muller method
    return sqrt(-2.0 * log(U)) * cos(2.0 * M_PI * V);
}

Pose2d add_noise(Pose2d in, double sigma)
{
    Pose2d out(in.x() + sigma * snrv(),
               in.y() + sigma * snrv(),
               in.t() + sigma * snrv());
    return out;
}

Point2d add_noise(Point2d in, double sigma)
{
    Point2d out(in.x() + sigma * snrv(),
                in.y() + sigma * snrv());
    return out;
}


void print_block_cov(const list<Node *> &nodes, const Covariances &covariances)
{
    Covariances::node_lists_t node_lists;
    for (auto & node : nodes)
    {
        node_lists.push_back(list<Node *>{node});
    }
    list<MatrixXd> cov_blocks = covariances.marginal(node_lists);
    int i = 1;
    for (list<MatrixXd>::iterator it = cov_blocks.begin(); it != cov_blocks.end(); it++, i++)
    {
        cout << "block " << i << ":" <<endl;
        cout << *it << endl;
    }
}

void print_cov(const Covariances &covariances)
{

}

// factor set definition
struct FactorHash
{
    size_t operator()(Factor * x) const
    {
        return std::hash<int>()(x->unique_id());
    }
};
struct FactorEqualTo
{
    bool operator()(Factor *factor1, Factor *factor2) const
    {
        return (factor1->unique_id() == factor2->unique_id());
    }
};

typedef std::unordered_set<Factor *, FactorHash, FactorEqualTo> FactorSet;

// ordered factor set definition
struct FactorLess
{
    bool operator()(Factor *factor1, Factor *factor2) const
    {
        return (factor1->unique_id() < factor2->unique_id());
    }
};

typedef std::set<Factor *, FactorLess> FactorOrderedSet;

void print_factor_set(const FactorSet &factor_set)
{
    FactorOrderedSet ordered_set;
    for(auto &factor:factor_set)
    {
        ordered_set.insert(factor);
    }
    for(auto &factor:ordered_set)
    {
        cout << *factor << "\t" << factor->unique_id() <<endl;
    }
}

int main(int argc, char** argv)
{
    if(argc != 5)
    {
        std::cout<<"./test save_path window_len loop_closure_times use_tree2sparsify\n";
        return 0;
    }
    std::string save_path = argv[1];
    int window_len = atoi(argv[2]);
    int loop_closure_times = atoi(argv[3]);
    int use_tree = atoi(argv[4]);

    // use tree to sparsify
    bool sparse;
    std::string sparse_type;
    if(use_tree)
    {
        sparse = true;
        sparse_type = "tree";
    }
    else
    {
        sparse = false;
        sparse_type = "dense";
    }


    // locally defined, as nodes and factors below are allocated on the
    // stack, and therefore will go out of scope after the end of this
    // function; note that you can (and typically have to) dynamically
    // allocate, but I used local variables here to make the example
    // simpler - see example.cpp for dynamic allocation.
    Slam slam;
    Properties properties;
    properties.mod_batch = 5;
    slam.set_properties(properties);

    Eigen::MatrixXd cov = 0.01 * eye(3);
    cov(2,2) = 0.000000001;
    Noise noise = Covariance(cov);

    std::vector<Node *> node_vec;
    FactorSet factor_set;

    Pose2d prior(0., 0., 0.);
    auto *x0 = new Pose2d_Node;
    node_vec.push_back(x0);
    slam.add_node(x0);
    auto *p_x0 =  new Pose2d_Factor(x0, prior, noise);
    slam.add_factor(p_x0);
    factor_set.insert(p_x0);

    Pose2d odo(1., 0., 0.);

    for (int i = 0; i < window_len-1; ++i)
    {
        auto *last_node = dynamic_cast<Pose2d_Node *>(node_vec.back());
        auto *current_node = new Pose2d_Node;
        node_vec.push_back(current_node);
        slam.add_node(current_node);
        auto *odom_factor = new Pose2d_Pose2d_Factor( last_node, current_node, odo, noise);
        slam.add_factor(odom_factor);
        factor_set.insert(odom_factor);
    }
//    cout << "initial graph state" << endl;
//    slam.print_graph();
    slam.batch_optimization();
//    slam.print_stats();
//    slam.print_graph();

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, window_len-2);

    std::stringstream stat_filename_s, save_filename_s;
    save_filename_s<<"time_isam_" << sparse_type <<"_win_"<<window_len<<"_loop_"<<loop_closure_times<<".json";

//    ofstream file2("graph_localization_optimization_statistics.txt");

    for (int i = 0; i < 10000; ++i)
    {
//        cout << "--------------- fixed lag smoothing: " << i << " ----------------" <<endl;

        // add new node and odometry factor
        {
            auto *last_node = dynamic_cast<Pose2d_Node *>(node_vec.back());
            auto *current_node = new Pose2d_Node;
            node_vec.push_back(current_node);
            slam.add_node(current_node);
            auto *odom_factor = new Pose2d_Pose2d_Factor( last_node, current_node, odo, noise);
            slam.add_factor(odom_factor);
            factor_set.insert(odom_factor);
        }

        // add prior factor for the new node
        {
            Pose2d prior_pose(i+window_len, 0., 0.);
            auto *current_node = dynamic_cast<Pose2d_Node *>(node_vec.back());
            auto *prior_factor_for_current_node =  new Pose2d_Factor(current_node, prior_pose, noise);
            slam.add_factor(prior_factor_for_current_node);
            factor_set.insert(prior_factor_for_current_node);
        }

        // add random loop closure edge between current node and node in the sliding window
        for (int j = 0; j < loop_closure_times; ++j)
        {
            // a uniform random node in the sliding window
            int node_index = dis(gen);
            auto *node_to_loop_with = dynamic_cast<Pose2d_Node *>(node_vec[node_index]);

            auto *current_node = dynamic_cast<Pose2d_Node *>(node_vec.back());
            Pose2d relative_pose( node_vec.size() - node_index - 1, 0., 0.);
            auto *loop_factor = new Pose2d_Pose2d_Factor( node_to_loop_with, current_node, relative_pose, noise);
            slam.add_factor(loop_factor);
            factor_set.insert(loop_factor);
        }


        {
            PZonePath("marg_and_opt_time", "default", save_path.c_str(), save_filename_s.str().c_str());
//            PZone("marg_and_opt_time");
            // marginalize to keep fixed length window
            {
                PZonePath("marg_time", "default", save_path.c_str(), save_filename_s.str().c_str());
//                PZone("marg_time");
                vector<Factor*> factor_elim = glc_elim_factors (node_vec.front()); //usefull for local managment of factors
                vector<Factor*> fnew = glc_remove_node (slam, node_vec.front(), sparse); // not root shifted

                for(auto &factor:fnew)
                {
                    factor_set.insert(factor);
                }

                // delete removed factors
                for (auto &factor:factor_elim)
                {
                    factor_set.erase(factor);
                    delete factor;
                }

                // delete removed node
                delete node_vec.front();
                node_vec.erase(node_vec.begin());
            }

            {
                PZonePath("opt_time", "default", save_path.c_str(), save_filename_s.str().c_str());
//                PZone("opt_time");
                slam.batch_optimization();
//                slam.update();
            }

            {
                {
                    PZonePath("current_pose_cov_calc_time", "default", save_path.c_str(), save_filename_s.str().c_str());
//                    PZone("current_pose_cov_calc_time");
                    const Covariances &covariances = slam.covariances().clone();

//                    cout << "covariances: " << endl;
                    list<Node *> nodes;
                    nodes.push_back(node_vec.back());
                    Covariances::node_lists_t node_lists;
                    node_lists.push_back(nodes);
                    // recover covariance matrix using R matrix, very efficient for most recent states
                    // as described in the iSAM paper.
                    list<MatrixXd> cov_blocks = covariances.marginal(node_lists);
                }

//                for (list<MatrixXd>::iterator it = cov_blocks.begin(); it != cov_blocks.end(); it++)
//                {
//                    cout << "current pose cov: " <<endl;
//                    cout << *it << endl;
//                }
            }
        }

//        slam.print_stats();
//        slam.print_graph();

//        auto R = slam.get_R();
//        double nnz = slam.get_R().nnz();
//        double max_per_col = slam.get_R().max_nz();
//        double dim = slam.get_nodes().size()*3;
//        double per_col = nnz / dim;
//        double fill_in = nnz / (dim * dim);

//        file2 << max_per_col << "\t" << per_col << "\t" << fill_in*0.01 << endl;

//        {
//            const Covariances &covariances = slam.covariances().clone();
//
//            // recovering the full covariance matrix
//            cout << "Full covariance matrix:" << endl;
//            MatrixXd cov_full = covariances.marginal(slam.get_nodes());
//            cout << cov_full << endl << endl;
//
//            // recovering the block-diagonals only of the full covariance matrix
//            cout << "Block-diagonals only:" << endl;
//            print_block_cov(node_vec, covariances);
//        }

        ProfilerWriteFileHeadToDisk();
    }

    // free all resources
    cout << "final factor_set size: " << factor_set.size() << endl;
    for (auto &factor:factor_set)
    {
        delete factor;
    }
    cout << "final node_vec size: " << node_vec.size() << endl;
    for(auto &node:node_vec)
    {
        delete node;
    }

    return 0;
}
