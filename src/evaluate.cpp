/*
 * evaluate.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: mazuran
 */

#include "evaluate.h"
#include "graph_wrapper_g2o.h"
#include "graph_wrapper_isam.h"
#include "compute_substitute_edge.h"
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include "stopwatch.h"
#include "utils.h"

/*
 * This is just a rough estimate, it absolutely doesn't mean that the final
 * memory usage will be constrained to maxMemory.
 */
const uint64_t maxMemory = 4ULL * 1024ULL * 1024ULL * 1024ULL; // Max 4GB
const uint64_t memoryMultiplier = 2;
const int sleepTime = 5;

void evaluate(GraphWrapper *gw, const EvaluateInfo &info)
{
    std::set<int> marginalized;
    GraphWrapper *incremental = gw->clonePortion(3);
    GraphWrapper *baseline    = gw->clonePortion(3);
    bool is2d = gw->vertex(1)->is2d();
    bool isg2o = true; //!info.useGLC;
    std::string saveext = (isg2o ? "g2o" : "graph");
    std::string alg = (info.algorithm == EvaluateInfo::NFR ?
            (is2d ? "se2" : "se3") :
            (info.algorithm == EvaluateInfo::GLC ? "glc" : "none")); //(isg2o ? (is2d ? "se2" : "se3") : "glc");
    std::string type, linpoint, profile, longtype, dsname;
    double kld = 0;

    switch(info.sparsityOptions.topology) {
        case SparsityOptions::Tree:            type = "tree";    longtype = "Tree";             break;
        case SparsityOptions::Subgraph:        type = "subgr";   longtype = "Subgraph";         break;
        case SparsityOptions::CliqueySubgraph: type = "clsubgr"; longtype = "Cliquey Subgraph"; break;
        case SparsityOptions::Dense:           type = "dense";   longtype = "Dense";            break;
        case SparsityOptions::CliqueyDense:    type = "cldense"; longtype = "Cliquey Dense";    break;
    }

    if(info.sparsityOptions.linPoint == SparsityOptions::Local) {
        linpoint = "l";
    } else {
        linpoint = "g";
    }

    if(info.decimate == onlineDecimate) {
        profile = "online";
    } else if(info.decimate == clusterDecimate) {
        profile = "cluster";
    } else if(info.decimate == globalDecimate) {
        profile = "global";
    } else {
        profile = "unknown";
    }

    int nstart = 0, nend = info.g2oname.length();
    size_t found1 = info.g2oname.rfind('/'), found2 = info.g2oname.rfind('.');
    if(found1 != std::string::npos) nstart = found1 + 1;
    if(found2 != std::string::npos) nend = found2;

    dsname = info.g2oname.substr(nstart, nend - nstart);

    std::string profiledir, sparsitydir, dsdir, savename;
    std::stringstream ss;
    ss << info.destdir + "/" + profile;
    profiledir = ss.str();
    ss << "/" << info.decimateOptions.sparsity;
    sparsitydir = ss.str();
    ss << "/" << dsname;
    dsdir = ss.str();
    ss << "/" << alg << "_" << type << "_" << linpoint;
    savename = ss.str();

    mkdir(info.destdir.c_str(), 0755);
    mkdir(profiledir.c_str(), 0755);
    mkdir(sparsitydir.c_str(), 0755);
    mkdir(dsdir.c_str(), 0755);
//#define G2S_DUMP_GRAPHS
#ifdef G2S_DUMP_GRAPHS
    mkdir(savename.c_str(), 0755);
#endif /* G2S_DUMP_GRAPHS */

    std::ofstream kldf((savename + ".kld").c_str()), f((savename + ".txt").c_str());
    int lastid = gw->vertices().back()->id();

    for(int i = 4; i <= lastid; i++) {
        auto latest = gw->vertex(i);

        incremental->addVertex(i, latest->estimate());
        baseline->addVertex(i, latest->estimate());

        auto edges = latest->edges();
        for(auto e: edges) {
            auto vertices = e->vertices();
            int from = vertices[0]->id(), to = vertices[1]->id();
            if(from > i || to > i) continue;
            int linkto = (vertices[0]->id() == latest->id() ? vertices[1] : vertices[0])->id();
            Eigen::MatrixXd info;
            IsometryXd meas;

            if(marginalized.count(linkto) > 0) {
                computeSubstituteEdge(gw, marginalized, i, from, to, meas, info);
            } else {
                info = e->information();
                meas = e->measurement();
            }
            incremental->addEdge(from, to, meas, info);
            baseline->addEdge(from, to, meas, info);
        }

        std::vector<int> which = info.decimate(i, lastid, info.decimateOptions);


        if(info.algorithm != EvaluateInfo::None && (
                    !which.empty() || i % info.kldPeriod == 0 || i == lastid)) {
            incremental->optimize();
            baseline->optimize();
        }

        if(!which.empty() && info.algorithm != EvaluateInfo::None) {
#ifdef G2S_VERBOSE
            std::cout << "Start sparsification (" << savename << ")" << std::endl;
#endif /* G2S_VERBOSE */
            Stopwatch s;
            s.start();
            incremental->marginalize(which, info.sparsityOptions);
            marginalized.insert(which.begin(), which.end());
            s.stop();
#ifdef G2S_VERBOSE
            std::cout << "Sparsification took " << s.time() << "s (" << savename << ")" << std::endl;
#endif /* G2S_VERBOSE */
        } else if(info.algorithm == EvaluateInfo::None) {
            marginalized.insert(which.begin(), which.end());
        }

        if(i % info.kldPeriod == 0 || i == lastid) {
            if(info.algorithm == EvaluateInfo::None) {
                baseline->optimize();
                if(info.useChi2) {
                    kld = baseline->chi2();
                    kldf << i << " " << kld << std::endl;
                } else {
                    kld = 0;
                    kldf << i << " 0" << std::endl;
                }
            } else if(info.useChi2) {
#ifdef G2S_VERBOSE
                std::cout << "Start delta chi2 computation (" << savename << ")" << std::endl;
#endif /* G2S_VERBOSE */
                Stopwatch s;
                s.start();
                kld = baseline->chi2(incremental) - baseline->chi2();
                kldf << i << " " << kld << std::endl;
                s.stop();
#ifdef G2S_VERBOSE
                std::cout << "Delta chi2 computation took " << s.time() << "s (" << savename << ")" << std::endl;
#endif /* G2S_VERBOSE */

            } else {
#ifdef G2S_VERBOSE
                std::cout << "Start KLD computation (" << savename << ")" << std::endl;
#endif /* G2S_VERBOSE */
                Stopwatch s;
                s.start();
                kld = baseline->kullbackLeibler(incremental);
                kldf << i << " " << kld << std::endl;
                s.stop();
#ifdef G2S_VERBOSE
                std::cout << "KLD computation took " << s.time() << "s (" << savename << ")" << std::endl;
#endif /* G2S_VERBOSE */
            }

#ifdef G2S_DUMP_GRAPHS
            char buf[512];
            sprintf(buf, "%s/baseline_%04d.%s", savename.c_str(), i, saveext.c_str());
            baseline->write(buf);
            sprintf(buf, "%s/marginal_%04d.%s", savename.c_str(), i, saveext.c_str());
            incremental->write(buf);
#endif /* G2S_DUMP_GRAPHS */
        }
    }

    std::transform(alg.begin(), alg.end(), alg.begin(), ::toupper);

    f << alg << " " << longtype;
    f << std::endl << "    baseline:     ";
    baseline->printStats(f);
    f << std::endl << "    marginalized: ";
    incremental->printStats(f);
    f << std::endl << "    last " << (info.useChi2 ? "chi2: " : "kld: ") << kld << std::endl;

//#define G2S_DUMP_INFOMAT
#ifdef G2S_DUMP_INFOMAT
    std::vector<Eigen::Triplet<double> > t = triplets(static_cast<GraphWrapperG2O *>(baseline)->sparseInformation());
    std::ofstream finfo((savename + ".info").c_str());
    for(const Eigen::Triplet<double> &triplet: t) {
        finfo << triplet.row() << " " << triplet.col() << " " << triplet.value() << std::endl;
    }
    finfo.close();
#endif /* G2S_DUMP_INFOMAT */

    kldf.close();
    f.close();

    delete incremental;
    delete baseline;
}

struct ThreadData
{
    std::vector<EvaluateInfo> joblist;
    std::vector<EvaluateInfo*> running;
    uint64_t memoryRequirement;
    std::map<std::string, int> freeVariables;
    pthread_mutex_t mutex;
};

int countFreeVariables(const std::string &fname)
{
    int variables = 0;
    std::vector<std::string> prefixes = { "VERTEX_SE2", "VERTEX_SE3" };
    std::vector<int> sizes = { 3, 6 };

    std::ifstream f(fname.c_str());
    while(f.good()) {
        std::string line;
        std::getline(f, line);
        for(int i = 0; i < prefixes.size(); i++) {
            if(line.compare(0, prefixes[i].length(), prefixes[i]) == 0) {
                variables += sizes[i];
                break;
            }
        }
    }
    return variables;
}

void printSelected(const EvaluateInfo *thisJob)
{
    std::string profile, topology, linPoint;

    if(thisJob->decimate == onlineDecimate) {
        profile = "online";
    } else if(thisJob->decimate == clusterDecimate) {
        profile = "clustered";
    } else {
        profile = "global";
    }

    if(thisJob->sparsityOptions.topology == SparsityOptions::Tree) {
        topology = "tree";
    } else if(thisJob->sparsityOptions.topology == SparsityOptions::Subgraph) {
        topology = "subgraph";
    } else if(thisJob->sparsityOptions.topology == SparsityOptions::CliqueySubgraph) {
        topology = "cliquey subgraph";
    } else if(thisJob->sparsityOptions.topology == SparsityOptions::Dense) {
        topology = "dense";
    } else {
        topology = "cliquey dense";
    }

    if(thisJob->sparsityOptions.linPoint == SparsityOptions::Local) {
        linPoint = "local";
    } else {
        linPoint = "global";
    }


    std::cout << "Running configuration\n";
    std::cout << "  name: " << thisJob->g2oname << "\n";
    std::cout << "  glc: " << (thisJob->algorithm == EvaluateInfo::GLC ? "true" : "false") << "\n";
    std::cout << "  profile: " << profile << "\n";
    std::cout << "  topology: " << topology << "\n";
    std::cout << "  linPoint: " << linPoint << "\n";
    std::cout << "  sparsity: " << thisJob->decimateOptions.sparsity << "\n";
    std::cout << "  eval: " << (thisJob->useChi2 ? "chi2" : "kld") << "\n";
    std::cout << std::flush;
}

void saveMissingRuns(const std::vector<EvaluateInfo> &joblist)
{
    std::stringstream ss;
    ss << "missing-" << getpid() << ".txt";
    std::ofstream f(ss.str().c_str());
    for(const EvaluateInfo &info: joblist) {
        if(info.algorithm == EvaluateInfo::GLC) {
            f << "glc ";
        } else if(info.algorithm == EvaluateInfo::NFR) {
            f << "sen ";
        } else {
            f << "none ";
        }

        f << info.g2oname;

        if(info.decimate == onlineDecimate) {
            f << " online";
        } else if(info.decimate == clusterDecimate) {
            f << " cluster";
        } else {
            f << " global";
        }

        if(info.sparsityOptions.topology == SparsityOptions::Tree) {
            f << " tree";
        } else if(info.sparsityOptions.topology == SparsityOptions::Subgraph) {
            f << " subgr";
        } else if(info.sparsityOptions.topology == SparsityOptions::CliqueySubgraph) {
            f << " clsubgr";
        } else if(info.sparsityOptions.topology == SparsityOptions::Dense) {
            f << " dense";
        } else {
            f << " cldense";
        }

        if(info.sparsityOptions.linPoint == SparsityOptions::Local) {
            f << " local";
        } else {
            f << " global";
        }

        f << " " << info.decimateOptions.sparsity << " " << info.kldPeriod << " " <<
                info.decimateOptions.clusterSize << std::endl;
    }
    f.close();
}

void *threadedRunner(void *arg)
{
    EvaluateInfo *thisJob = NULL;
    uint64_t memoryRequirement = 0;
    ThreadData *data = static_cast<ThreadData *>(arg);
    while(1) {
        pthread_mutex_lock(&data->mutex);
        data->memoryRequirement -= memoryRequirement;

        if(thisJob != NULL) {
            data->running.erase(std::find(data->running.begin(), data->running.end(), thisJob));
            delete thisJob;
            thisJob = NULL;

            std::vector<EvaluateInfo> joblistfull;
            for(const EvaluateInfo *info: data->running) {
                joblistfull.push_back(*info);
            }
            joblistfull.insert(joblistfull.end(), data->joblist.begin(), data->joblist.end());
            //FIXME: save missing runs only if needed
            //saveMissingRuns(data->joblist);
        }

        for(int i = 0; i < data->joblist.size(); i++) {
            EvaluateInfo job = data->joblist[i];
            uint64_t freeVariables;

            if(data->freeVariables.count(job.g2oname) > 0) {
                freeVariables = data->freeVariables[job.g2oname];
            } else {
                freeVariables = countFreeVariables(job.g2oname);
            }

            memoryRequirement = memoryMultiplier * sizeof(double) *
                    freeVariables * freeVariables /
                    job.decimateOptions.sparsity /
                    job.decimateOptions.sparsity;

            if(data->memoryRequirement + memoryRequirement <= maxMemory) {
                thisJob = new EvaluateInfo(job);
                data->joblist.erase(data->joblist.begin() + i);
                break;
            }
        }

        if(thisJob == NULL && data->joblist.size() == 0) {
            pthread_mutex_unlock(&data->mutex);
            return NULL;
        }

        if(thisJob) {
            printSelected(thisJob);
            data->memoryRequirement += memoryRequirement;
            data->running.push_back(thisJob);
        }

        pthread_mutex_unlock(&data->mutex);

        if(thisJob == NULL) {
            memoryRequirement = 0;
            sleep(sleepTime);
        } else {
            GraphWrapperG2O *g2o = new GraphWrapperG2O(
                    thisJob->g2oname.c_str(), true, thisJob->algorithm == EvaluateInfo::GLC);
            evaluate(g2o, *thisJob);
            delete g2o;
        }
    }
    return NULL;
}

void parallelEvaluate(const std::vector<EvaluateInfo> &joblist, int threads)
{
    if(threads <= 0) {
        threads = sysconf(_SC_NPROCESSORS_ONLN);
    }

    pthread_t handles[threads];
    ThreadData data = { joblist, std::vector<EvaluateInfo*>(), 0, std::map<std::string, int>(), pthread_mutex_t() };
    pthread_mutex_init(&data.mutex, NULL);

    for(int i = 0; i < threads; i++) {
        pthread_create(&handles[i], NULL, threadedRunner, &data);
    }

    for(int i = 0; i < threads; i++) {
        pthread_join(handles[i], NULL);
    }

    pthread_mutex_destroy(&data.mutex);

}
