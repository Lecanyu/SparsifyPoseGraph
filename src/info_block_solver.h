#ifndef INFO_BLOCK_SOLVER_H_
#define INFO_BLOCK_SOLVER_H_

#include <g2o/core/block_solver.h>

template <typename Traits>
class InfoBlockSolver : public g2o::BlockSolver<Traits>
{
public:
    static const int PoseDim = Traits::PoseDim;
    static const int LandmarkDim = Traits::LandmarkDim;
    typedef typename Traits::PoseMatrixType PoseMatrixType;
    typedef typename Traits::LandmarkMatrixType LandmarkMatrixType; 
    typedef typename Traits::PoseLandmarkMatrixType PoseLandmarkMatrixType;
    typedef typename Traits::PoseVectorType PoseVectorType;
    typedef typename Traits::LandmarkVectorType LandmarkVectorType;

    typedef typename Traits::PoseHessianType PoseHessianType;
    typedef typename Traits::LandmarkHessianType LandmarkHessianType;
    typedef typename Traits::PoseLandmarkHessianType PoseLandmarkHessianType;
    typedef typename Traits::LinearSolverType LinearSolverType;

public:
    InfoBlockSolver(LinearSolverType *linearSolver) : g2o::BlockSolver<Traits>(linearSolver)
    {
    }

    ~InfoBlockSolver()
    {
    }

    const g2o::SparseBlockMatrix<Eigen::MatrixXd> &information()
    {
        return *g2o::BlockSolver<Traits>::_Hpp;
    }
};

#endif /* INFO_BLOCK_SOLVER_H_ */
