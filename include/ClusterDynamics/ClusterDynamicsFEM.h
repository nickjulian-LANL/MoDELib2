/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsFEM_H_
#define model_ClusterDynamicsFEM_H_


#include <ClusterDynamicsParameters.h>
#include <DislocationDynamicsBase.h>
#include <EvalFunction.h>
#include <FixedDirichletSolver.h>
#include <DDconfigIO.h>
#include <MicrostructureBase.h>
#include <SecondOrderReaction.h>
#include <MicrostructureContainer.h>
namespace model
{

    template <int dim>
    struct FluxMatrix : public EvalFunction<FluxMatrix<dim>>
    {
        typedef typename DislocationDynamicsBase<dim>::ElementType ElementType;
        typedef Eigen::Matrix<double,dim+1,1> BaryType;
        constexpr static int mSize=ClusterDynamicsParameters<dim>::mSize;
        constexpr static int rows=dim*mSize;
        constexpr static int cols=rows;
        typedef Eigen::Matrix<double,rows,cols> MatrixType;

        const ClusterDynamicsParameters<dim>& cdp;
        
        FluxMatrix(const ClusterDynamicsParameters<dim>& cdp_in);
        const MatrixType operator() (const ElementType& elem, const BaryType& bary) const;
    };

    template<int dim>
    struct ClusterDynamicsFEM
    {
        typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
        typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
        typedef typename MicrostructureBase<dim>::ElementType ElementType;
        typedef typename MicrostructureBase<dim>::SimplexDim SimplexDim;
        typedef typename MicrostructureBase<dim>::NodeType NodeType;
        typedef typename MicrostructureBase<dim>::VectorMSize VectorMSize;
        typedef FiniteElement<ElementType> FiniteElementType;
        static constexpr int dVorder=4;
        typedef IntegrationDomain<FiniteElementType,0,dVorder,GaussLegendre> VolumeIntegrationDomainType;
        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;
        static constexpr int iSize=ClusterDynamicsParameters<dim>::iSize;
        typedef TrialFunction<'m',mSize,FiniteElementType> MobileTrialType;
        typedef TrialFunction<'i',iSize,FiniteElementType> ImmobileTrialType;
//        typedef TrialFunction<'z',dim,FiniteElementType> DiffusiveTrialType;
        typedef TrialGrad<MobileTrialType> MobileGradType;
        typedef TrialProd<FluxMatrix<dim>,MobileGradType> MobileFluxType;
        typedef BilinearForm<MobileGradType,TrialProd<Constant<double,1,1>,MobileFluxType>> MobileBilinearFormType;
        typedef BilinearWeakForm<MobileBilinearFormType,VolumeIntegrationDomainType> MobileBilinearWeakFormType;

        typedef TrialFunction<'d',mSize,FiniteElementType> MobileIncrementTrialType;
        typedef TrialGrad<MobileIncrementTrialType> MobileIncrementGradType;
        typedef TrialProd<FluxMatrix<dim>,MobileIncrementGradType> MobileIncrementFluxType;
        typedef BilinearForm<MobileIncrementGradType,TrialProd<Constant<double,1,1>,MobileIncrementFluxType>> MobileIncrementBilinearFormType;
        typedef BilinearWeakForm<MobileIncrementBilinearFormType,VolumeIntegrationDomainType> MobileIncrementBilinearWeakFormType;

        typedef Eigen::SparseMatrix<double> SparseMatrixType;
    #ifdef _MODEL_PARDISO_SOLVER_
        typedef Eigen::PardisoLLT<SparseMatrixType> DirectSPDSolverType;
        typedef Eigen::PardisoLU<SparseMatrixType> DirectSquareSolverType;
    #else
        typedef Eigen::SimplicialLLT<SparseMatrixType> DirectSPDSolverType;
        typedef Eigen::SparseLU<SparseMatrixType> DirectSquareSolverType;
    #endif
        typedef Eigen::ConjugateGradient<SparseMatrixType> IterativeSPDSolverType;
        typedef Eigen::BiCGSTAB<SparseMatrixType> IterativeSquareSolverType;


        const DislocationDynamicsBase<dim>& ddBase;
        const ClusterDynamicsParameters<dim>& cdp;
        
        MobileTrialType mobileClusters;
        MobileGradType mobileGrad;
        MobileFluxType mobileFlux;
        ImmobileTrialType immobileClusters;

        const int nodeListInternalExternal;
        MobileIncrementTrialType mobileClustersIncrement;
        VolumeIntegrationDomainType dV;
        MobileBilinearWeakFormType mBWF;
        MobileIncrementBilinearWeakFormType dmBWF;

//        FixedDirichletSolver<MobileBilinearWeakFormType> mSolver;
        FixedDirichletSolver mSolver;
        bool solverInitialized;

        const Eigen::VectorXd cascadeGlobalProduction;



                


        ClusterDynamicsFEM(const DislocationDynamicsBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in);
        void solveMobileClusters();
        void solveImmobileClusters();
        void solve();
//        void applyBoundaryConditions();
        void initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels);
        void initializeSolver();
        VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const;

    };
    
}
#endif

