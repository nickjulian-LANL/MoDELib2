/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_MicrostructureSpecification_CPP_
#define model_MicrostructureSpecification_CPP_

#include <MicrostructureSpecification.h>
namespace model
{
   MicrostructureSpecification::MicrostructureSpecification(
            const std::string& microstructureTypeIn,
            const std::string& styleIn,
            const std::string& tagIn,
            const DislocationDynamicsBase<3>::DislocationDynamicsBaseType& ddBaseIn,
            // periodicDipoles
            const std::vector<int>& periodicDipoleSlipSystemIDsIn,
            const std::vector<int>& periodicDipoleExitFaceIDsIn, //  #  each value in the vector is the ID of the face from which the dipole enters/exits
            const Eigen::Matrix<double, Eigen::Dynamic,3>& periodicDipolePointsIn, // # each row in the matrix is a the "center" of the dipole
            const std::vector<double>& periodicDipoleHeightsIn, // # each value in the vector is the height of the dipole in number of slip planes
            const std::vector<int>& periodicDipoleNodesIn,
            const std::vector<double>& periodicDipoleGlideStepsIn, // # [b], each value in the vector is the length of the dipole step on its glide plane
            // periodicLoops
            const std::vector<int>& periodicLoopSlipSystemIDsIn,
            const std::vector<double>& periodicLoopRadiiIn,
            const std::vector<long int>& periodicLoopSidesIn,
            const Eigen::Matrix<double, Eigen::Dynamic,3>& periodicLoopCentersIn, // # each row in the matrix is a the center of the loop
            const double& periodicLoopTargetDensityIn,
            const long int& periodicLoopSegmentCountIn,
            const double& periodicLoopRadiusDistributionMeanIn,
            const double& periodicLoopRadiusDistributionStdIn,
            const double& periodicDipoleTargetDensityIn,
            const double& prismaticLoopTargetDensityIn,
            const std::map<int, double>& periodicLoopTargetDensitiesPerSlipSystemIn,
            const std::vector<int>& prismaticLoopSlipSystemIDsIn,
            const Eigen::Matrix<double, Eigen::Dynamic,3>& prismaticLoopCentersIn, // # each row in the matrix is a the center of the loop
            const std::vector<double>& prismaticLoopRadiiIn,
            const std::vector<double>& prismaticLoopStepsIn,
            const std::map<int, double>& prismaticLoopTargetDensitiesPerSlipSystemIn,
            const double& prismaticLoopRadiiMeanIn,
            const double& prismaticLoopRadiiStdIn,
            const double& prismaticLoopStepMeanIn,
            const double& prismaticLoopStepStdIn
            ) :
         microstructureType( microstructureTypeIn)
         ,style( styleIn)
         ,tag( tagIn)
         ,ddBase( ddBaseIn)
         ,periodicFaceIDs( ddBase.simulationParameters.periodicFaceIDs)
         ,periodicDipoleSlipSystemIDs( periodicDipoleSlipSystemIDsIn)
         ,periodicDipoleExitFaceIDs( periodicDipoleExitFaceIDsIn)
         ,periodicDipolePoints( periodicDipolePointsIn)
         ,periodicDipoleHeights( periodicDipoleHeightsIn)
         ,periodicDipoleNodes( periodicDipoleNodesIn)
         ,periodicDipoleGlideSteps( periodicDipoleGlideStepsIn)
         ,periodicLoopSlipSystemIDs( periodicLoopSlipSystemIDsIn)
         ,periodicLoopRadii( periodicLoopRadiiIn)
         ,periodicLoopSides( periodicLoopSidesIn)
         ,periodicLoopCenters( periodicLoopCentersIn)
         ,periodicLoopTargetDensity( periodicLoopTargetDensityIn)
         ,periodicLoopSegmentCount( periodicLoopSegmentCountIn)
         ,periodicLoopRadiusDistributionMean( periodicLoopRadiusDistributionMeanIn)
         ,periodicLoopRadiusDistributionStd( periodicLoopRadiusDistributionStdIn)
         ,periodicDipoleTargetDensity( periodicDipoleTargetDensityIn)
         ,prismaticLoopTargetDensity( prismaticLoopTargetDensityIn)
         ,periodicLoopTargetDensitiesPerSlipSystem( periodicLoopTargetDensitiesPerSlipSystemIn)
         ,prismaticLoopSlipSystemIDs( prismaticLoopSlipSystemIDsIn)
         ,prismaticLoopCenters( prismaticLoopCentersIn)
         ,prismaticLoopRadii( prismaticLoopRadiiIn)
         ,prismaticLoopSteps( prismaticLoopStepsIn)
         ,prismaticLoopTargetDensitiesPerSlipSystem( prismaticLoopTargetDensitiesPerSlipSystemIn)
         ,prismaticLoopRadiusDistributionMean( prismaticLoopRadiiMeanIn)
         ,prismaticLoopRadiusDistributionStd( prismaticLoopRadiiStdIn)
         ,prismaticLoopStepDistributionMean( prismaticLoopStepMeanIn)
         ,prismaticLoopStepDistributionStd( prismaticLoopStepStdIn)
      {}
   //PeriodicDipoleIndividualSpecification::PeriodicDipoleIndividualSpecification(
   //   const DislocationDynamicsBase<3>::DislocationDynamicsBaseType& ddBaseIn,
   //   const std::string& tagIn,
   //   const std::vector<int>& periodicDipoleSlipSystemIDsIn,
   //   const std::vector<int>& periodicDipoleExitFaceIDsIn, //  #  each value in the vector is the ID of the face from which the dipole enters/exits
   //   const Eigen::Matrix<double, Eigen::Dynamic,3>& periodicDipolePointsIn, // # each row in the matrix is a the "center" of the dipole
   //   const std::vector<double>& periodicDipoleHeightsIn, // # each value in the vector is the height of the dipole in number of slip planes
   //   const std::vector<int>& periodicDipoleNodesIn,
   //   const std::vector<double>& periodicDipoleGlideStepsIn // # [b], each value in the vector is the length of the dipole step on its glide plane
   //   ) :
   //   MicrostructureSpecification( "PeriodicDipole", "individual", tagIn, ddBaseIn)
   //   ,periodicFaceIDs( ddBase.periodicFaceIDs)
   //   ,periodicDipoleSlipSystemIDs( periodicDipoleSlipSystemIDsIn)
   //   ,periodicDipoleExitFaceIDs( periodicDipoleExitFaceIDsIn)
   //   ,periodicDipolePoints( periodicDipolePointsIn)
   //   ,periodicDipoleHeights( periodicDipoleHeightsIn)
   //   ,periodicDipoleNodes( periodicDipoleNodesIn)
   //   ,periodicDipoleGlideSteps( periodicDipoleGlideStepsIn)
   //{
   //}
}
#endif
