/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_MicrostructureSpecification_H_
#define model_MicrostructureSpecification_H_

#include <string>
#include <set>
//#include <MicrostructureGeneratorBaseInMem.h>
#include <DislocationDynamicsBase.h>

namespace model
{
   class MicrostructureSpecification
   {
      public:
      //constexpr static int dim=MicrostructureGeneratorBase::dim;
      std::string microstructureType; // PeriodicDipole, PeriodicLoop, Inclusions, Irradiation, VTK
      std::string style; // individual, density
      std::string tag; // e.g.: incd0, incl0, pdd0, pdl0, plpd0, plp0
      const DislocationDynamicsBase<3>& ddBase;
      // the following might be put into their own class, if inheritance works
      const std::set<int> periodicFaceIDs; // inherit from ddBase
      const std::vector<int> periodicDipoleSlipSystemIDs;
      const std::vector<int> periodicDipoleExitFaceIDs; //  #  each value in the vector is the ID of the face from which the dipole enters/exits
      const Eigen::Matrix<double, Eigen::Dynamic,3> periodicDipolePoints; // # each row in the matrix is a the "center" of the dipole 
      const std::vector<double> periodicDipoleHeights; // # each value in the vector is the height of the dipole in number of slip planes
      const std::vector<int> periodicDipoleNodes;
      const std::vector<double> periodicDipoleGlideSteps;
      const std::vector<int> periodicLoopSlipSystemIDs;
      const std::vector<double> periodicLoopRadii;
      const std::vector<long int> periodicLoopSides; // count of segments
      const Eigen::Matrix<double, Eigen::Dynamic,3> periodicLoopCenters; // # each row in the matrix is a the center of the loop

      const double periodicLoopTargetDensity;
      const long int periodicLoopSegmentCount;//single value for density creation
      const double periodicLoopRadiusDistributionMean;
      const double periodicLoopRadiusDistributionStd;

      const double periodicDipoleTargetDensity;
      const double prismaticLoopTargetDensity;
      const std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
      const std::vector<int> prismaticLoopSlipSystemIDs;
      const Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
      const std::vector<double> prismaticLoopRadii;
      const std::vector<double> prismaticLoopSteps;
      const std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
      const double prismaticLoopRadiusDistributionMean;
      const double prismaticLoopRadiusDistributionStd;
      const double prismaticLoopStepDistributionMean;
      const double prismaticLoopStepDistributionStd;

      MicrostructureSpecification(
            const std::string& microstructureTypeIn,
            const std::string& styleIn,
            const std::string& tagIn,
            const DislocationDynamicsBase<3>::DislocationDynamicsBaseType& ddBaseIn,

            const std::vector<int>& periodicDipoleSlipSystemIDsIn,
            const std::vector<int>& periodicDipoleExitFaceIDsIn, //  #  each value in the vector is the ID of the face from which the dipole enters/exits
            const Eigen::Matrix<double, Eigen::Dynamic,3>& periodicDipolePointsIn, // # each row in the matrix is a the "center" of the dipole 
            const std::vector<double>& periodicDipoleHeightsIn, // # each value in the vector is the height of the dipole in number of slip planes
            const std::vector<int>& periodicDipoleNodesIn,
            const std::vector<double>& periodicDipoleGlideStepsIn, // # [b], each value in the vector is the length of the dipole step on its glide plane
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
            const std::vector<int>& prismaticLoopSlipSystemIDs,
            const Eigen::Matrix<double, Eigen::Dynamic,3>& prismaticLoopCentersIn, // # each row in the matrix is a the center of the loop
            const std::vector<double>& prismaticLoopRadiiIn,
            const std::vector<double>& prismaticLoopStepsIn,
            const std::map<int, double>& prismaticLoopTargetDensitiesPerSlipSystemIn,
            const double& prismaticLoopRadiiMeanIn,
            const double& prismaticLoopRadiiStdIn,
            const double& prismaticLoopStepMeanIn,
            const double& prismaticLoopStepStdIn
            );// :
      //   microstructureType( microstructureTypeIn)
      //   ,style( styleIn)
      //   ,tag( tagIn)
      //   ,ddBase( ddBaseIn)
      //   ,periodicDipoleSlipSystemIDs( periodicDipoleSlipSystemIDsIn)
      //   ,periodicDipoleExitFaceIDs( periodicDipoleExitFaceIDsIn)
      //   ,periodicDipolePoints( periodicDipolePointsIn)
      //   ,periodicDipoleHeights( periodicDipoleHeightsIn)
      //   ,periodicDipoleNodes( periodicDipoleNodesIn) //   ,periodicDipoleGlideSteps( periodicDipoleGlideStepsIn)
      //{} 
   }; // class MicrostructureSpecification

   //// TODO: consider creating classes for each type of defect specification. I tried onece, but PeriodicDipoleGeneratorInMem() would then be passed a generic MicrostructureSpecification without knowing its a PeriodicDipoleIndividualSpecification.
   //class PeriodicDipoleIndividualSpecification : public MicrostructureSpecification
   //{
   //   public:
   //   const std::set<int>& periodicFaceIDs; // inherit from ddBase
   //   const std::vector<int>& periodicDipoleSlipSystemIDs;
   //   const std::vector<int>& periodicDipoleExitFaceIDs; //  #  each value in the vector is the ID of the face from which the dipole enters/exits
   //   const Eigen::Matrix<double, Eigen::Dynamic,3>& periodicDipolePoints; // # each row in the matrix is a the "center" of the dipole 
   //   const std::vector<double>& periodicDipoleHeights; // # each value in the vector is the height of the dipole in number of slip planes
   //   const std::vector<int>& periodicDipoleNodes;
   //   const std::vector<double>& periodicDipoleGlideSteps;


   //   PeriodicDipoleIndividualSpecification(
   //      const DislocationDynamicsBase<3>::DislocationDynamicsBaseType& ddBase,
   //      const std::string& tag,

   //      const std::vector<int>& periodicDipoleSlipSystemIDs,
   //      const std::vector<int>& periodicDipoleExitFaceIDs, //  #  each value in the vector is the ID of the face from which the dipole enters/exits
   //      const Eigen::Matrix<double, Eigen::Dynamic,3>& periodicDipolePoints, // # each row in the matrix is a the "center" of the dipole 
   //      const std::vector<double>& periodicDipoleHeights, // # each value in the vector is the height of the dipole in number of slip planes
   //   const std::vector<int>& periodicDipoleNodesIn,
   //   const std::vector<double>& periodicDipoleGlideStepsIn // # [b], each value in the vector is the length of the dipole step on its glide plane
   //      );
   //}; // class PeriodicDipoleIndividualSpecification 
}
#endif
