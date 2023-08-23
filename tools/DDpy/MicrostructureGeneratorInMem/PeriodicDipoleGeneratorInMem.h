/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_PeriodicDipoleGeneratorInMem_H_
#define model_PeriodicDipoleGeneratorInMem_H_
#include <chrono>
#include <random>
#include <cmath>
#include <string>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGeneratorBaseInMem.h>
#include <MicrostructureGeneratorInMem.h>
#include <MicrostructureSpecification.h>
namespace model
{
   class PeriodicDipoleGeneratorInMem : public MicrostructureGeneratorBaseInMem
   {
      public:
      static void generateSingle(
            MicrostructureGeneratorInMem& mg,
            const int& rSS,
            const VectorDimD& dipolePoint,
            const int& exitFaceID,
            const int& dipoleHeight,
            const int& dipoleNodes,
            const double& glideStep
            );
      void generateIndividual( MicrostructureGeneratorInMem& mg) override;
      void generateDensity( MicrostructureGeneratorInMem& mg) override;
      void generateDensitiesPerSlipSystem( MicrostructureGeneratorInMem& mg) override;

      PeriodicDipoleGeneratorInMem(
          const std::shared_ptr<MicrostructureSpecification>& microSpec
          //const std::shared_ptr<PeriodicDipoleIndividualSpecification>& microSpec
          );
   }; // class PeriodicDipoleGeneratorInMem

} // namespace model
#endif
