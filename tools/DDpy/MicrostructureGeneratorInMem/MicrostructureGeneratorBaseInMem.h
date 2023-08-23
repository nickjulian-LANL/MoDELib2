/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_MicrostructureGeneratorBaseInMem_H_
#define model_MicrostructureGeneratorBaseInMem_H_

#include <iostream>
#include <string>

#include <TerminalColors.h>

#include <string>
#include <Eigen/Dense>

#include <LatticeModule.h>
#include <GlidePlaneModule.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DislocationDynamicsBase.h>
#include <MicrostructureSpecification.h>
//#include <MicrostructureGeneratorInMem.h>
#include <MicrostructureGeneratorBase.h> // contains struct PolyPoint


namespace model
{
   class MicrostructureGeneratorInMem;

   //struct PolyPoint // use MicrostructureGeneratorBase::PolyPoint instead
   //{
   //   std::shared_ptr<PeriodicPlanePatch<3>> periodicPlanePatch() const;
   //};

   struct MicrostructureGeneratorBaseInMem
   {
      constexpr static int dim=3;
      typedef Eigen::Matrix<double,dim,1> VectorDimD;
      typedef Eigen::Matrix<long int,dim,1> VectorDimI;
      typedef LatticeDirection<dim> LatticeDirectionType;
      typedef DislocationLoopIO<dim>::DislocationLoopType DislocationLoopType;
      typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
      typedef Eigen::Matrix<long int,dim,dim> MatrixDimI;

      const std::shared_ptr<MicrostructureSpecification>& microstructureSpecification;
      const std::string type;
      const std::string style;
      const std::string tag;

      MicrostructureGeneratorBaseInMem(
         const std::shared_ptr<MicrostructureSpecification>& microSpec
      );

      virtual void generateIndividual( MicrostructureGeneratorInMem&) =0;
      virtual void generateDensity( MicrostructureGeneratorInMem&) =0;
      virtual void generateDensitiesPerSlipSystem( MicrostructureGeneratorInMem&) =0;
    };
} // namespace model
#endif
