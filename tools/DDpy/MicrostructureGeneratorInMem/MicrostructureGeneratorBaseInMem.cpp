/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_MicrostructureGeneratorBaseInMem_cpp_
#define model_MicrostructureGeneratorBaseInMem_cpp_

#include <iostream>
#include <string>

#include <TerminalColors.h>


#include <MicrostructureGeneratorBaseInMem.h>
#include <MicrostructureGeneratorBase.h>

namespace model
{

    //std::shared_ptr<PeriodicPlanePatch<3>> PolyPoint::periodicPlanePatch() const
    //{
    //    return nullptr;
    //}

    MicrostructureGeneratorBaseInMem::MicrostructureGeneratorBaseInMem(
          const std::shared_ptr<MicrostructureSpecification>&
            microstructureSpecificationIn
          //const DislocationDynamicsBase<3>& ddBase
          ) :
    /* init */ microstructureSpecification( microstructureSpecificationIn)
    /* init */,type( microstructureSpecification->microstructureType)
    /* init */,style( microstructureSpecification->style)
    /* init */,tag( microstructureSpecification->tag)
   {
       std::cout << "called constructor for MicrostructureGeneratorBaseInMem" << std::endl;
        std::cout<<greenColor<<tag<<" -> "<<type<<" "<<style<<defaultColor<<std::endl; // debug
   }
}
#endif
