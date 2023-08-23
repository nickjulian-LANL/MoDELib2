/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_PrismaticLoopGeneratorInMem_cpp_
#define model_PrismaticLoopGeneratorInMem_cpp_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>
#include <random>
#include <iomanip>

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
#include <MicrostructureGeneratorInMem.h>
#include <PrismaticLoopGeneratorInMem.h>
#include <PlaneLineIntersection.h>

namespace model
{

    PrismaticLoopGeneratorInMem::PrismaticLoopGeneratorInMem(
          const std::shared_ptr<MicrostructureSpecification>& microSpec
          ) :
    /* init */ MicrostructureGeneratorBaseInMem( microSpec)
    {
       std::cout << "called constructor for PrismaticLoopGeneratorInMem" << std::endl;
    }

double PrismaticLoopGeneratorInMem::generateSingle(
      MicrostructureGeneratorInMem& mg,
      const int& rSS,
      const VectorDimD& guessCenter,
      const double& radius,
      const double& L
      )
{
    std::pair<bool,const Simplex<dim,dim>*> found(mg.ddBase.mesh.search(guessCenter));
    if(found.first)
    {
        const int grainID(found.second->region->regionID);
        assert(mg.ddBase.poly.grains.size()==1 && "Prismatic dislocations only supported for single crystals");
        const auto& grain(mg.ddBase.poly.grain(grainID));

        if(rSS>=0 && rSS<int(grain.singleCrystal->slipSystems().size()))
        {
            
            // sort slip systems in the zone axis with incrasing angle from a reference
            const RationalLatticeDirection<3>  b(grain.singleCrystal->slipSystems()[rSS]->s);
            const VectorDimD unitAxis(b.cartesian().normalized());
            const VectorDimD refNormal(Plane<3>::getL2G(unitAxis).col(0));
            std::map<double,ReciprocalLatticeDirection<3>> ssMap;
            for(const auto& ss : grain.singleCrystal->slipSystems())
            {
                if((ss->s-b).squaredNorm()<FLT_EPSILON)
                {
                    const double cosAngle(ss->unitNormal.dot(refNormal));
                    const double sinAngle(ss->unitNormal.dot(unitAxis.cross(refNormal)));
                    const double refAngle1(std::atan2(sinAngle,cosAngle));
                    const double refAngle2(std::atan2(-sinAngle,-cosAngle));
                    
                    if(std::fabs(refAngle1)<std::fabs(refAngle2))
                    {
                        ssMap.emplace(refAngle1,ss->n*(+1));
                    }
                    else
                    {
                        ssMap.emplace(refAngle2,ss->n*(-1));
                    }
                }
            }
            
            // Recenter and define sessile plane
            const VectorDimD center(grain.singleCrystal->snapToLattice(guessCenter).cartesian());
            const ReciprocalLatticeDirection<3> axis(grain.singleCrystal->reciprocalLatticeDirection(b.cartesian()));
            const int sessilePlaneHeight(axis.planeIndexOfPoint(center));
            const LatticePlaneKey sessilePlaneKey(axis,sessilePlaneHeight,grain.grainID);
            const auto sessilePlane(mg.ddBase.periodicGlidePlaneFactory.getFromKey(sessilePlaneKey));
            
            // Define prismatic planes.
            std::vector<LatticePlane> planes;
            
            for(const auto& pair : ssMap)
            {// first half of the prism
                const int planeHeight(pair.second.closestPlaneIndexOfPoint(center+pair.second.cartesian().normalized()*radius));
                planes.emplace_back(planeHeight,pair.second);
            }
            for(const auto& pair : ssMap)
            {// second half of the prism
                const int planeHeight(pair.second.closestPlaneIndexOfPoint(center-pair.second.cartesian().normalized()*radius));
                planes.emplace_back(planeHeight,pair.second);
            }
            
            // Intersect two consecutive prismatic planes and sessile plane to find polygon points
            std::vector<VectorDimD> sessileLoopNodePos;
            for(size_t k1=0;k1<planes.size();++k1)
            {
                size_t k2(k1==planes.size()-1? 0 : k1+1);
                const auto& plane1(planes[k1]);
                const auto& plane2(planes[k2]);
                
                Eigen::Matrix<double,3,3> N;
                N.col(0)=sessilePlane->referencePlane->unitNormal;
                N.col(1)=plane1.n.cartesian().normalized();
                N.col(2)=plane2.n.cartesian().normalized();
                
                Eigen::Matrix<double,3,3> P;
                P.col(0)=sessilePlane->referencePlane->P;
                P.col(1)=plane1.planeOrigin();
                P.col(2)=plane2.planeOrigin();
                
                PlanesIntersection<3> pi(N,P,FLT_EPSILON);
                const std::pair<bool,VectorDimD> intPair(pi.snap(VectorDimD::Zero()));
                if(intPair.first)
                {
                    sessileLoopNodePos.push_back(intPair.second);
                }
                else
                {
                    throw std::runtime_error("Cannot determine prismatic polygon.");
                }
                
            }
            
            double loopLength=0.0;
            for(size_t k1=0;k1<sessileLoopNodePos.size();++k1)
            {
                size_t k2(k1==sessileLoopNodePos.size()-1? 0 : k1+1);
                loopLength+=(sessileLoopNodePos[k2]-sessileLoopNodePos[k1]).norm();
            }
            
            
            mg.insertJunctionLoop(sessileLoopNodePos,sessilePlane,
                                  b.cartesian(),sessilePlane->referencePlane->unitNormal,
                                  center,grain.grainID,DislocationLoopIO<dim>::SESSILELOOP);
            
            
            //Create glissile loops on prism planes
            const auto periodicShifts(mg.ddBase.mesh.periodicBasis());
            Eigen::Matrix<double,3,3> box(Eigen::Matrix<double,3,3>::Zero());
            Eigen::Matrix<double,3,3> invBox(Eigen::Matrix<double,3,3>::Zero());
            
            if(periodicShifts.size()==3)
            {
                for(int k=0;k<3;++k)
                {
                    box.col(k)=periodicShifts[k];
                }
                invBox=box.inverse();
            }
            else
            {
                throw std::runtime_error("Cannot determine periodic box size.");
            }
            
            //const double L(50.0);
            const VectorDimD step(L*b.cartesian());
            for(size_t k1=0;k1<planes.size();++k1)
            {
                size_t k2(k1==planes.size()-1? 0 : k1+1);
                
                const Eigen::Matrix<double,3,1> boxCoord((invBox*(sessileLoopNodePos[k2]-mg.ddBase.mesh.xMin())).array().floor().matrix());
                const Eigen::Matrix<double,3,1> shift(box*boxCoord);
                
                if(mg.ddBase.mesh.searchRegion(grainID,sessileLoopNodePos[k2]-shift).first)
                {
                    std::vector<VectorDimD> loopNodePos;
                    loopNodePos.push_back(sessileLoopNodePos[k2]-shift);
                    loopNodePos.push_back(sessileLoopNodePos[k1]-shift);
                    loopNodePos.push_back(sessileLoopNodePos[k1]+step-shift);
                    loopNodePos.push_back(sessileLoopNodePos[k2]+step-shift);
                    
                    const int planeHeight(planes[k2].n.planeIndexOfPoint(loopNodePos.front()));
                    const LatticePlaneKey planeKey(planes[k2].n,planeHeight,grain.grainID);
                    const auto periodicGlidePlane(mg.ddBase.periodicGlidePlaneFactory.getFromKey(planeKey));
                    
                    mg.insertJunctionLoop(loopNodePos,periodicGlidePlane,
                                          b.cartesian(),periodicGlidePlane->referencePlane->unitNormal,
                                          periodicGlidePlane->referencePlane->P,grain.grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                }
                else
                {
                    throw std::runtime_error("Mesh does not contain shifted point.");
                }
            }
            return loopLength;
            
        }
        else
        {
            if(rSS<0)
            {
                std::cout<<"Skipping slip system "<<rSS<<std::endl;
            }
            else
            {
                throw std::runtime_error("slipSystem "+std::to_string(rSS)+" not found, skipping.");
            }
        }
    }
    else
    {
        std::cout<<"Center outside mesh, skipping prismatic loop"<<std::endl;
    }
    return 0.0;
}

void PrismaticLoopGeneratorInMem::generateDensity( MicrostructureGeneratorInMem& mg)
{
    std::cout<<magentaBoldColor<<"Generating prismatic loop density"<<defaultColor<<std::endl;
    const double targetDensity(
              microstructureSpecification->prismaticLoopTargetDensity);
    if(targetDensity>0.0)
    {
        const double radiusDistributionMean(
            microstructureSpecification->prismaticLoopRadiusDistributionMean);
        const double radiusDistributionStd(
            microstructureSpecification->prismaticLoopRadiusDistributionStd);
        std::normal_distribution<double> radiusDistribution(radiusDistributionMean/mg.ddBase.poly.b_SI,radiusDistributionStd/mg.ddBase.poly.b_SI);
        std::mt19937 generator;
        double density=0.0;
        while(density<targetDensity)
        {
            const std::pair<LatticeVector<dim>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
            const LatticeVector<dim> L0=rp.first;
            const size_t grainID=rp.second;
            std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID).singleCrystal->slipSystems().size()-1);
            const int rSS(ssDist(generator)); // a random SlipSystem
            const double radius(radiusDistribution(generator));
            try
            {
                
                density+=generateSingle(mg,rSS,L0.cartesian(),radius,50.0)/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                std::cout<<"prismatic loop density="<<density<<std::endl;
            }
            catch(const std::exception& e)
            {
                
            }
        }
    }
}

void PrismaticLoopGeneratorInMem::generateIndividual(MicrostructureGeneratorInMem& mg)
{
    const std::vector<int> prismaticLoopSlipSystemIDs(
              microstructureSpecification->prismaticLoopSlipSystemIDs
              );
    
    if(prismaticLoopSlipSystemIDs.size())
    {
        std::cout<<magentaBoldColor<<"Generating individual prismatic loops"<<defaultColor<<std::endl;
        const std::vector<double> prismaticLoopRadii(
               microstructureSpecification->prismaticLoopRadii
               );
        const Eigen::Matrix<double,Eigen::Dynamic,dim>
            prismaticLoopCenters(
                  microstructureSpecification->prismaticLoopCenters
                  );
        const std::vector<double> prismaticLoopSteps(
                  microstructureSpecification->prismaticLoopSteps
                  );
        
        if(prismaticLoopSlipSystemIDs.size()!=prismaticLoopRadii.size())
        {
            throw std::runtime_error("prismaticLoopSlipSystemIDs.size()="+std::to_string(prismaticLoopSlipSystemIDs.size())+" NOT EQUAL TO prismaticLoopRadii.size()="+std::to_string(prismaticLoopRadii.size()));
        }
        if(int(prismaticLoopSlipSystemIDs.size())!=prismaticLoopCenters.rows())
        {
            throw std::runtime_error("prismaticLoopSlipSystemIDs.size()="+std::to_string(prismaticLoopSlipSystemIDs.size())+" NOT EQUAL TO prismaticLoopCenters.rows()="+std::to_string(prismaticLoopCenters.rows()));
        }
        if(prismaticLoopSlipSystemIDs.size()!=prismaticLoopSteps.size())
        {
            throw std::runtime_error("prismaticLoopSlipSystemIDs.size()="+std::to_string(prismaticLoopSlipSystemIDs.size())+" NOT EQUAL TO prismaticLoopSteps.size()="+std::to_string(prismaticLoopSteps.size()));
        }
        
        for(size_t k=0;k<prismaticLoopSlipSystemIDs.size();++k)
        {
            generateSingle(mg,prismaticLoopSlipSystemIDs[k],prismaticLoopCenters.row(k),prismaticLoopRadii[k]/mg.ddBase.poly.b_SI,prismaticLoopSteps[k]);
        }
    }
    
}

void PrismaticLoopGeneratorInMem::generateDensitiesPerSlipSystem(
      MicrostructureGeneratorInMem& mg
      )
{
   std::cout<<magentaBoldColor<<"Generating prismatic loop densities per slip system"<<defaultColor<<std::endl;

   const std::map<int, double> ssDensities(
      microstructureSpecification->prismaticLoopTargetDensitiesPerSlipSystem);
   //the following instantiations could be made slipsystem specific
   const double radiusDistributionMean(
       microstructureSpecification->prismaticLoopRadiusDistributionMean);
   const double radiusDistributionStd(
       microstructureSpecification->prismaticLoopRadiusDistributionStd);
   const double stepDistributionMean(
       microstructureSpecification->prismaticLoopStepDistributionMean);
   const double stepDistributionStd(
       microstructureSpecification->prismaticLoopStepDistributionStd);
   std::normal_distribution<double> radiusDistribution(radiusDistributionMean/mg.ddBase.poly.b_SI,radiusDistributionStd/mg.ddBase.poly.b_SI);
   std::normal_distribution<double> stepDistribution(stepDistributionMean,stepDistributionStd);
   std::mt19937 generator;

   // iterate over ssDensities
   for ( const auto& ssDensitiesItr : ssDensities)
   {
      //std::cout << "ssDensitiesItr.second: " << ssDensitiesItr.second << ", ";
      size_t targetSlipSystemID( ssDensitiesItr.first);
      double targetDensity( ssDensitiesItr.second);
      if( targetDensity > 0.0)
      {
         double density=0.0;
         while( density < targetDensity)
         {
            const std::pair<LatticeVector<dim>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
            const LatticeVector<dim> L0=rp.first;
            const size_t grainID=rp.second;

            // check validity of targetSlipSystemID
            bool slipSystemIDGood = false;
            for ( const auto& ss : mg.ddBase.poly.grain(grainID).singleCrystal->slipSystems())
            {
              if ( targetSlipSystemID == ss->sID)
              {
                 slipSystemIDGood = true;
              }
            }

            if ( ! slipSystemIDGood)
            {
               std::cout << "warning: "
                  << "from PrismaticLoopGeneratorInMem::generateDensityPerSlipSystem()"
                  << " parameter microstructureSpecification->prismaticLoopTargetDensitiesPerSlipSystem "
                  << " target slip system ID "
                  << targetSlipSystemID << " is invalid" << std::endl;
               //continue; // find another point
               return; // fail 
            }

            std::cout << "adding a prismatic loop. density " << density << ", targetDensity " << targetDensity << std::endl; // debug
            try
            {
               const double radius(radiusDistribution(generator));
               const int steps( static_cast<int>(round(stepDistribution(generator))));
               std::cout << "radius " << radius << ", steps " << steps << std::endl; // debug
               density+=generateSingle(mg,targetSlipSystemID,L0.cartesian(),radius,steps)/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
               // Note: density assigned 0 when switching slip systems
               std::cout << "slip system " << targetSlipSystemID
                  << " prismatic loop density=" << density << std::endl;
            }
            catch(const std::exception& e)
            {
                    
            }
         }
      }
   }
}
} // namespace model
#endif
