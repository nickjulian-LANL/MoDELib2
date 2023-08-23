/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_PeriodicLoopGeneratorInMem_cpp_
#define model_PeriodicLoopGeneratorInMem_cpp_


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
#include <PeriodicLoopGeneratorInMem.h>
#include <PlaneLineIntersection.h>

namespace model
{

    PeriodicLoopGeneratorInMem::PeriodicLoopGeneratorInMem(
          //const PeriodicLoopIndividualSpecification& microSpec
          //const std::shared_ptr<PeriodicLoopIndividualSpecification>& microSpec
          const std::shared_ptr<MicrostructureSpecification>& microSpec
          //const std::shared_ptr<PeriodicLoopIndividualSpecification>& microSpec
          ) :
    /* init */ MicrostructureGeneratorBaseInMem( microSpec)
    {
       std::cout << "called constructor for PeriodicLoopGeneratorInMem" << std::endl;
    }


    void PeriodicLoopGeneratorInMem::generateIndividual(
          MicrostructureGeneratorInMem& mg
          )
    {
      std::cout << "microstructureSpecification->periodicLoopSlipSystemIDs.size(): "
         << microstructureSpecification->periodicLoopSlipSystemIDs.size() 
              << std::endl; // debug
      const std::vector<int> periodicLoopSlipSystemIDs(
              microstructureSpecification->periodicLoopSlipSystemIDs
              );
        
      if(periodicLoopSlipSystemIDs.size())
      {
         std::cout<<magentaBoldColor<<"Generating individual periodic loop"<<defaultColor<<std::endl;
         const std::vector<double> periodicLoopRadii(
               microstructureSpecification->periodicLoopRadii
               );
         const Eigen::Matrix<double,Eigen::Dynamic,dim>
            periodicLoopCenters(
                  microstructureSpecification->periodicLoopCenters
                  );
         const std::vector<long int> periodicLoopSides(
               microstructureSpecification->periodicLoopSides
               );
         
         if(periodicLoopSlipSystemIDs.size()!=periodicLoopRadii.size())
         {
             throw std::runtime_error("periodicLoopSlipSystemIDs.size()="+std::to_string(periodicLoopSlipSystemIDs.size())+" NOT EQUAL TO periodicLoopRadii.size()="+std::to_string(periodicLoopRadii.size()));
         }
         if(int(periodicLoopSlipSystemIDs.size())!=periodicLoopCenters.rows())
         {
             throw std::runtime_error("periodicLoopSlipSystemIDs.size()="+std::to_string(periodicLoopSlipSystemIDs.size())+" NOT EQUAL TO periodicLoopCenters.rows()="+std::to_string(periodicLoopCenters.rows()));
         }
         if(periodicLoopSlipSystemIDs.size()!=periodicLoopSides.size())
         {
             throw std::runtime_error("periodicLoopSlipSystemIDs.size()="+std::to_string(periodicLoopSlipSystemIDs.size())+" NOT EQUAL TO periodicLoopSides.size()="+std::to_string(periodicLoopSides.size()));
         }
         
         for(size_t k=0; k < periodicLoopSlipSystemIDs.size(); ++k)
         {
            generateSingle(
                  mg,
                  periodicLoopSlipSystemIDs[k],
                  periodicLoopCenters.row(k),
                  periodicLoopRadii[k]/mg.ddBase.poly.b_SI,
                  periodicLoopSides[k]
                  );
         }
      }
   }

   void PeriodicLoopGeneratorInMem::generateDensity(
         MicrostructureGeneratorInMem& mg
         )
   {
        std::cout<<magentaBoldColor<<"Generating periodic loop density"<<defaultColor<<std::endl;
        const double targetDensity(
              microstructureSpecification->periodicLoopTargetDensity);
        if(targetDensity>0.0)
        {
            const int numberOfSides(
              microstructureSpecification->periodicLoopSegmentCount);
            const double radiusDistributionMean(
                  microstructureSpecification->periodicLoopRadiusDistributionMean);
            const double radiusDistributionStd(
                  microstructureSpecification->periodicLoopRadiusDistributionStd);
            std::normal_distribution<double> radiusDistribution(radiusDistributionMean/mg.ddBase.poly.b_SI,radiusDistributionStd/mg.ddBase.poly.b_SI);
            std::mt19937 generator;
            double density=0.0;
            while( density < targetDensity)
            {
                const std::pair<LatticeVector<dim>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
                const LatticeVector<dim> L0=rp.first;
                const size_t grainID=rp.second;
                std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID).singleCrystal->slipSystems().size()-1);
                const int rSS(ssDist(generator)); // a random SlipSystem
                const double radius(radiusDistribution(generator));
                try
                {
                    generateSingle(mg,rSS,L0.cartesian(),radius,numberOfSides);
                    density+=2.0*std::numbers::pi*radius/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                    std::cout<<"periodic loop density="<<density<<std::endl;
                }
                catch(const std::exception& e)
                {
                    
                }
            }
        }
    }

   void PeriodicLoopGeneratorInMem::generateDensitiesPerSlipSystem(
         MicrostructureGeneratorInMem& mg
         )
   {
      std::cout<<magentaBoldColor<<"Generating periodic loop density"<<defaultColor<<std::endl;

      //const std::map<int, double> ssDensities(
      //   microstructureSpecification->periodicLoopTargetDensitiesPerSlipSystem);
      //the following five instantiations could be made slipsystem specific
      const int numberOfSides(
         microstructureSpecification->periodicLoopSegmentCount);
      const double radiusDistributionMean(
         microstructureSpecification->periodicLoopRadiusDistributionMean);
      const double radiusDistributionStd(
         microstructureSpecification->periodicLoopRadiusDistributionStd);
      std::normal_distribution<double> radiusDistribution(radiusDistributionMean/mg.ddBase.poly.b_SI,radiusDistributionStd/mg.ddBase.poly.b_SI);
      std::mt19937 generator;

      // iterate over ssDensities
      for ( const auto& ssDensitiesItr : microstructureSpecification->periodicLoopTargetDensitiesPerSlipSystem)
      {
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
                     << "from PeriodicLoopGeneratorInMem::generateDensityPerSlipSystem()"
                     << " parameter microstructureSpecification->periodicLoopTargetDensitiesPerSlipSystem "
                     << " target slip system ID "
                     << targetSlipSystemID << " is invalid" << std::endl;
                  //continue; // find another point
                  return; // fail 
               }

               const double radius(radiusDistribution(generator));
               try
               {
                  generateSingle(mg,targetSlipSystemID,L0.cartesian(),radius,numberOfSides);
                  // Note: density assigned 0 when switching slip systems
                  density+=2.0*std::numbers::pi*radius/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                  std::cout << "slip system " << targetSlipSystemID
                     << " periodic loop density=" << density << std::endl;
               }
               catch(const std::exception& e)
               {
                       
               }
            }
         }
      }
   }

//    void PeriodicDipoleGeneratorInMem::insertJunctionLoop(MicrostructureGeneratorInMem& mg,
//                                                     std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>>& uniqueNetworkNodeMap,
//                                                     const std::vector<VectorDimD>& loopNodePos,
//                                                     const std::shared_ptr<PeriodicGlidePlane<3>>& periodicPlane,
//                                                     const VectorDimD& b,
//                                                     const VectorDimD& unitNormal,
//                                                     const VectorDimD& P0,
//                                                     const size_t& grainID,
//                                                     const DislocationLoopIO<dim>::DislocationLoopType& loopType)
//    {
//        std::vector<PolyPoint> dummyPolyPoints;
//        std::vector<std::pair<VectorDimD, const PolyPoint *const>> loopNodePosTemp;
//        for(const auto& pos : loopNodePos)
//        {
//            dummyPolyPoints.push_back(PolyPoint());
//            loopNodePosTemp.emplace_back(pos, &dummyPolyPoints.back());
//        }
//
//        const auto ppi(periodicPlane->polygonPatchIntersection(loopNodePosTemp));
//        const size_t loopID(mg.insertLoop(b,unitNormal,P0,grainID,loopType));
//        std::vector<size_t> loopNodeIDs;
//        for(const auto &tup : ppi)
//        {
//            const VectorDimD loopNodePos(periodicPlane->referencePlane->globalPosition(std::get<0>(tup)));
//            const VectorDimD networkNodePos(loopNodePos+std::get<1>(tup));
//            const auto networkNodeIter(uniqueNetworkNodeMap.find(networkNodePos));
//            if(networkNodeIter==uniqueNetworkNodeMap.end())
//            {// no NetworkNode found at current position
//                uniqueNetworkNodeMap.emplace(networkNodePos,mg.insertNetworkNode(networkNodePos)); // insert NetworkNode and store its ID
//            }
//            loopNodeIDs.push_back(mg.insertLoopNode(loopID,loopNodePos,uniqueNetworkNodeMap.at(networkNodePos),std::get<1>(tup),std::get<2>(tup))); // insert LoopNode and store its ID
//        }
//        mg.insertLoopLinks(loopID,loopNodeIDs);
//    }


    void PeriodicLoopGeneratorInMem::generateSingle(
          MicrostructureGeneratorInMem& mg,
          const int& rSS,
          const VectorDimD& center,
          const double& radius,
          const int& sides
          )
    {
        std::pair<bool,const Simplex<dim,dim>*> found(mg.ddBase.mesh.search(center));
        if(!found.first)
        {
            std::cout<<"Point "<<center.transpose()<<" is outside mesh. EXITING."<<std::endl;
            exit(EXIT_FAILURE);
        }
        
        const int grainID(found.second->region->regionID);
        assert(mg.ddBase.poly.grains.size()==1 && "Periodic dislocations only supported for single crystals");
        const auto& grain(mg.ddBase.poly.grain(grainID));
        
        if(rSS>=0 && rSS<int(grain.singleCrystal->slipSystems().size()))
        {
            const auto& slipSystem(*grain.singleCrystal->slipSystems()[rSS]);
            
            
            const long int planeIndex(slipSystem.n.closestPlaneIndexOfPoint(center));
            GlidePlaneKey<3> glidePlaneKey(planeIndex, slipSystem.n);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
            const VectorDimD P0(glidePlane->referencePlane->snapToPlane(center));
            
            std::vector<VectorDimD> loopNodePos;
            for(int k=0;k< sides;++k)
            {
                loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*std::numbers::pi/sides,slipSystem.unitNormal)*slipSystem.s.cartesian().normalized()*radius);
            }
            
            
//            std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>> uniqueNetworkNodeMap; // networkNodePosition->networkNodeID

            mg.insertJunctionLoop(loopNodePos,glidePlane,
                               slipSystem.s.cartesian(),glidePlane->referencePlane->unitNormal,
                                  P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP);
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

}
#endif
