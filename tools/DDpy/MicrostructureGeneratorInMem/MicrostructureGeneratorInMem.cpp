/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGeneratorInMem_cpp_
#define model_MicrostructureGeneratorInMem_cpp_

#include <fstream>
#include <filesystem>


#include <MicrostructureGeneratorInMem.h>
#include <PeriodicDipoleGeneratorInMem.h>
#include <PeriodicLoopGeneratorInMem.h>
#include <PrismaticLoopGeneratorInMem.h>
//#include <SphericalInclusionsGenerator.h>
//#include <PolyhedronInclusionsGenerator.h>
//#include <IrradiationDefectsGenerator.h>

namespace model
{

    /**********************************************************************/
    MicrostructureGeneratorInMem::MicrostructureGeneratorInMem(
         DislocationDynamicsBase<3>::DislocationDynamicsBaseType& ddBaseIn,
         const std::list<std::shared_ptr<MicrostructureSpecification>> microstructureSpecifications) :
    ///* init*/ traitsIO( ddBase.simulationParameters.traitsIO) // moved to DislocationDynamicsBase::simulationParameters.traitsIO
    /* init */ configIO( ddBaseIn.simulationParameters.traitsIO.evlFolder)
    /* init */,auxIO( ddBaseIn.simulationParameters.traitsIO.auxFolder)
    /* init */,ddBase( ddBaseIn)
    /* init */,outputBinary(TextFileParser( ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("outputBinary",true)) // TODO: remove // outputBinary might not exist yet
    //,outputBinary( false)
    ///* init */,periodicFaceIDs(TextFileParser( ddBase.simulationParameters.traitsIO.polyFile).template readSet<int>("periodicFaceIDs",true)) // TODO: remove
    /* init*/,minSize(0.1*std::min(ddBase.mesh.xMax(0)-ddBase.mesh.xMin(0),std::min(ddBase.mesh.xMax(1)-ddBase.mesh.xMin(1),ddBase.mesh.xMax(2)-ddBase.mesh.xMin(2))))
    /* init*/,maxSize(std::max(ddBase.mesh.xMax(0)-ddBase.mesh.xMin(0),std::max(ddBase.mesh.xMax(1)-ddBase.mesh.xMin(1),ddBase.mesh.xMax(2)-ddBase.mesh.xMin(2))))
    ///* init*/,poly( ddBase.poly)
    ///* init*/,glidePlaneFactory(poly)
    ///* init*/,periodicGlidePlaneFactory(poly, glidePlaneFactory)
    {
        
        std::cout<<greenBoldColor<<"Generating microstructure for "<< ddBase.simulationParameters.traitsIO.simulationFolder <<defaultColor<<std::endl;
        
    /* debug **************************************************************/
        std::cout << "microstructureSpecifications types:" << std::endl; // debug
        for ( const auto& spec : microstructureSpecifications) // debug
           std::cout << spec->microstructureType << std::endl; // debug
    /* debug **************************************************************/
        
        // Some sanity checks
        if(ddBase.mesh.volume()<FLT_EPSILON)
        {
            throw std::runtime_error("mesh "+ddBase.simulationParameters.traitsIO.meshFile+" is empty.");
        }
        
        for(const auto& spec : microstructureSpecifications)
        {
            const std::string microstructureType( spec->microstructureType);
            std::cout << "microstructureType: " << microstructureType << std::endl; // debug
            std::cout << "spec->periodicDipoleSlipSystemIDs.size(): "
               << spec->periodicDipoleSlipSystemIDs.size() << std::endl; // debug
            const std::string tag( spec->tag);
            bool success(false);
            if(microstructureType=="PeriodicDipole")
            {
                success=this->emplace(tag,
                      new PeriodicDipoleGeneratorInMem( spec)).second;
                if ( success) 
                   std::cout << "success" << std::endl; // debug
            }
            else if(microstructureType=="PeriodicLoop")
            {
                success=this->emplace(tag,
                      new PeriodicLoopGeneratorInMem( spec)).second;
            }
            else if(microstructureType=="PrismaticLoop")
            {
                success=this->emplace(tag,
                      new PrismaticLoopGeneratorInMem( spec)).second;
            }
            //else if(microstructureType=="SphericalInclusions")
            //{
            //    success=this->emplace(tag,new SphericalInclusionsGeneratorInMem( spec)).second;
            //}
            //else if(microstructureType=="PolyhedronInclusions")
            //{
            //    success=this->emplace(tag,new PolyhedronInclusionsGeneratorInMem( spec)).second;
            //}
            //else if(microstructureType=="Irradiation")
            //{
            //    success=this->emplace(tag,new IrradiationDefectsGeneratorInMem( spec)).second;
            //}
            //else if(microstructureType=="VTK")
            //{
            //    success=this->emplace(tag,new VTKGeneratorInMem( spec)).second;
            //}
            else
            {
               std::cout<<"unkown microstructure type "<<microstructureType<<std::endl;
            }
            if(!success)
            {
                throw std::runtime_error("Duplicate microstructure tag "+tag+".");
            }
        }
        
        for(auto& gen : *this)
        {
            if(gen.second->style=="individual")
            {
                gen.second->generateIndividual(*this);
            }
            else if(gen.second->style=="density")
            {
                gen.second->generateDensity(*this);
            }
            else if(gen.second->style=="densitiesPerSlipSystem")
            {
                gen.second->generateDensitiesPerSlipSystem(*this);
            }
            else
            {
                throw std::runtime_error("Uknown style for generator "+gen.second->tag);
            }
        }
        
        // Call individual generators
        //            addStraightDislocations();
        //            addFrankReadSources();
        //            addSingleArmDislocations();
        //            addPrismaticLoops();
        //            addIndividualStraightDislocations();
        //            addFrankLoops();
        //            addNonPlanarLoops();
        //            // addPeriodicLoops();
        //            addStatisticallyHomegeneousPeriodicLoops();
        //    addIndividualPeriodicStraighDislocation();
        //            addStatisticallyHomegeneousPlanarDipolarLoops();
        //            addPeriodicJunctionLoops();
        //            addIrradiationLoops();
        //            addStackingFaultTetrahedra();
        //            addEshelbyInclusions();
        writeConfigFiles(0); // write the configuration to evl/evl_0.txt
        
    }

    const DDtraitsIO& MicrostructureGeneratorInMem::traits() const
    {
        return ddBase.simulationParameters.traitsIO;
    }


    const DDconfigIO<3>& MicrostructureGeneratorInMem::config() const
    {
        return configIO;
    }

    const DDauxIO<3>& MicrostructureGeneratorInMem::aux() const
    {
        return auxIO;
    }

    void MicrostructureGeneratorInMem::insertJunctionLoop(const std::vector<VectorDimD>& loopNodePos,
                                                 const std::shared_ptr<PeriodicGlidePlane<3>>& periodicPlane,
                                                 const VectorDimD& b,
                                                 const VectorDimD& unitNormal,
                                                 const VectorDimD& P0,
                                                 const size_t& grainID,
                                                 const DislocationLoopIO<3>::DislocationLoopType& loopType)
{
  std::vector<PolyPoint> dummyPolyPoints;
  std::vector<std::pair<VectorDimD, const PolyPoint *const>> loopNodePosTemp;
  for(const auto& pos : loopNodePos)
  {
      dummyPolyPoints.push_back(PolyPoint());
      loopNodePosTemp.emplace_back(pos, &dummyPolyPoints.back());
  }
  
  const auto ppi(periodicPlane->polygonPatchIntersection(loopNodePosTemp,true));
  const size_t loopID(insertLoop(b,unitNormal,P0,grainID,loopType));
  std::vector<size_t> loopNodeIDs;
  for(const auto &tup : ppi)
  {
      const VectorDimD loopNodePos(periodicPlane->referencePlane->globalPosition(std::get<0>(tup)));
      const VectorDimD networkNodePos(loopNodePos+std::get<1>(tup));
      const auto networkNodeIter(uniqueNetworkNodeMap.find(networkNodePos));
      if(networkNodeIter==uniqueNetworkNodeMap.end())
      {// no NetworkNode found at current position
          uniqueNetworkNodeMap.emplace(networkNodePos,insertNetworkNode(networkNodePos)); // insert NetworkNode and store its ID
      }
      loopNodeIDs.push_back(insertLoopNode(loopID,loopNodePos,uniqueNetworkNodeMap.at(networkNodePos),std::get<1>(tup),std::get<2>(tup))); // insert LoopNode and store its ID
  }
  insertLoopLinks(loopID,loopNodeIDs);
}

    size_t MicrostructureGeneratorInMem::insertLoop(const VectorDimD& b,const VectorDimD& unitNormal,const VectorDimD& P0,const size_t& grainID,const DislocationLoopType& loopType)
    {
        const size_t loopID(configIO.loops().size());
        configIO.loops().emplace_back(loopID, b,unitNormal,P0,grainID,loopType);
        return loopID;
    }

    size_t MicrostructureGeneratorInMem::insertLoopNode(const size_t& loopID,const VectorDimD& loopNodePos,const size_t& networkNodeID,const VectorDimD& loopNodeShift,const std::pair<short int,short int>& periodicEdgeIDs)
    {
        const size_t loopNodeID(configIO.loopNodes().size());
        configIO.loopNodes().emplace_back(loopNodeID,loopID,loopNodePos,networkNodeID,loopNodeShift,periodicEdgeIDs);
        return loopNodeID;
    }

    std::vector<size_t> MicrostructureGeneratorInMem::insertLoopLinks(const size_t& loopID,const std::vector<size_t>& loopNodeIDs)
    {
        std::vector<size_t> temp;
        for(size_t k=0;k<loopNodeIDs.size();++k)
        {
            const size_t k1=(k+1)<loopNodeIDs.size()? k+1 : 0;
            temp.push_back(configIO.loopLinks().size());
            const size_t sourceNodeID(loopNodeIDs[k ]);
            const size_t   sinkNodeID(loopNodeIDs[k1]);
            const auto sourceNode(configIO.loopNodes()[sourceNodeID]);
            const auto   sinkNode(configIO.loopNodes()[sinkNodeID]);
            configIO.loopLinks().emplace_back(loopID,sourceNodeID,sinkNodeID,(sourceNode.P-sinkNode.P).norm()>FLT_EPSILON,0);
        }
        return temp;
    }


    size_t MicrostructureGeneratorInMem::insertNetworkNode(const VectorDimD& networkNodePos)
    {
        const size_t networkNodeID(configIO.nodes().size());
        configIO.nodes().emplace_back(networkNodeID,networkNodePos,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
        return networkNodeID;
    }

    size_t MicrostructureGeneratorInMem::insertInclusion(const VectorDimD& pos,const double& R, const Eigen::Matrix<double,dim,dim>& eT, const double& vrc,const int&type)
    {
        const size_t inclusionID(configIO.sphericalInclusions().size()+configIO.polyhedronInclusions().size());
        configIO.sphericalInclusions().emplace_back(inclusionID,pos,R,eT,vrc,type);
        return inclusionID;
    }

size_t MicrostructureGeneratorInMem::insertInclusion(const std::map<size_t,Eigen::Vector3d>& polyNodes,const std::map<size_t,std::vector<size_t>>& faceMap, const Eigen::Matrix<double,dim,dim>& eT, const double& vrc,const int&type)
{
    const size_t inclusionID(configIO.sphericalInclusions().size()+configIO.polyhedronInclusions().size());
    configIO.polyhedronInclusions().emplace_back(inclusionID,eT,vrc,type);
    const size_t startNodeID(configIO.polyhedronInclusionNodes().size());
    size_t nodeCounter(0);
    for(const auto& node : polyNodes)
    {
        configIO.polyhedronInclusionNodes().emplace_back(startNodeID+nodeCounter,node.second);
        nodeCounter++;
    }

    for(const auto& pair : faceMap)
    {
        const size_t& faceID(pair.first);
        for(size_t k=0;k<pair.second.size();++k)
        {
            const size_t k1(k<pair.second.size()-1? k+1 : 0);
            const auto sourceIter(polyNodes.find(pair.second[k]));
            const auto sinkIter(polyNodes.find(pair.second[k1]));
            if(sourceIter!=polyNodes.end() && sinkIter!=polyNodes.end())
            {                
                const size_t sourceID(startNodeID +std::distance(polyNodes.begin(),sourceIter));
                const size_t sinkID(startNodeID +std::distance(polyNodes.begin(),sinkIter));
                configIO.polyhedronInclusionEdges().emplace_back(inclusionID,faceID,sourceID,sinkID);
            }
        }
    }
    
    return inclusionID;
}



    /**********************************************************************/
    void MicrostructureGeneratorInMem::writeConfigFiles(const size_t& fileID)
    {
        
        //const int outputGlidePlanes(TextFileParser(traitsIO.ddFile).readScalar<int>("outputGlidePlanes",true));
        //DC->DN->outputGlidePlanes // Note: DC and DN are not existing yet

        //if(outputGlidePlanes)
        //{
        //    for(const auto& loop : configIO.loops())
        //    {
        //        GlidePlaneKey<3> loopPlaneKey(loop.P, poly.grain(loop.grainID).singleCrystal->reciprocalLatticeDirection(loop.N));
        //        auxIO.glidePlanes().emplace_back(loopPlaneKey);
        //    }
        //}
                
        if(outputBinary)
        {
            std::cout<<greenBoldColor<<"Writing configuration to "<<configIO.getBinFilename(fileID)<<defaultColor<<std::endl;
            configIO.writeBin(fileID);
            auxIO.writeBin(fileID);
        }
        else
        {
            std::cout<<greenBoldColor<<"Writing configuration to "<<configIO.getTxtFilename(fileID)<<defaultColor<<std::endl;
            configIO.writeTxt(fileID); // writes to {configIO.folderName}/evl_{fileID}.txt
            auxIO.writeTxt(fileID); // writes to {auxIO.folderName}/ddAux_{fileID}.txt
        }
    }

bool MicrostructureGeneratorInMem::allPointsInGrain(const std::vector<VectorDimD>& points,const int& grainID)
{
    bool temp=true;
    for(const auto& point : points)
    {
        temp *= ddBase.mesh.searchRegion(grainID,point).first;
        if(!temp)
        {
            break;
        }
    }
    return temp;
}

}
#endif
