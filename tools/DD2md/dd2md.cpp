#ifndef dd2md_readLmp_CPP_
#define dd2md_readLmp_CPP_

#include <dd2md.h>
#include <LmpReader.h>

/**********************************************************************/


//model::AtomDisplacementGenerator::fieldPointsOutputType getFieldPoints()
//{
//   return fieldPointsOutput;
//}

void model::AtomDisplacementGenerator::readBurgersMagnitude(
      const std::string& materialPath)
{
   burgersMagnitude = model::TextFileParser( materialPath).readScalar<double>("b_SI",true);
   std::cout << "burgersMagnitude: " << burgersMagnitude << std::endl; // debug
   return;
}

void model::AtomDisplacementGenerator::readLammpsConfigurationFile(
      const std::string& lammpsDataFilePath
      )
{
   LmpReader lammpsReader(
          lammpsDataFilePath,
          //1.0,
          1e-10/burgersMagnitude, // b_SI:2.556e-10 #scaleFactor for atom positions
          lmpDeformationMatrix,
          debugFlag);
   if ( lammpsReader.readLmpStream(
            lammpsBoxBounds,
            masses,
            atomIDs,
            atomTypes,
            fieldPoints, // atomPositions
            pointIDs
            ) != EXIT_SUCCESS)
   {
      std::cout << "error: readLammpsConfigurationFile() failed to read "
         << lammpsDataFilePath
         << std::endl;
      //if ( periodicFlag)
      //   DC->DN->simulationParameters.simulationType = 2; // PERIODIC_IMAGES
      return;
   }
   // TODO: shift atoms so that modelib slip systems and lammps slip
   //  systems coincide
   shift_atoms();

   return;
}

void model::AtomDisplacementGenerator::resetStaticIDs()
{
   model::StaticID<model::EshelbyInclusionBase<3>>::set_count(0);
   model::StaticID<model::Lattice<3>>::set_count(0);
   model::StaticID<model::DislocationMobilityBase>::set_count(0);
   model::StaticID<model::PeriodicPlaneNode<3>>::set_count(0);
   model::StaticID<model::PeriodicPlanePatch<3>>::set_count(0);
   model::StaticID<model::PeriodicPatchBoundary<3>>::set_count(0);
   //model::StaticID<model::Derived>::set_count(0);
   model::StaticID<model::SlipSystem>::set_count(0);
   //model::StaticID<model::LoopPathClipperNode>::set_count(0);
   model::StaticID<model::PlanarMeshFace<3>>::set_count(0);
   //model::StaticID<model::LoopPathClipperNode>::set_count(0);
   model::StaticID<model::NetworkNode<model::DislocationNode<3,0>>>::set_count(0);
   model::StaticID<model::NetworkLink<model::DislocationSegment<3,0>>>::set_count(0);
   model::StaticID<model::Loop<model::DislocationLoop<3,0>>>::set_count(0);
   model::StaticID<model::LoopNode<model::DislocationLoopNode<3,0>>>::set_count(0);
   model::StaticID<model::MeshPlane<3>>::set_count(0);
   model::StaticID<model::SimplexBase<3,0>>::set_count(0);
   model::StaticID<model::SimplexBase<3,1>>::set_count(0);
   model::StaticID<model::SimplexBase<3,2>>::set_count(0);
   model::StaticID<model::SimplexBase<3,3>>::set_count(0);
   //model::StaticID<model::SequentialOutputFile<prefix,auto>>::set_count(0);
   return;
}

void model::AtomDisplacementGenerator::readddBase()
{
   if ( ddBase != nullptr) resetStaticIDs(); // TODO: is this necessary?
   ddBase =  std::unique_ptr<DislocationDynamicsBaseType>(
         new DislocationDynamicsBaseType( dddPath)
         );
   configFields = std::unique_ptr<DDconfigFields<3>>(
         new DDconfigFields<3>( *ddBase, configIO));
   return;
}

void model::AtomDisplacementGenerator::readDefectiveCrystal()
{
   if ( ddBase == nullptr) readddBase();
   DC = std::unique_ptr<DefectiveCrystalType>(
         new DefectiveCrystalType( *ddBase)
         );
   return;
}

void model::AtomDisplacementGenerator::regenerateMicrostructure()
{
   if ( ddBase == nullptr) readddBase();
   // instantiate a MicrostructureGenerator
   model::MicrostructureGenerator mg( *ddBase);
   mg.readMicrostructureFile();
   mg.writeConfigFiles(0);
   //microstructureGenerator.readMicrostructureFile();
   //microstructureGenerator.writeConfigFiles(0);
   if ( debugFlag)
      std::cout << "finished call to MicrostructureGenerator" << std::endl;

   //setCurrentStep(0);// DefectiveCrystal will use runID to read some things

   // instantiate a DefectiveCrystalType DC( ddBase)
   DC = std::unique_ptr<DefectiveCrystalType>(
         new DefectiveCrystalType( *ddBase)
         );
   if ( debugFlag)
      std::cout << "finished call to DefectiveCrystalType" << std::endl;

   return;
}

void model::AtomDisplacementGenerator::readLammpsBounds(
      const std::string& lammpsFilePath
      )
{
   LmpReader lammpsReader(
          lammpsFilePath,
          //1.0,
          1e-10/burgersMagnitude, // b_SI:2.556e-10 #scaleFactor for atom positions
          Eigen::Matrix<double,3,3>::Identity(),
          debugFlag);
   if ( lammpsReader.readLmpStreamBounds( lammpsBoxBounds) != EXIT_SUCCESS)
   {
      std::cout << "error: cannot read boundaries from lammps data file."
         << std::endl;
   }
   return;
}

//void model::AtomDisplacementGenerator::applyDisplacement( VectorDim& x)
//{
//   VectorDim temp( VectorDim::Zero());
//   //double temp( 0.0);
//   for ( const auto& patch : loopPatches())
//   {
//      const auto& loop( configIO.loop( patch.first));
//      for ( const auto& shift : periodicShifts)
//      {
//         temp += patch.second.solidAngle( x + shift)
//                  /4.0/std::numbers::pi*loop.B;
//      }
//   }
//   return;
//}

const model::DDconfigIO<3>& model::AtomDisplacementGenerator::config(
      ) const
{
   return configIO;
}

model::DDconfigIO<3>& model::AtomDisplacementGenerator::config()
{
   return configIO;
}

typename model::AtomDisplacementGenerator::VectorDim
model::AtomDisplacementGenerator::dislocationPlasticDisplacement(
      const VectorDim& x) const
{
    return configFields->dislocationPlasticDisplacement( x);
}

typename model::AtomDisplacementGenerator::VectorDim
model::AtomDisplacementGenerator::dislocationPlasticDisplacement(
      const double& x,
      const double& y,
      const double& z) const
{
    return dislocationPlasticDisplacement((VectorDim()<<x,y,z).finished());
}

void model::AtomDisplacementGenerator::readConfiguration(
      const size_t& runID)
{
   //std::cout << "configIO.getTxtFilename(0): " << configIO.getTxtFilename(0) << std::endl; // debug
   configIO.read( runID);
   configFields->updateConfiguration();
   return;
}

double model::AtomDisplacementGenerator::solidAngle(
      const VectorDim& x) const
{
   return configFields->solidAngle( x);
}

double model::AtomDisplacementGenerator::solidAngle(
            const double& x,
            const double& y,
            const double& z
      ) const
{
   return solidAngle((VectorDim()<<x,y,z).finished());
}

void model::AtomDisplacementGenerator::computeDisplacements()
      //const std::string& lammpsFilePath)
{
   if ( DC == nullptr)
   {
      std::cout << "error: DefectiveCrystal not yet instantiated"
         << std::endl; 
      return;
   }
   if ( DC->DN == nullptr)
   {
      std::cout << "error: DefectiveCrystal->DislocationNetwork "
         << "not yet instantiated"
         << std::endl; 
      return;
   }

   if ( fieldPoints.size() <= 0)
   {
      std::cout << "error: computDisplacements() called while fieldPoints"
         << " is empty" << std::endl;
      return;
   }

   if ( debugFlag) std::cout << "Computing DislocationDisplacement at field points..." << std::endl;
   //DC->DN->displacement( fieldPoints);
   if (displacements.size() != fieldPoints.size())
   {
      displacements.resize( fieldPoints.size());
   }
   for (unsigned int k=0; k < fieldPoints.size(); ++k)
   {
      displacements[k] = DC->DN->displacement( fieldPoints[k]);
   }

   //if ( debugFlag) std::cout << "Applying external elastic field ..." << std::endl;
   //applyExternalElasticField();  // updates elasticDisplacements member
   if ( debugFlag) std::cout << "Finished computeDisplacements ..." << std::endl;

   //// return simulationParameters to periodic if it was altered
   //if ( periodicFlag)
   //   DC->DN->simulationParameters.simulationType = 2; // PERIODIC_IMAGES
   return;
} // void model::AtomDisplacementGenerator::computeDisplacements()

void model::AtomDisplacementGenerator::computeDisplacements2()
     // const std::string& lammpsFilePath)
{

   if ( debugFlag) std::cout << "Computing dislocationPlasticDisplacement() at field points..." << std::endl;
   // Iterate over fieldPoints and modify by applying
   //  dislocationPlasticDisplacement( ); to each.
   for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   {
      //VectorDim temp( VectorDim::Zero());
      //temp += dislocationPlasticDisplacement( (burgersMagnitude/(1e-10))*fieldPoints[ ii]);
      //temp += dislocationPlasticDisplacement(fieldPoints[ ii]);
      fieldPoints[ ii] = dislocationPlasticDisplacement(
           fieldPoints[ ii]
           );
      //fieldPoints[ ii] = ((1e-10)/burgersMagnitude)*(fieldPoints[ ii] + temp);
   }
   //DC->DN->displacement( fieldPoints);

   //if ( debugFlag) std::cout << "Applying external elastic field ..." << std::endl;
   //applyExternalElasticField();  // updates elasticDisplacements member
   if ( debugFlag) std::cout << "Finished computeDisplacements2 ..." << std::endl;

   return;
} // void model::AtomDisplacementGenerator::computeDisplacements2()

void model::AtomDisplacementGenerator::applyDisplacements()
      //const std::string& lammpsFilePath)
{
   if ( DC == nullptr)
   {
      std::cout << "error: DefectiveCrystal not yet instantiated"
         << std::endl; 
      return;
   }
   if ( DC->DN == nullptr)
   {
      std::cout << "error: DefectiveCrystal->DislocationNetwork "
         << "not yet instantiated"
         << std::endl; 
      return;
   }

   if ( fieldPoints.size() <= 0)
   {
      std::cout << "error: applyDisplacements() called while fieldPoints"
         << " is empty" << std::endl;
      return;
   }

   if (displacements.size() != fieldPoints.size())
   {
      std::cout << "error: applyDisplacements() called while displacements.size() " << displacements.size() << " != fieldPoints.size()" << fieldPoints.size() << std::endl;
      return;
   }
   if ( debugFlag) std::cout << "Applying displacements to fieldPoints..." << std::endl;
   for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   {
      fieldPoints[ii] += displacements[ii];
   }
   return;
}

void model::AtomDisplacementGenerator::applyExternalElasticField()
{
   // Elastic displacement container 
   if ( debugFlag) std::cout<<"applyExternalElasticField(): Computing the displacement associated with the external stressZX...\n"<<std::flush;
   const double gamma_ZX = 8.e8/DC->DN->poly.mu_SI;
   VectorDim box_size( ddBase->mesh.xMax() - ddBase->mesh.xMin());

   for (unsigned int k=0; k < fieldPoints.size(); ++k)
   {
      const double dis = gamma_ZX * (
                              fieldPoints[ k] - ddBase->mesh.xMin()
                           ).dot( VectorDim::UnitZ());
      elasticDisplacements.push_back( dis * (VectorDim::UnitX()));
   }
   if ( debugFlag) std::cout<<"applyExternalElasticField(): finished computing the displacement associated with the external stressZX...\n"<<std::flush;
   return;
}


//typedef std::tuple<size_t, > patchGlidePlaneType;

std::tuple<
   pybind11::array_t<double, pybind11::array::c_style>, // planeNormals
   pybind11::array_t<double, pybind11::array::c_style>, // Burger's vectors
   std::map< int, // size_t and long int were disallowed during compilation
      std::map< int,
         pybind11::array_t<double, pybind11::array::c_style>
         >
      > // polygons
   >
model::AtomDisplacementGenerator::getPatchGlidePlanes()
{
   // NOTES:
   // std::map<std::shared_ptr<PeriodicPlanePatch<_dim>>,std::vector<Eigen::Matrix<double,_dim-1,1>>> _patches;
   // patch.first points to a PeriodicPlanePatch<3>
   //  patch.first->patchBoundary->referencePlane points to a GlidePlane<3>
   //  which is a LatticePlane and MeshPlane<3>
   //  LatticePlane inherits from LatticePlaneBase 
   //   which inherits from ReciprocalLatticeDirection<3>
   //   which inherits from Eigen::Matrix<long int, dim, 1>
   //  MeshPlane<3> inherits from Plane<3> which contains VectorDim unitNormal
   // patch.second is a vector of 2x1 matrices ( the boundary points)
   //std::map<size_t, VectorDim> ; // one per loop
   pybind11::array_t<double, pybind11::array::c_style> planeNormals; // per loop
   pybind11::array_t<double, pybind11::array::c_style> burgersVectors; // per loop
   //pybind11::array_t<double, pybind11::array::c_style> polygons; // (#loops, #patches per loop, #vertices,3)
   std::map< int,  // one element per loop
      std::map< int, // one element per patch
         pybind11::array_t<double, pybind11::array::c_style> // vertices 3-D
      >
   > polygons;

   int loopCount; //loopCount = DC->DN->loops().size();
   std::vector< int > patchCounts; // one patch count per loop
   std::vector<std::vector< int>> vertexCounts; // one count per patch

   // count the dimensions required for the output arrays
   int loopNumber; loopNumber = 0;
   for( const auto& loop : DC->DN->loops())
   {
      if( loop.second.lock()->glidePlane)
      {
         if(loop.second.lock()->getSlippedArea() > FLT_EPSILON)
         {// a right-handed normal for the loop can be determined
            //patchCounts.emplace_back( loop.second.lock()->_patches.size());
            int patchNumber; patchNumber = 0;
            vertexCounts.emplace_back( std::vector<int>());
            // inspect all _patches 
            for ( const auto& patch : loop.second.lock()->_patches.globalPatches())
            {
               ++patchNumber;
               vertexCounts.back().emplace_back(
                     patch.second.size() // number of vertices on the patch
                     );
            }
            ++loopNumber;
            patchCounts.emplace_back( patchNumber);
         }
      }
   }
   loopCount = loopNumber;

   planeNormals.resize( {loopCount, 3});
   burgersVectors.resize( {loopCount, 3});

   pybind11::buffer_info planeNormalsBuf = planeNormals.request();
   //std::cout << "planeNormalsBuf.ndim: " << planeNormalsBuf.ndim  << std::endl;
   //std::cout << "planeNormalsBuf.shape[0]: " << planeNormalsBuf.shape[0] << std::endl;
   //std::cout << "planeNormalsBuf.shape[1]: " << planeNormalsBuf.shape[0] << std::endl;

   pybind11::buffer_info burgersVectorsBuf = burgersVectors.request();
   //std::cout << "burgersVectorsBuf.ndim: " << burgersVectorsBuf.ndim  << std::endl;
   //std::cout << "burgersVectorsBuf.shape[0]: " << burgersVectorsBuf.shape[0] << std::endl;
   //std::cout << "burgersVectorsBuf.shape[1]: " << burgersVectorsBuf.shape[0] << std::endl;


   double* planeNormalsPtr
      = static_cast<double*>( planeNormalsBuf.ptr);
   double* burgersVectorsPtr
      = static_cast<double*>( burgersVectorsBuf.ptr);

   // inspect all loops
   loopNumber = 0;
   for( const auto& loop : DC->DN->loops())
   {
      if( loop.second.lock()->glidePlane)
      {
         if(loop.second.lock()->getSlippedArea() > FLT_EPSILON)
         {
            planeNormalsPtr[ 0 + 3*loopNumber]
                  = loop.second.lock()->rightHandedUnitNormal()[0];
            planeNormalsPtr[ 1 + 3*loopNumber]
                  = loop.second.lock()->rightHandedUnitNormal()[1];
            planeNormalsPtr[ 2 + 3*loopNumber]
                  = loop.second.lock()->rightHandedUnitNormal()[2];

            burgersVectorsPtr[ 0 + 3*loopNumber]
               = loop.second.lock()->burgers()[0] * burgersMagnitude;
            burgersVectorsPtr[ 1 + 3*loopNumber]
               = loop.second.lock()->burgers()[1] * burgersMagnitude;
            burgersVectorsPtr[ 2 + 3*loopNumber]
               = loop.second.lock()->burgers()[2] * burgersMagnitude;

            polygons[ loopNumber]
               = std::map< int,
                  pybind11::array_t<double, pybind11::array::c_style> >();

            // inspect all _patches 
            int patchNumber; patchNumber = 0;
            // loop.second.lock()->_patches.localPatches() // 2-D positions
            // loop.second.lock()->_patches.globalPatches() // 3-D positions
            for ( const auto& patch : loop.second.lock()->_patches.globalPatches())
            {
               polygons[ loopNumber][ patchNumber]
                  =  pybind11::array_t<double, pybind11::array::c_style>();
               polygons[ loopNumber][ patchNumber].resize(
                     {vertexCounts[ loopNumber][ patchNumber], 3} );
               pybind11::buffer_info polygonBuf
                  = polygons[ loopNumber][ patchNumber].request();
               //std::cout << "loop " << loopNumber << ", patch " << patchNumber << ", polygonBuf.ndim: " << polygonBuf.ndim  << std::endl;
               //std::cout << "loop " << loopNumber << ", patch " << patchNumber << ", polygonBuf.shape[0]: " << polygonBuf.shape[0] << std::endl;
               //std::cout << "loop " << loopNumber << ", patch " << patchNumber << ", polygonBuf.shape[1]: " << polygonBuf.shape[0] << std::endl;

               double* polygonPtr = static_cast<double*>( polygonBuf.ptr);

               // coordinates of the polygon:
               const auto patchGlidePlane(
                     patch.first->patchBoundary->referencePlane
                     );
               // patch.second is a vector of 3x1 matrices ( the boundary points)
               for (int vv=0;
                     vv < vertexCounts[ loopNumber][ patchNumber]; ++vv)
               {
                  VectorDim tempVec( patch.second[ vv]);
                  polygonPtr[ 0 + 3*vv] = tempVec[0] * burgersMagnitude/1e-10;
                  polygonPtr[ 1 + 3*vv] = tempVec[1] * burgersMagnitude/1e-10;
                  polygonPtr[ 2 + 3*vv] = tempVec[2] * burgersMagnitude/1e-10;
               }
               ++patchNumber;
            }
            ++loopNumber;
         }
      }
   }
   return std::tuple<
      pybind11::array_t<double, pybind11::array::c_style>, // planeNormals
      pybind11::array_t<double, pybind11::array::c_style>, // burgersVectors
      std::map< int,
         std::map< int,
            pybind11::array_t<double, pybind11::array::c_style>
            >
         > // polygons
      >
      ( planeNormals, // one per loop
        burgersVectors, // one per loop
        polygons // indexed as: [loop#][patch#][vertex#]
      );
}

//std::tuple<
//   std::map<size_t, model::AtomDisplacementGenerator::VectorDim>,
//   std::map<size_t, model::AtomDisplacementGenerator::VectorDim>,
//   std::map< std::pair<size_t,size_t>,
//      std::vector< model::AtomDisplacementGenerator::VectorDim>>>
//model::AtomDisplacementGenerator::getPatchGlidePlanes()
//{
//   // NOTES:
//   // std::map<std::shared_ptr<PeriodicPlanePatch<_dim>>,std::vector<Eigen::Matrix<double,_dim-1,1>>> _patches;
//   // patch.first points to a PeriodicPlanePatch<3>
//   //  patch.first->patchBoundary->referencePlane points to a GlidePlane<3>
//   //  which is a LatticePlane and MeshPlane<3>
//   //  LatticePlane inherits from LatticePlaneBase 
//   //   which inherits from ReciprocalLatticeDirection<3>
//   //   which inherits from Eigen::Matrix<long int, dim, 1>
//   //  MeshPlane<3> inherits from Plane<3> which contains VectorDim unitNormal
//   // patch.second is a vector of 2x1 matrices ( the boundary points)
//   std::map<size_t, VectorDim> planeNormals; // one per loop
//   std::map<size_t, VectorDim> burgersVectors; // one per loop
//   std::map<std::pair<size_t,size_t>, std::vector<VectorDim>> polygons;// several per loop
//   // inspect all loops
//   size_t loopNumber; loopNumber = 0;
//   for( const auto& loop : DC->DN->loops())
//   {
//      if( loop.second.lock()->glidePlane)
//      {
//         if(loop.second.lock()->_slippedArea > FLT_EPSILON)
//         {
//            planeNormals[ loopNumber]
//                  = loop.second.lock()->rightHandedUnitNormal();
//            //planeNormals.emplace_back(
//            //   patch.first->patchBoundary->referencePlane->unitNormal
//            //   );
//
//            burgersVectors[ loopNumber] = loop.second.lock()->burgers();
//
//            // inspect all _patches 
//            for ( const auto& patch : loop.second.lock()->_patches)
//            {
//               size_t patchNumber; patchNumber = 0;
//               // coordinates of the polygon:
//               const auto patchGlidePlane(
//                     patch.first->patchBoundary->referencePlane
//                     );
//               // patch.second is a vector of 2x1 matrices ( the boundary points)
//               std::vector<VectorDim> polygonPoints;
//               for ( size_t ii=0; ii < patch.second.size(); ++ii)
//               {
//                  polygonPoints.emplace_back(
//                        patchGlidePlane->globalPosition( patch.second[ ii]) 
//                        );
//               }
//               polygons[
//                     std::pair<size_t,size_t>( loopNumber, patchNumber)
//                  ] = polygonPoints;
//               ++patchNumber;
//            }
//            ++loopNumber;
//         }
//      }
//   }
//   return std::make_tuple<
//   std::map<size_t, VectorDim>,
//   std::map<size_t, VectorDim>,
//   std::map< std::pair<size_t,size_t>, std::vector<VectorDim>>
//      >( planeNormals, burgersVectors, polygons);
//
//
//   //for(const auto& pair : _patches)
//   //{
//   //  const auto patchGlidePlane(pair.first->patchBoundary->referencePlane);
//   //
//   //// pair.first->patchBoundary->referencePlane is a:
//   ////  std::shared_ptr<GlidePlane<dim>>
//   //
//   //  std::vector<std::pair<VectorDim,VectorDim>> segments;
//   //  for(size_t k=0;k<pair.second.size();++k)
//   //  {
//   //     const size_t k1(k+1==pair.second.size()? 0 : k+1);
//   //     segments.emplace_back(
//   //        patchGlidePlane->globalPosition(pair.second[k]),
//   //        patchGlidePlane->globalPosition(pair.second[k1])
//   //        );
//   //     // coordinates of the polygon:
//   //     //  (patchGlidePlane->globalPosition(pair.second[k]),
//   //     //   patchGlidePlane->globalPosition(pair.second[k1])
//   //  }
//   //  temp += planarSolidAngle(x,patchGlidePlane->P,rightHandedUnitNormal(),segments);
//
//   return;
//}

double model::AtomDisplacementGenerator::getBurgersMagnitude()
{
   return DC->DN->poly.b_SI/1e-10;
}

void model::AtomDisplacementGenerator::writeConfigurationToFile(
      const std::string& outputFilePath)
{
   if ( debugFlag)
      std::cout << "Outputing configuration to file " << outputFilePath
               << " ..." << std::flush;
   const auto t2= std::chrono::system_clock::now();

   // open output file and truncate its existing contents 
   std::ofstream outputFile( outputFilePath,
              std::ofstream::out | std::ofstream::trunc);

   if ( fieldPoints.size() <= 0)
   {
      std::cout << "error "
         << "AtomDisplacementGenerator::writeConfigurationToFile"
         << " called while fieldPoints is empty" << std::endl;
      return;
   }
   if ( fieldPoints.size() != atomTypes.size()) 
   {
      std::cout << "error "
         << "AtomDisplacementGenerator::writeConfigurationToFile"
         << " fieldPoints.size() != atomType.size()) " << std::endl;
      return;
   }

   outputFile << "# LAMMPS data file written by DD2MD" << std::endl << std::endl;
   outputFile << fieldPoints.size() << " atoms" << std::endl << std::endl;
   outputFile << masses.size() << " atom types" << std::endl << std::endl;
   outputFile
      // << std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[0] << " "
      //<< std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[1]
      << " xlo xhi" << std::endl
      //<< std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[2] << " "
      //<< std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[3]
      << " ylo yhi" << std::endl
      //<< std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[4] << " "
      //<< std::setw(16)
      << std::setprecision(12)  << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[5]
      << " zlo zhi" << std::endl << std::endl;
   outputFile << "Masses" << std::endl << std::endl;
   for ( const auto& mm : masses)
   {
      outputFile << mm.first << " " << mm.second << std::endl;
   }

   outputFile << std::endl << "Atoms # atomic" << std::endl << std::endl;
   for (unsigned int k=0; k < fieldPoints.size(); ++k)
   {
       //outputFile << fieldPoints[k].pointID << " "
       outputFile << pointIDs[k] << " "
                    //<< std::setw(8)
                    << atomTypes[ k] << " "
                    << std::setiosflags( std::ios::fixed );
       //if(applyStress)
       //outputFile //<< std::setw(16)
       //   << std::setprecision(12)
       //   << fieldPoints[k].transpose() * DC->DN->poly.b_SI/1e-10
       //      + fieldPoints[k].transpose() * DC->DN->poly.b_SI/1e-10
       //      + elasticDisplacements[k].transpose() * DC->DN->poly.b_SI/1e-10
       //   << "\n";
       // if elasticDisplacements shouldn't be applied, then the following line is desired, not the above
       //VectorDim tempPoint( lmpDeformationMatrixInverse *fieldPoints[k]);
       outputFile //<< std::setw(15)
        << std::setprecision(12)
        << (fieldPoints[k]).transpose() * DC->DN->poly.b_SI /1e-10
          //+ (fieldPoints[k]).transpose() * DC->DN->poly.b_SI /1e-10
        //<< (lmpDeformationMatrixInverse * fieldPoints[k]).transpose() * DC->DN->poly.b_SI /1e-10
        //  + (lmpDeformationMatrixInverse * fieldPoints[k]).transpose() * DC->DN->poly.b_SI /1e-10
        << "\n";
   }
   if ( debugFlag)
      std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<std::endl;

   return;
}

void model::AtomDisplacementGenerator::writeDisplacementsToFile(
      const std::string& outputFilePath)
{
   if ( debugFlag)
      std::cout << "Outputting displacements to file " << outputFilePath
               << " ..." << std::flush;
   const auto t2= std::chrono::system_clock::now();

   // open output file and truncate its existing contents 
   std::ofstream outputFile( outputFilePath,
              std::ofstream::out | std::ofstream::trunc);

   if ( displacements.size() <= 0)
   {
      std::cout << "error "
         << "AtomDisplacementGenerator::writeDisplacementsToFile"
         << " called while displacements is empty" << std::endl;
      return;
   }
   if ( displacements.size() != atomTypes.size()) 
   {
      std::cout << "error "
         << "AtomDisplacementGenerator::writeDisplacementsToFile"
         << " displacements.size() != atomType.size()) " << std::endl;
      return;
   }

   outputFile << "# LAMMPS data file written by DD2MD" << std::endl << std::endl;
   outputFile << displacements.size() << " atoms" << std::endl << std::endl;
   outputFile << masses.size() << " atom types" << std::endl << std::endl;
   outputFile
      // << std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[0] << " "
      //<< std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[1]
      << " xlo xhi" << std::endl
      //<< std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[2] << " "
      //<< std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[3]
      << " ylo yhi" << std::endl
      //<< std::setw(16)
      << std::setprecision(12) << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[4] << " "
      //<< std::setw(16)
      << std::setprecision(12)  << (burgersMagnitude/1e-10) * lammpsDeformedBoxBounds[5]
      << " zlo zhi" << std::endl << std::endl;
   outputFile << "Masses" << std::endl << std::endl;
   for ( const auto& mm : masses)
   {
      outputFile << mm.first << " " << mm.second << std::endl;
   }

   outputFile << std::endl << "Atoms # atomic" << std::endl << std::endl;
   for (unsigned int k=0; k < displacements.size(); ++k)
   {
       //outputFile << displacements[k].pointID << " "
       outputFile << pointIDs[k] << " "
                    //<< std::setw(8)
                    << atomTypes[ k] << " "
                    << std::setiosflags( std::ios::fixed );
       //if(applyStress)
       //outputFile //<< std::setw(16)
       //   << std::setprecision(12)
       //   << displacements[k].transpose() * DC->DN->poly.b_SI/1e-10
       //      + displacements[k].transpose() * DC->DN->poly.b_SI/1e-10
       //      + elasticDisplacements[k].transpose() * DC->DN->poly.b_SI/1e-10
       //   << "\n";
       // if elasticDisplacements shouldn't be applied, then the following line is desired, not the above
       //VectorDim tempPoint( lmpDeformationMatrixInverse *displacements[k]);
       outputFile //<< std::setw(15)
        << std::setprecision(12)
        << (displacements[k]).transpose() * DC->DN->poly.b_SI /1e-10
          //+ (displacements[k]).transpose() * DC->DN->poly.b_SI /1e-10
        //<< (lmpDeformationMatrixInverse * displacements[k]).transpose() * DC->DN->poly.b_SI /1e-10
        //  + (lmpDeformationMatrixInverse * displacements[k]).transpose() * DC->DN->poly.b_SI /1e-10
        << "\n";
   }
   if ( debugFlag)
      std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<std::endl;

   return;
}


void model::AtomDisplacementGenerator::setCurrentStep( const long int& step)
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot setCurrentStep until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }

   ddBase->simulationParameters.runID = step;
   DC->externalLoadController->update( DC->plasticStrain());
   ddBase->simulationParameters.manageRestart();
   return;
}


void model::AtomDisplacementGenerator::shift_atoms(
      //fieldPointsType& fieldPoints
      )
{
   // iterate over all positions to find the lowest
   if ( fieldPoints.size() <= 0)
   {
      std::cout << "AtomDisplacementGenerator::shift_atoms(): "
         << "passed an empty set of atoms ..."
         << std::endl;
      return;
   }
   //double oneOverRootThree( 1.0/sqrt(3.0));
   //VectorDim unitVector( VectorDim()<<oneOverRootThree);

   // TODO: parallelize each of the following loops
   // find the lowest extents of atom positions
   VectorDim lowestExtent( fieldPoints[0]);
   for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   {
      if ( fieldPoints[ii](0) < lowestExtent(0))
         lowestExtent(0) = fieldPoints[ii](0);
      if ( fieldPoints[ii](1) < lowestExtent(1))
         lowestExtent(1) = fieldPoints[ii](1);
      if ( fieldPoints[ii](2) < lowestExtent(2))
         lowestExtent(2) = fieldPoints[ii](2);
   }
   // translate atoms so that their lowest extent is (0,0,0)
   for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   {

      VectorDim temp;
      temp << 
         fieldPoints[ii](0) - lowestExtent(0),
         fieldPoints[ii](1) - lowestExtent(1),
         fieldPoints[ii](2) - lowestExtent(2);
      fieldPoints[ii] = temp;
   }
   // translate atoms to ensure lammps and modelib slip planes coincide
   VectorDim displacementVector;
   double latticeConstant;
   if (! lattice.compare("bcc"))
   {
      //latticeConstant = (burgersMagnitude/1e-10) * 2.0/sqrt(3.0); 
      latticeConstant = 2.0/sqrt(3.0); 
      displacementVector << 0.5, 0.25, 0.125;
   }
   else if (! lattice.compare("fcc"))
   {
      //latticeConstant = (burgersMagnitude/1e-10) * sqrt(2.0); 
      latticeConstant = sqrt(2.0); 
      displacementVector << 0.5, 0.0, 0.0;
   }
   displacementVector *= latticeConstant;
   displacementVector = (C2G.inverse())*displacementVector; 
   std::cout << "displacementVector: \n" << displacementVector << std::endl; // debug
   // TODO: incorporate atomic crystallographic orientation rather than use this
   for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   {
      fieldPoints[ii] += displacementVector;
   }
   // translate atoms whose positions exceed the periodic boundaries
   for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   {
      while ( fieldPoints[ii](0) < lammpsDeformedBoxBounds[0]) // xlo
         fieldPoints[ii](0) += lammpsDeformedBoxDimensions[0];
      while ( fieldPoints[ii](0) >= lammpsDeformedBoxBounds[1]) // xhi
         fieldPoints[ii](0) -= lammpsDeformedBoxDimensions[0];

      while ( fieldPoints[ii](1) < lammpsDeformedBoxBounds[2]) // ylo
         fieldPoints[ii](1) += lammpsDeformedBoxDimensions[1];
      while ( fieldPoints[ii](1) >= lammpsDeformedBoxBounds[3]) // yhi
         fieldPoints[ii](1) -= lammpsDeformedBoxDimensions[1];

      while ( fieldPoints[ii](2) < lammpsDeformedBoxBounds[4]) // zlo
         fieldPoints[ii](2) += lammpsDeformedBoxDimensions[2];
      while ( fieldPoints[ii](2) >= lammpsDeformedBoxBounds[5]) // zhi
         fieldPoints[ii](2) -= lammpsDeformedBoxDimensions[2];
   }
   return;
}

void model::AtomDisplacementGenerator::regeneratePolycrystalFile(
            const pybind11::array_t<double,
               pybind11::array::c_style | pybind11::array::forcecast>
               c2gIn,
            //const std::string& lammpsDataFilePath,
            const std::string& latticeIn,
            const std::string& material,
            const std::string& meshFilePath
            )

{
   auto c2gNp = c2gIn.unchecked<2>(); // for reading C2G input np.array
   if (! ((c2gNp.shape(0) == 3) && (c2gNp.shape(1) == 3)))
   {
      std::cout << "error: regeneratePolycrystalFile( C2G) requires C2G "
         << "to be a 3x3 numpy array consisting of rows of normalized "
         << "basis vectors of the crystal lattice." << std::endl;
      return;
   }
   bool latticeIsAcceptable = false;
   bool materialIsAcceptable = false;
   for ( const auto& lat : acceptableLattices)
   {
      if ( latticeIn == std::string(lat)) latticeIsAcceptable = true;
   }
   for ( const auto& mat : acceptableMaterials)
   {
      if ( material == std::string(mat)) materialIsAcceptable = true;
   }
   if ( ! ( latticeIsAcceptable && materialIsAcceptable))
   {
      std::cout << "error: unacceptable lattice or material: "
         << lattice << ", " << material << std::endl;
      std::cout << "  acceptable lattice or materials are: ";
      for ( const auto& mat : acceptableMaterials) std::cout << mat << ", ";
      for ( const auto& lat : acceptableLattices) std::cout << lat << ", ";
      return;
   }

   lattice = latticeIn; // assign member of AtomDisplacementGenerator

   std::cout << "lattice: " << lattice << ", material: " << material << std::endl; // debug

   C2G << c2gNp(0,0), c2gNp(0,1), c2gNp(0,2),
         c2gNp(1,0), c2gNp(1,1), c2gNp(1,2),
         c2gNp(2,0), c2gNp(2,1), c2gNp(2,2);
   std::cout << "C2G:\n" << C2G << std::endl; // debug
   std::string outputFilePath = dddPath + "/inputFiles/polycrystal.txt";

   // TODO: detect and move any existing polycrystal.txt file
   std::ofstream outputFile( outputFilePath,
              std::ofstream::out | std::ofstream::trunc);

   std::cout << "lammps boundaries: " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug

   double deltaX, deltaY, deltaZ; 
   deltaX = abs( lammpsBoxBounds[1] - lammpsBoxBounds[0]);
   deltaY = abs( lammpsBoxBounds[3] - lammpsBoxBounds[2]);
   deltaZ = abs( lammpsBoxBounds[5] - lammpsBoxBounds[4]);
   std::cout << "lammps: deltaX " << deltaX << ", deltaY " << deltaY
      << ", deltaZ " << deltaZ << std::endl; // debug
   
   MatrixDim AA; // AA scales the mesh: y=A(x-x0), where x is the input mesh
   if (! lattice.compare("bcc")) // .compare() returns 0 if they're equal
   {
      AA << -1.0, 1.0, 1.0,
         1.0, -1.0, 1.0,
         1.0, 1.0, -1.0;
      AA /= sqrt(3.0);
   }
   else if (! lattice.compare("fcc"))
   {
      AA << 0.0, 1.0, 1.0,
          1.0, 0.0, 1.0,
          1.0, 1.0, 0.0;
      AA /= sqrt(2.0);
   }
   else
   {
      std::cout << "error: lattice type not recognized while trying to "
         << "create matrix A" << std::endl;
      return;
   }
    
   AA = C2G * AA;

   MatrixDim f12, f31, f23;
   f12 << 1.0, skew, 0.0,
       0.0, 1.0, 0.0,
       0.0, 0.0, 1.0;
   f31 << 1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       skew, 0.0, 1.0;
   f23 << 1.0, 0.0, 0.0,
       0.0, 1.0, skew,
       0.0, 0.0, 1.0;
   MatrixDim deformingMatrix;
   deformingMatrix = f12 * (f23 * f31);

   MatrixDim scalingMatrix;
   scalingMatrix << deltaX, 0, 0,
      0, deltaY, 0,
      0, 0, deltaZ;

   deformingMatrix = ( AA.inverse()) * (deformingMatrix * scalingMatrix);
   // Create a transformation to be applied to the atoms to align
   //  their strained atomic planes to the perfect slip systems of DDD.
   // Atomic planes will still need to be shifted to make atomic slip
   //  planes coincide with DDD slip planes.
   MatrixDim deformingMatrixRounded( deformingMatrix.array().round());

   // deform AA to be written to polycrystal.txt
   AA = AA * deformingMatrixRounded;
   //std::cout  // debug
   //     << std::setw(22) << std::setprecision(15) << AA(0,0)
   //     << std::setw(22) << std::setprecision(15) << AA(0,1)
   //     << std::setw(22) << std::setprecision(15) << AA(0,2)
   //     << std::endl
   //     << std::setw(22) << std::setprecision(15) << AA(1,0)
   //     << std::setw(22) << std::setprecision(15) << AA(1,1)
   //     << std::setw(22) << std::setprecision(15) << AA(1,2)
   //     << std::endl
   //     << std::setw(22) << std::setprecision(15) << AA(2,0)
   //     << std::setw(22) << std::setprecision(15) << AA(2,1)
   //     << std::setw(22) << std::setprecision(15) << AA(2,2)
   //     << std::endl; // debug

   lmpDeformationMatrix = deformingMatrixRounded * (deformingMatrix.inverse()) ;
   lmpDeformationMatrixInverse =  lmpDeformationMatrix.inverse();
   std::cout << "lmpDeformationMatrix :\n" << lmpDeformationMatrix // debug
      << std::endl; // debug
   std::cout << "lmpDeformationMatrixInverse :\n" // debug
      << lmpDeformationMatrixInverse // debug
      << std::endl; // debug

   ////////////////////////////////////////////////
   /// use the local lmpDeformationMatrix to calculate lammpsDeformedBox
   // NOTE: assume skew == 0
   VectorDim lmpAxis00(
         lammpsBoxBounds[0], // xlo
         0.0,
         0.0);
   VectorDim lmpAxis01(
         lammpsBoxBounds[1], // xhi
         0.0,
         0.0);
   VectorDim lmpAxis10(
         0.0,
         lammpsBoxBounds[2], // ylo
         0.0);
   VectorDim lmpAxis11(
         0.0,
         lammpsBoxBounds[3], // yhi
         0.0);
   VectorDim lmpAxis20(
         0.0,
         0.0,
         lammpsBoxBounds[4]); // zlo
   VectorDim lmpAxis21(
         0.0,
         0.0,
         lammpsBoxBounds[5]); // zhi

   // transform lammps box boundary coordinates
   lammpsDeformedBoxBounds.clear();
   lammpsDeformedBoxBounds.resize(6);
   lammpsDeformedBoxBounds[0] // xlo
      = (lmpDeformationMatrix * lmpAxis00)[0];
   lammpsDeformedBoxBounds[1] // xhi
      = (lmpDeformationMatrix * lmpAxis01)[0];

   lammpsDeformedBoxBounds[2] // ylo
      = (lmpDeformationMatrix * lmpAxis10)[1];
   lammpsDeformedBoxBounds[3] // yhi
      = (lmpDeformationMatrix * lmpAxis11)[1];

   lammpsDeformedBoxBounds[4] // zlo
      = (lmpDeformationMatrix * lmpAxis20)[2];
   lammpsDeformedBoxBounds[5] // zhi
      = (lmpDeformationMatrix * lmpAxis21)[2];


   lammpsDeformedBoxDimensions.clear();
   lammpsDeformedBoxDimensions.resize(3);
   lammpsDeformedBoxDimensions[0] = lammpsDeformedBoxBounds[1]
                                     - lammpsDeformedBoxBounds[0];
   lammpsDeformedBoxDimensions[1] = lammpsDeformedBoxBounds[3]
                                     - lammpsDeformedBoxBounds[2];
   lammpsDeformedBoxDimensions[2] = lammpsDeformedBoxBounds[5]
                                     - lammpsDeformedBoxBounds[4];
   std::cout << "lammpsDeformedBoxBounds: (" 
      << lammpsDeformedBoxBounds[0] << ", "
      << lammpsDeformedBoxBounds[1] << ", "
      << lammpsDeformedBoxBounds[2] << ", "
      << lammpsDeformedBoxBounds[3] << ", "
      << lammpsDeformedBoxBounds[4] << ", "
      << lammpsDeformedBoxBounds[5] << ") " << std::endl; // debug
   std::cout << "lammpsDeformedBoxDimensions: (" 
      << lammpsDeformedBoxDimensions[0] << ", "
      << lammpsDeformedBoxDimensions[1] << ", "
      << lammpsDeformedBoxDimensions[2] << ")" << std::endl; // debug
   //////////////////////////////////////

   double x0x, x0y, x0z;
   x0x = lammpsBoxBounds[0]/deltaX; // lammpsBoxBounds[0] is xlo, deltaX = xhi -xlo
   x0y = lammpsBoxBounds[2]/deltaY; // lammpsBoxBounds[2] is ylo
   x0z = lammpsBoxBounds[4]/deltaZ; // lammpsBoxBounds[4] is zlo
   VectorDim x0;
   x0 << x0x, x0y, x0z; // scaling of the mesh is y=A(x-x0)
   std::cout << "x0: " << x0x << ", " << x0y << ", " << x0z << std::endl; // debug

   outputFile << "materialFile=" << material + ".txt;" << std::endl;
   outputFile << "enablePartials=0;"
      << std::endl;
   outputFile << "absoluteTemperature = 300; # [K] simulation temperature"
      << std::endl;
   outputFile << "meshFile=" << meshFilePath << ";"
      << std::endl;
   outputFile << "C2G1="
        << std::setw(22) << std::setprecision(15) << C2G(0,0)
        << std::setw(22) << std::setprecision(15) << C2G(0,1)
        << std::setw(22) << std::setprecision(15) << C2G(0,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << C2G(1,0)
        << std::setw(22) << std::setprecision(15) << C2G(1,1)
        << std::setw(22) << std::setprecision(15) << C2G(1,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << C2G(2,0)
        << std::setw(22) << std::setprecision(15) << C2G(2,1)
        << std::setw(22) << std::setprecision(15) << C2G(2,2)
        << ";" << std::endl;
   //outputFile <<  c2g << ";" <<  std::endl; // precision is too low
   outputFile << std::endl;
   outputFile << "A=" 
        << std::setw(22) << std::setprecision(15) << AA(0,0)
        << std::setw(22) << std::setprecision(15) << AA(0,1)
        << std::setw(22) << std::setprecision(15) << AA(0,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << AA(1,0)
        << std::setw(22) << std::setprecision(15) << AA(1,1)
        << std::setw(22) << std::setprecision(15) << AA(1,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << AA(2,0)
        << std::setw(22) << std::setprecision(15) << AA(2,1)
        << std::setw(22) << std::setprecision(15) << AA(2,2)
        << ";" << std::endl;
   outputFile << std::endl << std::endl;
   outputFile << "x0="
      << std::setw(21) << std::setprecision(15) << x0(0)
      << std::setw(21) << std::setprecision(15) << x0(1)
      << std::setw(21) << std::setprecision(15) << x0(2)
      << ";" << std::endl;
   outputFile << "periodicFaceIDs= 0 1 2 3 4 5 " << ";" << std::endl;
   outputFile << std::endl;

   outputFile << "solidSolutionNoiseMode=" << solidSolutionNoiseMode
      << "; # 0=no noise, 1= read noise, 2=compute noise" << std::endl;
   outputFile << "stackingFaultNoiseMode=" << stackingFaultNoiseMode
      << ";" << std::endl;
   if (! lattice.compare("bcc")) // .compare() returns 0 if they're equal
   {
      outputFile << "dislocationMobilityType='BCC';" << std::endl;
   }
   else if (! lattice.compare("fcc"))
   {
      outputFile << "dislocationMobilityType='FCC';" << std::endl;
   }
   outputFile << std::endl;
   //outputFile << "solidSolutionNoiseFile_xz=" << ";" << std::endl;
   //outputFile << "solidSolutionNoiseFile_yz=" << ";" << std::endl;
      
   return;
}

void model::AtomDisplacementGenerator::runGlideSteps( const size_t& stepsToRun)
{
   if ( DC == nullptr) readDefectiveCrystal();
   if ( DC->DN == nullptr)
   {
      std::cout << "error: dd2md::AtomDisplacementGenerator::runGlideSteps()"
         << " DefectiveCrystal::DefectiveNetwork not yet instantiated."
         << std::endl;
      return;
   }
   if ( DC->externalLoadController == nullptr)
   {
      std::cout << "error: dd2md::AtomDisplacementGenerator::runGlideSteps()"
        << " DC->externalLoadController isn't instantiated" << std::endl;
      return;
   }

   ddBase->simulationParameters.Nsteps
      = ddBase->simulationParameters.runID + stepsToRun;
   DC->runGlideSteps();
   return;
}


PYBIND11_MODULE( dd2md, m) {
   namespace py = pybind11;
   m.doc() = "TODO: revise m.doc() in dd2md.cpp";
   py::class_<model::AtomDisplacementGenerator>( m, "AtomDisplacementGenerator")
      .def( py::init([](
                  const std::string& modelibFolderPath
                  ){
               // lambda function that returns an instantiation
               return std::unique_ptr< model::AtomDisplacementGenerator>(
                  new model::AtomDisplacementGenerator(
                     modelibFolderPath
                     )
                  );
            }),
            py::arg("dddFolderPath").none(false)
            )
      .def("computeDisplacements",
            &model::AtomDisplacementGenerator::computeDisplacements
            //py::arg("lammpsDataFilePath").none(false)
            )
      .def("applyDisplacements",
            &model::AtomDisplacementGenerator::applyDisplacements
            )
      .def("computeDisplacements2",
            &model::AtomDisplacementGenerator::computeDisplacements2
            //py::arg("lammpsDataFilePath").none(false)
            )
      .def("getBurgersMagnitude",
            &model::AtomDisplacementGenerator::getBurgersMagnitude
          )
      .def("readBurgersMagnitude",
            &model::AtomDisplacementGenerator::readBurgersMagnitude,
            py::arg("materialFilePath").none(false)
          )
      .def("readLammpsConfiguration",
            &model::AtomDisplacementGenerator::readLammpsConfigurationFile,
            py::arg("lammpsDataFilePath").none(false)
          )
      .def("writeDisplacementsToFile",
            &model::AtomDisplacementGenerator::writeDisplacementsToFile,
            py::arg("outputFile").none(false)
            )
      .def("writeConfigurationToFile",
            &model::AtomDisplacementGenerator::writeConfigurationToFile,
            py::arg("outputFile").none(false)
            )
      .def("getPatchGlidePlanes",
            &model::AtomDisplacementGenerator::getPatchGlidePlanes
          )
      .def("regeneratePolycrystalFile",
            &model::AtomDisplacementGenerator::regeneratePolycrystalFile,
            py::kw_only(),
            py::arg( "c2g").none(false),
            py::arg( "lattice").none(false),
            py::arg( "material").none(false),
            py::arg( "meshFilePath").none(false)
          )
      .def("regenerateMicrostructure",
            &model::AtomDisplacementGenerator::regenerateMicrostructure
          )
      .def("readDefectiveCrystal",
            &model::AtomDisplacementGenerator::readDefectiveCrystal
          )
      .def("readLammpsBounds",
            &model::AtomDisplacementGenerator::readLammpsBounds,
            py::arg("lammpsDataFilePath").none(false)
          )
      .def("readConfiguration",
            &model::AtomDisplacementGenerator::readConfiguration,
            py::arg("runID").none(false)
          )
      .def("runGlideSteps",
            &model::AtomDisplacementGenerator::runGlideSteps,
            py::arg("stepsToRun").none(false)
          )
      //.def("dislocationPlasticDisplacement",
      //      static_cast<typename model::AtomDisplacementGenerator::VectorDim (model::AtomDisplacementGenerator::*)(const double&,const double&,const double&) const>(&model::AtomDisplacementGenerator::dislocationPlasticDisplacement)
      //      )
      //.def("getDisplacedPoints",
      //      &model::AtomDisplacementGenerator::getFieldPoints
      //      )
      ;
   //PYBIND11_NUMPY_DTYPE( );
}
#endif
