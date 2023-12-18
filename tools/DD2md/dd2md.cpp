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
          //Eigen::Matrix<double,3,3>::Identity(),
          debugFlag);
   std::cout << "lmpDeformationMatrix used for lammpsReader.readLmpStream():\n" << lmpDeformationMatrix << std::endl; // debug
   std::cout << "lammps boundaries (before readLmpStream): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
   if ( lammpsReader.readLmpStream(
            lammpsBoxBounds,
            lammpsTiltFactors,
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
   std::cout << "lammps boundaries (after readLmpStream): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
   // Shift atoms so that modelib slip systems and lammps slip systems
   //  coincide.
   shift_atoms();
   enforce_periodicity();

   std::cout << "lammps boundaries (readLammpsConfigurationFile): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
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
   return;
}


void model::AtomDisplacementGenerator::readConfigIO()
{
   if ( ddBase == nullptr)
   {
      std::cout << "error model::AtomDisplacementGenerator::readConfigIO()"
        << " called while ddBase == nullptr" << std::endl;
      return;
   }
   configIO = std::make_unique<DDconfigIO<3>>(
               ddBase->simulationParameters.traitsIO.evlFolder
               );
   configFields = std::make_unique<DDconfigFields<3>>(
                                     *ddBase, *configIO);
   return;
}

void model::AtomDisplacementGenerator::readDefectiveCrystal()
{
   if ( ddBase == nullptr) readddBase();
   readConfigIO();

   DC = std::unique_ptr<DefectiveCrystalType>(
         new DefectiveCrystalType( *ddBase)
         );
   return;
}

void model::AtomDisplacementGenerator::specifyLoops(
         const std::string& tag,
         const std::vector<int>& periodicLoopSlipSystemIDsIn,
         const std::vector<double>& periodicLoopRadii,
         const std::vector<long int>& loopSegmentCountsIn,
         const pybind11::array_t<double,
                  pybind11::array::c_style | pybind11::array::forcecast>&
                  loopCentersIn
      )
{
   if ( ddBase == nullptr) readddBase();

   // copy Nx3 numpy array of points in 3-D to Eigen equivalents
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   periodicLoopCenters.resize( loopCentersIn.shape()[0], loopCentersIn.shape()[1]);
   auto loopCentersInBuf = loopCentersIn.unchecked<2>();
   for ( ssize_t ii=0; ii < loopCentersIn.shape()[0]; ++ii)
      for ( ssize_t jj=0; jj < loopCentersIn.shape()[1]; ++jj)
         periodicLoopCenters( ii, jj) = loopCentersInBuf( ii, jj);

   // cast vectors of python integers to vector<int> (python int isn't int)
   std::vector<int> periodicLoopSlipSystemIDs( periodicLoopSlipSystemIDsIn.size(), 0);
   for ( size_t ii=0; ii < periodicLoopSlipSystemIDsIn.size(); ++ii)
      periodicLoopSlipSystemIDs[ ii] = static_cast<int>( periodicLoopSlipSystemIDsIn[ ii]);

   std::vector<long int> periodicLoopSides( loopSegmentCountsIn.size(), 0);
   for ( size_t ii=0; ii < loopSegmentCountsIn.size(); ++ii)
      periodicLoopSides[ ii] = static_cast<long int>( loopSegmentCountsIn[ ii]);

   // numpy loopRadii input is automatically seen as vector<double>

   // dummy instantiations
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   std::vector<double> prismaticLoopRadii;
   std::vector<double> prismaticLoopSteps;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PeriodicLoop", // type
               "individual", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn, // double
               0, //periodicLoopSegmentCountIn, // long int
               0.0, //periodicLoopRadiusDistributionMeanIn // double
               0.0, //periodicLoopRadiusDistributionStdIn // double
               0.0, //periodicDipoleTargetDensityIn, // double
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopSteps,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, // prismaticLoopRadiiMean
               1.0, // prismaticLoopRadiiStd
               1.0, // prismaticLoopStepMean
               1.0 // prismaticLoopStepStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}

void model::AtomDisplacementGenerator::setOutputPath( const std::string& outputPath)
{ // TODO: change return value to allow error checking
   if ( ddBase == nullptr) readddBase();
   std::cout << "prior evlFolder: "
      << ddBase->simulationParameters.traitsIO.evlFolder << std::endl; // debug
   std::cout << "assigning " << outputPath + "/evl" << " to evlFolder" << std::endl;
   ddBase->simulationParameters.traitsIO.evlFolder = outputPath + "/evl";
   ddBase->simulationParameters.traitsIO.auxFolder = outputPath + "/evl";
   ddBase->simulationParameters.traitsIO.fFolder = outputPath + "/F";
   ddBase->simulationParameters.traitsIO.fFile = outputPath + "/F/F_0.txt";
   ddBase->simulationParameters.traitsIO.flabFile = outputPath + "/F/F_labels.txt";
   std::cout << "results:" << std::endl
      << " ddBase->simulationParameters.traitsIO.evlFolder: " << ddBase->simulationParameters.traitsIO.evlFolder
      << std::endl
      << " ddBase->simulationParameters.traitsIO.auxFolder: "
      << ddBase->simulationParameters.traitsIO.auxFolder
      << std::endl
      << " ddBase->simulationParameters.traitsIO.fFolder: "
      << ddBase->simulationParameters.traitsIO.fFolder
      << std::endl;

   // output folders are at:
   //  ddBase->simulationParameters.traitsIO.evlFolder
   //  ddBase->simulationParameters.traitsIO.auxFolder
   //  ddBase->simulationParameters.traitsIO.fFolder
   configIO->clear();
   configIO = std::make_unique<DDconfigIO<3>>(outputPath + "/evl");
   return;
}

void model::AtomDisplacementGenerator::specifyLoopDensitiesPerSlipSystem(
         const std::string& tag,
         const std::map<int, double>& periodicLoopTargetDensitiesPerSlipSystemIn,
         const long int& periodicLoopSegmentCountIn,
         const double& periodicLoopRadiusDistributionMeanIn,
         const double& periodicLoopRadiusDistributionStdIn
      )
{
   if ( ddBase == nullptr) readddBase();
   // instantiate empty placeholders
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   std::vector<long int> loopSegmentCountsIn;
   //pybind11::array_t<double,
   //         pybind11::array::c_style | pybind11::array::forcecast>
   //         periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<double> prismaticLoopSteps;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PeriodicLoop", // type
               "densitiesPerSlipSystem", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, // periodicLoopTargetDensityIn,
               periodicLoopSegmentCountIn,
               periodicLoopRadiusDistributionMeanIn,
               periodicLoopRadiusDistributionStdIn,
               0.0, //periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystemIn,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopSteps,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, // prismaticLoopRadiiMean
               1.0, // prismaticLoopRadiiStd
               1.0, // prismaticLoopStepMean
               1.0 // prismaticLoopStepStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug


   return;
}

void model::AtomDisplacementGenerator::specifyLoopDensity(
         const std::string& tag,
         const double& periodicLoopTargetDensityIn,
         const long int& periodicLoopSegmentCountIn,
         const double& periodicLoopRadiusDistributionMeanIn,
         const double& periodicLoopRadiusDistributionStdIn
      )
{
   if ( ddBase == nullptr) readddBase();

   // instantiate empty placeholders
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   std::vector<long int> loopSegmentCountsIn;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   //pybind11::array_t<double,
   //         pybind11::array::c_style | pybind11::array::forcecast>
   //         periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<double> prismaticLoopSteps;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PeriodicLoop", // type
               "density", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               periodicLoopTargetDensityIn,
               periodicLoopSegmentCountIn,
               periodicLoopRadiusDistributionMeanIn,
               periodicLoopRadiusDistributionStdIn,
               0.0, //periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopSteps,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, // prismaticLoopRadiiMean
               1.0, // prismaticLoopRadiiStd
               1.0, // prismaticLoopStepMean
               1.0 // prismaticLoopStepStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}
/**********************************************************************/
void model::AtomDisplacementGenerator::specifyPrismaticLoops(
         const std::string& tag,
         const std::vector<int>& prismaticLoopSlipSystemIDsIn,
         const std::vector<double>& prismaticLoopRadiiIn,
         const pybind11::array_t<double,
                  pybind11::array::c_style | pybind11::array::forcecast>&
                  prismaticLoopCentersIn,
         const std::vector<double>& prismaticLoopStepsIn
      )
{
   if ( ddBase == nullptr) readddBase();

   // copy Nx3 numpy array of points in 3-D to Eigen equivalents
   Eigen::Matrix<double, Eigen::Dynamic, 3> prismaticLoopCenters;
   prismaticLoopCenters.resize( prismaticLoopCentersIn.shape()[0], prismaticLoopCentersIn.shape()[1]);
   auto prismaticLoopCentersInBuf = prismaticLoopCentersIn.unchecked<2>();
   for ( ssize_t ii=0; ii < prismaticLoopCentersIn.shape()[0]; ++ii)
      for ( ssize_t jj=0; jj < prismaticLoopCentersIn.shape()[1]; ++jj)
         prismaticLoopCenters( ii, jj) = prismaticLoopCentersInBuf( ii, jj);

   // cast vectors of python integers to vector<int> (python int isn't int)
   std::vector<int> prismaticLoopSlipSystemIDs( prismaticLoopSlipSystemIDsIn.size(), 0);
   for ( size_t ii=0; ii < prismaticLoopSlipSystemIDsIn.size(); ++ii)
      prismaticLoopSlipSystemIDs[ ii] = static_cast<int>( prismaticLoopSlipSystemIDsIn[ ii]);

   // numpy loopRadii input is automatically seen as vector<double>

   // dummy instantiations
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<long int> periodicLoopSides;
   std::vector<double> periodicLoopRadii;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> periodicLoopCenters;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PrismaticLoop", // type
               "individual", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn, // double
               0, //periodicLoopSegmentCountIn, // long int
               0.0, //periodicLoopRadiusDistributionMeanIn // double
               0.0, //periodicLoopRadiusDistributionStdIn // double
               0.0, //periodicDipoleTargetDensityIn, // double
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadiiIn,
               prismaticLoopStepsIn,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, // prismaticLoopRadiiMean
               1.0, // prismaticLoopRadiiStd
               1.0, // prismaticLoopStepMean
               1.0 // prismaticLoopStepStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}

void model::AtomDisplacementGenerator::specifyPrismaticLoopDensitiesPerSlipSystem(
         const std::string& tag,
         const std::map<int, double>& prismaticLoopTargetDensitiesPerSlipSystemIn,
         const double& prismaticLoopRadiiMeanIn,
         const double& prismaticLoopRadiiStdIn,
         const double& prismaticLoopStepMeanIn,
         const double& prismaticLoopStepStdIn
      )
{
   if ( ddBase == nullptr) readddBase();
   // instantiate empty placeholders
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   std::vector<long int> loopSegmentCountsIn;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   //pybind11::array_t<double,
   //         pybind11::array::c_style | pybind11::array::forcecast>
   //         periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<double> prismaticLoopSteps;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PrismaticLoop", // type
               "densitiesPerSlipSystem", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, // periodicLoopTargetDensityIn,
               0, //periodicLoopSegmentCountIn,
               1.0, //periodicLoopRadiusDistributionMeanIn,
               1.0, //periodicLoopRadiusDistributionStdIn,
               0.0, //periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopSteps,
               prismaticLoopTargetDensitiesPerSlipSystemIn,
               prismaticLoopRadiiMeanIn,
               prismaticLoopRadiiStdIn,
               prismaticLoopStepMeanIn,
               prismaticLoopStepStdIn
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug


   return;
}

void model::AtomDisplacementGenerator::specifyPrismaticLoopDensity(
         const std::string& tag,
         const double& prismaticLoopTargetDensityIn,
         const double& prismaticLoopRadiiMeanIn,
         const double& prismaticLoopRadiiStd,
         const double& prismaticLoopStepMeanIn,
         const double& prismaticLoopStepStd
      )
{
   if ( ddBase == nullptr) readddBase();

   // instantiate empty placeholders
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   std::vector<long int> loopSegmentCountsIn;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   //pybind11::array_t<double,
   //         pybind11::array::c_style | pybind11::array::forcecast>
   //         periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<double> prismaticLoopSteps;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PrismaticLoop", // type
               "density", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn,
               0, // periodicLoopSegmentCountIn,
               0.0, //periodicLoopRadiusDistributionMeanIn,
               0.0, //periodicLoopRadiusDistributionStdIn,
               0.0, //periodicDipoleTargetDensityIn,
               prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopSteps,
               prismaticLoopTargetDensitiesPerSlipSystem,
               prismaticLoopRadiiMeanIn,
               prismaticLoopRadiiStd,
               prismaticLoopStepMeanIn,
               prismaticLoopStepStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}

/**********************************************************************/
void model::AtomDisplacementGenerator::specifyDipoleDensity(
               const std::string& tag,
               const double& periodicDipoleTargetDensityIn
      )
{
   if ( ddBase == nullptr) readddBase();

   // dummy instantiations. TODO: find another way to implement this class
   Eigen::Matrix<double, Eigen::Dynamic, 3> points;
   std::vector<double> heights;
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> exitFaceIDs;//( exitFaceIDsIn.size(), 0);
   std::vector<int> nodes;//( nodesIn.size(), 0);
   std::vector<double> glideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   Eigen::Matrix<double,Eigen::Dynamic,3> periodicLoopCenters;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<double> prismaticLoopSteps;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
         std::shared_ptr<model::MicrostructureSpecification>(
            //new model::PeriodicDipoleIndividualSpecification(
            new model::MicrostructureSpecification(
               "PeriodicDipole", // type
               "density", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs, // std::vector<int>
               // dipoles
               exitFaceIDs, // std::vector<int>
               points, // Eigen::Matrix
               heights, // std::vector<double>
               nodes, // std::vector<int>
               glideSteps, // std::vector<double>
               // loops
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn,
               0, //periodicLoopSegmentCountIn,
               0.0, //periodicLoopRadiusDistributionMeanIn
               0.0, //periodicLoopRadiusDistributionStdIn
               periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopSteps,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, //prismaticLoopRadiiMean // not used
               1.0, //prismaticLoopRadiiStd // not used
               1.0, //prismaticLoopStepMean // not used
               1.0 //prismaticLoopStepStd // not used
               )
            )
         );

   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}


void model::AtomDisplacementGenerator::specifyDipoles(
      const std::string& tag,
      const std::vector<int>& periodicDipoleSlipSystemIDsIn,
      const std::vector<int>& exitFaceIDsIn,
      const pybind11::array_t<double,
               pybind11::array::c_style | pybind11::array::forcecast>&
               pointsIn,
      const std::vector<double>& heights,
      const std::vector<int>& nodesIn,
      const std::vector<double>& glideSteps
      )
{
   if ( ddBase == nullptr) readddBase();

   std::cout << "pointsIn.shape()[0]: " << pointsIn.shape()[0] << std::endl; // debug
   std::cout << "pointsIn.shape()[1]: " << pointsIn.shape()[1] << std::endl; // debug
   std::cout << "periodicDipoleSlipSystemIDsIn.size(): " << periodicDipoleSlipSystemIDsIn.size() << std::endl; // debug
   std::cout << "exitFaceIDsIn.size(): " << exitFaceIDsIn.size() << std::endl; // debug

   Eigen::Matrix<double, Eigen::Dynamic, 3> points;
   points.resize( pointsIn.shape()[0], pointsIn.shape()[1]);
   auto pointsInBuf = pointsIn.unchecked<2>();
   for ( ssize_t ii=0; ii < pointsIn.shape()[0]; ++ii)
      for ( ssize_t jj=0; jj < pointsIn.shape()[1]; ++jj)
         points( ii, jj) = pointsInBuf( ii, jj);
   std::cout << "points.size(): " << points.size() << std::endl; // debug

   std::vector<int> periodicDipoleSlipSystemIDs( periodicDipoleSlipSystemIDsIn.size(), 0);
   for ( size_t ii=0; ii < periodicDipoleSlipSystemIDsIn.size(); ++ii)
      periodicDipoleSlipSystemIDs[ ii] = static_cast<int>( periodicDipoleSlipSystemIDsIn[ ii]);

   std::vector<int> exitFaceIDs( exitFaceIDsIn.size(), 0);
   for ( size_t ii=0; ii < exitFaceIDsIn.size(); ++ii)
      exitFaceIDs[ ii] = static_cast<int>( exitFaceIDsIn[ ii]);

   std::vector<int> nodes( nodesIn.size(), 0);
   for ( size_t ii=0; ii < nodesIn.size(); ++ii)
      nodes[ ii] = static_cast<int>( nodesIn[ ii]);

   // dummy instantiations. TODO: find another way to implement this class
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   Eigen::Matrix<double,Eigen::Dynamic,3> periodicLoopCenters;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<double> prismaticLoopSteps;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
         std::shared_ptr<model::MicrostructureSpecification>(
            //new model::PeriodicDipoleIndividualSpecification(
            new model::MicrostructureSpecification(
               "PeriodicDipole", // type
               "individual", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDsIn, // std::vector<int>
               // dipoles
               exitFaceIDs, // std::vector<int>
               points, // Eigen::Matrix
               heights, // std::vector<double>
               nodes, // std::vector<int>
               glideSteps, // std::vector<double>
               // loops
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn,
               0, //periodicLoopSegmentCountIn,
               0.0, //periodicLoopRadiusDistributionMeanIn
               0.0, //periodicLoopRadiusDistributionStdIn
               0.0, //periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopSteps,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, //prismaticLoopRadiiMean // not used
               1.0, // prismaticLoopRadiiStd
               1.0, //prismaticLoopStepMean // not used
               1.0 // prismaticLoopStepStd
               )
            )
         );

   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug
   return;
}

void model::AtomDisplacementGenerator::generateMicrostructure()
{
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "generating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << std::endl; // debug

   if ( microstructureSpecifications.size() == 0)
   {
      std::cout << "error: generateMicrostructure() requires defects to be specified via respective functions" << std::endl;
      return;
   }
   for (const auto& spec : microstructureSpecifications) // debug
   { // debug
      std::cout << "spec->tag: " << spec->tag << std::endl; // debug
      std::cout << "spec->microstructureType: " // debug
         << spec->microstructureType<< std::endl; // debug
      std::cout << "spec->style: " << spec->style<< std::endl; // debug
   } // debug
   model::MicrostructureGeneratorInMem mg( *ddBase, microstructureSpecifications);
   return;
}

void model::AtomDisplacementGenerator::clearMicrostructureSpecifications()
{
   microstructureSpecifications.clear();
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

   std::cout << "lammps boundaries (regenerateMicrostructure): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
   return;
}

void model::AtomDisplacementGenerator::readLammpsBox(
      const std::string& lammpsFilePath
      )
{
   LmpReader lammpsReader(
          lammpsFilePath,
          //1.0,
          1e-10/burgersMagnitude, // b_SI:2.556e-10 #scaleFactor for atom positions
          Eigen::Matrix<double,3,3>::Identity(),
          debugFlag);
   if ( lammpsReader.readLmpStreamBox(
            lammpsBoxBounds,
            lammpsTiltFactors
            ) != EXIT_SUCCESS)
   {
      std::cout << "error: cannot read boundaries from lammps data file."
         << std::endl;
   }
   std::cout << "lammpsFilePath: " << lammpsFilePath << std::endl; // debug
   std::cout << "lammpsBoxBounds: "; // debug
   for ( size_t ii=0; ii < lammpsBoxBounds.size(); ++ii) // debug
   { // debug
      std::cout << lammpsBoxBounds[ii] << ", "; // debug
   } // debug
   std::cout << std::endl; // debug
   std::cout << "lammpsTiltFactors: "; // debug
   for ( size_t ii=0; ii < lammpsTiltFactors.size(); ++ii) // debug
   { // debug
      std::cout << lammpsTiltFactors[ii] << ", "; // debug
   } // debug
   std::cout << std::endl; // debug
   std::cout << "lammps boundaries (readLammpsBox): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
   // shift lammps box boundaries so that the origin is at (0,0,0)
   //  as required by modelib
   lammpsBoxBounds[1] = lammpsBoxBounds[1] - lammpsBoxBounds[0];
   lammpsBoxBounds[0] = 0.0;
   lammpsBoxBounds[3] = lammpsBoxBounds[3] - lammpsBoxBounds[2];
   lammpsBoxBounds[2] = 0.0;
   lammpsBoxBounds[5] = lammpsBoxBounds[5] - lammpsBoxBounds[4];
   lammpsBoxBounds[4] = 0.0;
   std::cout << "shifted lammps boundaries (readLammpsBox): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug

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
   return *configIO;
}

model::DDconfigIO<3>& model::AtomDisplacementGenerator::config()
{
   return *configIO;
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
   configIO->read( runID);
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

// TODO: create a version of this that takes arguments specifying a subset of all loops and a collection of points
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
      std::cout << "error: computeDisplacements() called while fieldPoints"
         << " is empty" << std::endl;
      return;
   }

   if ( debugFlag) std::cout << "Computing DislocationDisplacement at field points..." << std::endl;
   //DC->DN->displacement( fieldPoints);
   if (displacements.size() != 0) displacements.clear();
   if (displacements.size() != fieldPoints.size())
   {
      displacements.resize( fieldPoints.size());
   }
   for (unsigned int k=0; k < fieldPoints.size(); ++k)
   {
      //if (
      //      std::isnan( fieldPoints[k](0))
      //      || std::isnan( fieldPoints[k](1))
      //      || std::isnan( fieldPoints[k](2))
      //   )
      //{
      //   throw  std::runtime_error("error computeDisplacements(): fieldpoints["+std::to_string(k)+"] is nan");
      //}
      displacements[k] = DC->DN->displacement( fieldPoints[k]);
      //if (
      //      std::isnan( displacements[k](0))
      //      || std::isnan( displacements[k](1))
      //      || std::isnan( displacements[k](2))
      //   )
      //{
      //   throw  std::runtime_error("error computeDisplacements(): displacements["+std::to_string(k)+"] is nan");
      //}
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
   if (displacements.size() != 0) displacements.clear();
   if (displacements.size() != fieldPoints.size())
   {
      displacements.resize( fieldPoints.size());
   }
   for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   {
      displacements[ ii] = dislocationPlasticDisplacement(
           fieldPoints[ ii]
           );
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
   enforce_periodicity();
   std::cout << "lammps boundaries (applyDisplacements): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
   return;
}

void model::AtomDisplacementGenerator::applyExternalElasticField()
{
   // Elastic displacement container
   if ( debugFlag) std::cout<<"applyExternalElasticField(): Computing the displacement associated with the external stressZX...\n"<<std::flush;
   const double gamma_ZX = 8.e8/ddBase->poly.mu_SI;
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
                  =  py::array_t<double, pybind11::array::c_style>();
               polygons[ loopNumber][ patchNumber].resize(
                     {vertexCounts[ loopNumber][ patchNumber], 3} );
               py::buffer_info polygonBuf
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
   if ( ddBase == nullptr) readddBase();
   std::cout << "ddBase->poly.b_SI: " << ddBase->poly.b_SI << std::endl; // debug
   return ddBase->poly.b_SI/1e-10;
}

void model::AtomDisplacementGenerator::writeConfigurationToFile(
      const std::string& outputFilePath)
{
   std::cout << "lammps boundaries (writeConfigurationToFile pre): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
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
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[0] << " "
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[1]
      //<< std::setprecision(16) << "0 "
      //<< std::setprecision(16) << (burgersMagnitude/1e-10) * (lammpsBoxBounds[1] - lammpsBoxBounds[0])
      << " xlo xhi" << std::endl
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[2] << " "
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[3]
      //<< std::setprecision(16) << "0 "
      //<< std::setprecision(16) << (burgersMagnitude/1e-10) * (lammpsBoxBounds[3] - lammpsBoxBounds[2])
      << " ylo yhi" << std::endl
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[4] << " "
      << std::setprecision(16)  << (burgersMagnitude/1e-10) * lammpsBoxBounds[5]
      //<< std::setprecision(16) << "0 "
      //<< std::setprecision(16) << (burgersMagnitude/1e-10) * (lammpsBoxBounds[5] - lammpsBoxBounds[4])
      << " zlo zhi" << std::endl;

   if (( lammpsTiltFactors.size() >= 3)
         && (
            (lammpsTiltFactors[0] != 0)
            ||
            (lammpsTiltFactors[1] != 0)
            ||
            (lammpsTiltFactors[2] != 0)
            ))
   {
      outputFile
         << std::setprecision(16)
         << (burgersMagnitude/1e-10) * lammpsTiltFactors[0]
         << " " << std::setprecision(16)
         << (burgersMagnitude/1e-10) * lammpsTiltFactors[1]
         << " " << std::setprecision(16)
         << (burgersMagnitude/1e-10) * lammpsTiltFactors[2]
         << " xy xz yz" << std::endl;
   }

   outputFile << std::endl << "Masses" << std::endl << std::endl;
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
        << std::setprecision(16)
        << (fieldPoints[k]).transpose() * ddBase->poly.b_SI /1e-10
          //+ (fieldPoints[k]).transpose() * DC->DN->poly.b_SI /1e-10
        //<< (lmpDeformationMatrixInverse * fieldPoints[k]).transpose() * DC->DN->poly.b_SI /1e-10
        //  + (lmpDeformationMatrixInverse * fieldPoints[k]).transpose() * DC->DN->poly.b_SI /1e-10
        //
        << "\n";
   }
   if ( debugFlag)
      std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<std::endl;

   std::cout << "lammps boundaries (writeConfigurationToFile post): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
   return;
}

py::array_t<double, py::array::c_style>
   model::AtomDisplacementGenerator::getDisplacementsNumpy()
{
   py::array_t<double, py::array::c_style> displacementsTmp;
   int atomCount( displacements.size());
   displacementsTmp.resize( { atomCount, 3});

   py::buffer_info displacementsBuf = displacementsTmp.request();
   double* displacementsPtr = static_cast<double*>( displacementsBuf.ptr);

   for ( int ii=0; ii < atomCount; ++ii)
   {
      displacementsPtr[0 + 3*ii]
         = displacements[ii](0) * ddBase->poly.b_SI /1e-10;
      displacementsPtr[1 + 3*ii]
         = displacements[ii](1) * ddBase->poly.b_SI /1e-10;
      displacementsPtr[2 + 3*ii]
         = displacements[ii](2) * ddBase->poly.b_SI /1e-10;
   }
   return displacementsTmp;
}

py::array_t<double, py::array::c_style>
   model::AtomDisplacementGenerator::getLammpsBoxPBCVectors()
{
   py::array_t<double, py::array::c_style> lammpsBoxPBCVectorsTmp;
   if ( lammpsBoxPBCVectors.size() == 0)
   {
      std::cout << "error, getLammpsBoxPBCVectors() was called before lammpsBoxPBCVectors were evaluated" << std::endl;
      return lammpsBoxPBCVectorsTmp;
   }

   int vectorCount( lammpsBoxPBCVectors.size());
   lammpsBoxPBCVectorsTmp.resize( {vectorCount, 3});

   py::buffer_info lammpsBoxPBCVectorsBuf = lammpsBoxPBCVectorsTmp.request();
   double* lammpsBoxPBCVectorsPtr = static_cast<double*>( lammpsBoxPBCVectorsBuf.ptr);
   for ( int ii=0; ii < vectorCount; ++ii)
   {
      lammpsBoxPBCVectorsPtr[0 + 3*ii] = lammpsBoxPBCVectors[ii](0);
      lammpsBoxPBCVectorsPtr[1 + 3*ii] = lammpsBoxPBCVectors[ii](1);
      lammpsBoxPBCVectorsPtr[2 + 3*ii] = lammpsBoxPBCVectors[ii](2);
   }

   return lammpsBoxPBCVectorsTmp;
}

py::array_t<double, py::array::c_style>
   model::AtomDisplacementGenerator::getLammpsBoxVertices()
{
   py::array_t<double, py::array::c_style> lammpsBoxVerticesTmp;
   if ( lammpsBoxVertices.size() == 0)
   {
      std::cout << "error, getLammpsBoxVertices() was called before lammpsBoxVertices were evaluated" << std::endl;
      return lammpsBoxVerticesTmp;
   }

   int vertexCount( lammpsBoxVertices.size());
   lammpsBoxVerticesTmp.resize( {vertexCount, 3});

   py::buffer_info lammpsBoxVerticesBuf = lammpsBoxVerticesTmp.request();
   double* lammpsBoxVerticesPtr = static_cast<double*>( lammpsBoxVerticesBuf.ptr);
   for ( int ii=0; ii < vertexCount; ++ii)
   {
      lammpsBoxVerticesPtr[0 + 3*ii] = lammpsBoxVertices[ii](0);
      lammpsBoxVerticesPtr[1 + 3*ii] = lammpsBoxVertices[ii](1);
      lammpsBoxVerticesPtr[2 + 3*ii] = lammpsBoxVertices[ii](2);
   }

   return lammpsBoxVerticesTmp;
}

std::map< size_t, model::AtomDisplacementGenerator::VectorDim> //pybind11::array_t<double, pybind11::array::c_style> >
   model::AtomDisplacementGenerator::getDisplacementsMap()
{
   std::map< size_t, model::AtomDisplacementGenerator::VectorDim>
      displacementsTmp;
   displacementsTmp.clear();

   for ( size_t ii=0; ii < displacements.size(); ++ii)
   {
      displacementsTmp[ pointIDs[ii]]
         = displacements[ii] * ddBase->poly.b_SI /1e-10;
   }

   return displacementsTmp;
}

void model::AtomDisplacementGenerator::writeDisplacementsToLammpsFile(
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
         << "AtomDisplacementGenerator::writeDisplacementsToLammpsFile"
         << " called while displacements is empty" << std::endl;
      return;
   }
   if ( displacements.size() != atomTypes.size())
   {
      std::cout << "error "
         << "AtomDisplacementGenerator::writeDisplacementsToLammpsFile"
         << " displacements.size() != atomType.size()) " << std::endl;
      return;
   }

   outputFile << "# LAMMPS data file written by DD2MD" << std::endl << std::endl;
   outputFile << displacements.size() << " atoms" << std::endl << std::endl;
   outputFile << masses.size() << " atom types" << std::endl << std::endl;
   outputFile
      // << std::setw(16)
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[0] << " "
      //<< std::setw(16)
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[1]
      << " xlo xhi" << std::endl
      //<< std::setw(16)
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[2] << " "
      //<< std::setw(16)
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[3]
      << " ylo yhi" << std::endl
      //<< std::setw(16)
      << std::setprecision(16) << (burgersMagnitude/1e-10) * lammpsBoxBounds[4] << " "
      //<< std::setw(16)
      << std::setprecision(16)  << (burgersMagnitude/1e-10) * lammpsBoxBounds[5]
      << " zlo zhi" << std::endl;

   // write the lammps triclinic box tilt factors
   if ((lammpsTiltFactors.size() == 3)
         && (
            lammpsTiltFactors[0] != 0.0
            ||
            lammpsTiltFactors[1] != 0.0
            ||
            lammpsTiltFactors[2] != 0.0
            )
         )
   {
      outputFile
         << std::setprecision(16)  << (burgersMagnitude/1e-10) * lammpsTiltFactors[0] << " "
         << std::setprecision(16)  << (burgersMagnitude/1e-10) * lammpsTiltFactors[1] << " "
         << std::setprecision(16)  << (burgersMagnitude/1e-10) * lammpsTiltFactors[2] << " "
         << " xy xz yz" << std::endl;
   }
   std::cout << std::endl;
   outputFile << "Masses" << std::endl << std::endl;
   for ( const auto& mm : masses)
   {
      outputFile << mm.first << " " << mm.second << std::endl;
   }

   outputFile << std::endl << "Atoms # atomic" << std::endl << std::endl;

   //for (unsigned int k=0; k < displacements.size(); ++k) // debug
   //{ // debug
   //   if ( std::isnan( displacements[k](0)) // debug
   //         || std::isnan( displacements[k](1)) // debug
   //         || std::isnan( displacements[k](2)) // debug
   //         ) // debug
   //   { // debug
   //      throw std::runtime_error("error writeDisplacementsToLammpsFile: displacements["+std::to_string(k)+"] is nan"); // debug
   //   } // debug
   //} // debug

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
        << std::setprecision(16)
        << (displacements[k]).transpose() * ddBase->poly.b_SI /1e-10
          //+ (displacements[k]).transpose() * DC->DN->poly.b_SI /1e-10
        //<< (lmpDeformationMatrixInverse * displacements[k]).transpose() * DC->DN->poly.b_SI /1e-10
        //  + (lmpDeformationMatrixInverse * displacements[k]).transpose() * DC->DN->poly.b_SI /1e-10
        << "\n";
   }
   if ( debugFlag)
      std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<std::endl;

   std::cout << "lammps boundaries (writeDisplacementsToLammpsFilepost): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
   return;
}

long int model::AtomDisplacementGenerator::getCurrentStep()
{
   if ( ddBase == nullptr)
   {
      std::cout << "error, model::AtomDisplacementGenerator::getCurrentStep(): ddBase == nullptr. returning -1" << std::endl;
      return -1;
   }
   return ddBase->simulationParameters.runID;
}

void model::AtomDisplacementGenerator::setCurrentStep( const long int& step)
{
   //if ( DC == nullptr)
   //{
   //   std::cout << "error: cannot setCurrentStep until "
   //      << "DefectiveCrystal is instantiated" << std::endl;
   //   return;
   //}

   ddBase->simulationParameters.runID = step;
   if ( DC != nullptr)
   {
      std::cout << "model::AtomDisplacementGenerator::setCurrentStep(): DC != nullptr, updating load and calling managerRestart()" << std::endl; // debug
      DC->externalLoadController->update( DC->plasticStrain());
      ddBase->simulationParameters.manageRestart();
   }
   return;
}


void model::AtomDisplacementGenerator::shift_atoms(
      //fieldPointsType& fieldPoints
      )
{
   std::cout << "lammps boundaries (shift_atoms init): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
   // iterate over all positions to find the lowest
   if ( fieldPoints.size() <= 0)
   {
      std::cout << "AtomDisplacementGenerator::shift_atoms(): "
         << "passed an empty set of atoms ..."
         << std::endl;
      return;
   }

   //for (unsigned int k=0; k < fieldPoints.size(); ++k) // debug
   //{ // debug
   //   if ( std::isnan( fieldPoints[k](0)) // debug
   //         || std::isnan( fieldPoints[k](1)) // debug
   //         || std::isnan( fieldPoints[k](2)) // debug
   //         ) // debug
   //   { // debug
   //      throw std::runtime_error("error shift_atoms(): received fieldPoints["+std::to_string(k)+"] is nan"); // debug
   //   } // debug
   //} // debug

   //if ((lammpsTiltFactors.size() != 3)
   //      || !(
   //         lammpsTiltFactors[0] != 0.0
   //         ||
   //         lammpsTiltFactors[1] != 0.0
   //         ||
   //         lammpsTiltFactors[2] != 0.0
   //         )
   //      )
   //{

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
   // shift lammps box boundaries
   lammpsBoxBounds[1] = lammpsBoxBounds[1] - lammpsBoxBounds[0];
   lammpsBoxBounds[0] = 0.0;
   lammpsBoxBounds[3] = lammpsBoxBounds[3] - lammpsBoxBounds[2];
   lammpsBoxBounds[2] = 0.0;
   lammpsBoxBounds[5] = lammpsBoxBounds[5] - lammpsBoxBounds[4];
   lammpsBoxBounds[4] = 0.0;

   //}

   //for (unsigned int k=0; k < fieldPoints.size(); ++k) // debug
   //{ // debug
   //   if ( std::isnan( fieldPoints[k](0)) // debug
   //         || std::isnan( fieldPoints[k](1)) // debug
   //         || std::isnan( fieldPoints[k](2)) // debug
   //         ) // debug
   //   { // debug
   //      throw std::runtime_error("error shift_atoms(): turned fieldPoints["+std::to_string(k)+"] into nan"); // debug
   //   } // debug
   //} // debug

   // else
   // {
   // TODO: figure out how to handle atoms exceeding triclinic box bounds
   // }

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
   displacementVector = (C2Ginv)*displacementVector;
   std::cout << "displacementVector: \n" << displacementVector << std::endl; // debug

   // TODO: incorporate atomic crystallographic orientation rather than use this
   for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   {
      //if ( std::isnan( displacementVector(0)) // debug
      //      || std::isnan( displacementVector(1)) // debug
      //      || std::isnan( displacementVector(2)) // debug
      //      ) // debug
      //{ // debug
      //   throw std::runtime_error("error shift_atoms(): displacementVector for ["+std::to_string(ii)+"] is nan"); // debug
      //} // debug
      fieldPoints[ii] += displacementVector;
   }

   //// translate atoms whose positions exceed the periodic boundaries
   //for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   //{
   //   while ( fieldPoints[ii](0) < lammpsDeformedBoxBounds[0]) // xlo
   //      fieldPoints[ii](0) += lammpsDeformedBoxDimensions[0];
   //   while ( fieldPoints[ii](0) >= lammpsDeformedBoxBounds[1]) // xhi
   //      fieldPoints[ii](0) -= lammpsDeformedBoxDimensions[0];

   //   while ( fieldPoints[ii](1) < lammpsDeformedBoxBounds[2]) // ylo
   //      fieldPoints[ii](1) += lammpsDeformedBoxDimensions[1];
   //   while ( fieldPoints[ii](1) >= lammpsDeformedBoxBounds[3]) // yhi
   //      fieldPoints[ii](1) -= lammpsDeformedBoxDimensions[1];

   //   while ( fieldPoints[ii](2) < lammpsDeformedBoxBounds[4]) // zlo
   //      fieldPoints[ii](2) += lammpsDeformedBoxDimensions[2];
   //   while ( fieldPoints[ii](2) >= lammpsDeformedBoxBounds[5]) // zhi
   //      fieldPoints[ii](2) -= lammpsDeformedBoxDimensions[2];
   //}
   std::cout << "lammps boundaries (shift_atoms post): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
   return;
}

void model::AtomDisplacementGenerator::enforce_periodicity()
{
   // iterate over all positions to find the lowest
   if ( fieldPoints.size() <= 0)
   {
      std::cout << "AtomDisplacementGenerator::enforce_periodicity(): "
         << "passed an empty set of atoms ..."
         << std::endl;
      return;
   }

   // Iterate over all fieldPoints and translate them if they exceed
   //  lammpsBoxBounds.
   double periodX, periodY, periodZ;
   periodX = lammpsBoxBounds[1] - lammpsBoxBounds[0];
   periodY = lammpsBoxBounds[3] - lammpsBoxBounds[2];
   periodZ = lammpsBoxBounds[5] - lammpsBoxBounds[4];

   for ( size_t ii=0; ii < fieldPoints.size(); ++ii)
   {
      // translate when points are below lower bounds
      while ( fieldPoints[ii](0) < lammpsBoxBounds[0])
      {
         fieldPoints[ii](0) += periodX;
      }
      while ( fieldPoints[ii](1) < lammpsBoxBounds[2])
      {
         fieldPoints[ii](1) += periodY;
      }
      while ( fieldPoints[ii](2) < lammpsBoxBounds[4])
      {
         fieldPoints[ii](2) += periodZ;
      }
      // translate when points are above upper bounds
      while ( fieldPoints[ii](0) > lammpsBoxBounds[1])
      {
         fieldPoints[ii](0) -= periodX;
      }
      while ( fieldPoints[ii](1) > lammpsBoxBounds[3])
      {
         fieldPoints[ii](1) -= periodY;
      }
      while ( fieldPoints[ii](2) > lammpsBoxBounds[5])
      {
         fieldPoints[ii](2) -= periodZ;
      }
   }
   return;
}

void model::AtomDisplacementGenerator::regeneratePolycrystalFile(
            const py::array_t<double,
               py::array::c_style | py::array::forcecast>
               grain1globalX1In, // vector in crystal coordinates to align with the x1 global axis (e.g. [11-1] for bcc)
            const py::array_t<double,
               py::array::c_style | py::array::forcecast>
               grain1globalX3In, // vector in crystal coordinates to align with the x3 global axis (e.g. [101] for bcc)
            const py::array_t<int,
               py::array::c_style | py::array::forcecast>
               boxScalingIn, // number of unit cells per box direction
            const py::array_t<int,
               py::array::c_style | py::array::forcecast>
               boxEdges1In, // in crystal coordinates
            const py::array_t<int,
               py::array::c_style | py::array::forcecast>
               boxEdges2In, // in crystal coordinates
            const py::array_t<int,
               py::array::c_style | py::array::forcecast>
               boxEdges3In, // in crystal coordinates
            const double& TT,
            const int& enablePartials,
            const std::string& latticeIn,
            const std::string& materialIn,
            const std::string& meshFilePathIn,
            //const py::array_t<double,
            //   py::array::c_style | py::array::forcecast>
            //   X0In, // mesh nodes are mapped to x=F*(X-X0)
            const py::array_t<double,
               py::array::c_style | py::array::forcecast>
               periodicFaceIDsIn
            )
{
   if ( DC != nullptr) DC = nullptr;
   if ( ddBase != nullptr) ddBase = nullptr;

   std::cout << "lammps boundaries (regeneratePolycrystalFile): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug

   auto grain1globalX1Np = grain1globalX1In.unchecked<1>();
   auto grain1globalX3Np = grain1globalX3In.unchecked<1>();
   //auto X0Np = X0In.unchecked<1>();
   auto periodicFaceIDsNp = periodicFaceIDsIn.unchecked<1>();
   auto boxScalingNp = boxScalingIn.unchecked<1>();
   auto boxEdges1Np = boxEdges1In.unchecked<1>();
   auto boxEdges2Np = boxEdges2In.unchecked<1>();
   auto boxEdges3Np = boxEdges3In.unchecked<1>();
   if (! ((grain1globalX1Np.shape(0) == 3) ))
   {
      std::cout << "error: regeneratePolycrystalFile( ) requires "
         << "grain1globalX1 to be a 3 element numpy array consisting of "
         << "a vector in crystal coordinates to align with the first "
         << "axis of the MoDELib box"
         << std::endl;
      return;
   }
   if (! ((grain1globalX3Np.shape(0) == 3) ))
   {
      std::cout << "error: regeneratePolycrystalFile( ) requires "
         << "grain1globalX3 to be a 3 element numpy array consisting of "
         << "a vector in crystal coordinates to align with the third "
         << "axis of the MoDELib box"
         << std::endl;
      return;
   }
   VectorDim grain1globalX1;
   grain1globalX1 << grain1globalX1Np[0], grain1globalX1Np[1], grain1globalX1Np[2];
   VectorDim grain1globalX3;
   grain1globalX3 << grain1globalX3Np[0], grain1globalX3Np[1], grain1globalX3Np[2];
   //if (! ((X0Np.shape(0) == 3) ))
   //{
   //   std::cout << "error: regeneratePolycrystalFile( ) requires "
   //      << "X0 to be a 3 element numpy array consisting of "
   //      << "a vector in mesh coordinates to traslate the mesh before "
   //      << "it is rescaled to the MoDELib box (x=F(X-X0))"
   //      << std::endl;
   //   return;
   //}
   if (! ((boxEdges1Np.shape(0) == 3) && (boxEdges2Np.shape(0) == 3)
            && (boxEdges3Np.shape(0) == 3)
            ))
   {
      std::cout << "error: regeneratePolycrystalFile( ) requires "
         << "boxEdges1, boxEdges2, boxEdges3, to each be a 3 element "
         << "numpy array consisting of vectors in crystal "
         << "coordinates parallel to the MoDELib "
         << "box edges"
         << std::endl;
      return;
   }

   if (! ((periodicFaceIDsNp.shape(0) == 6) ))
   {
      std::cout << "error: regeneratePolycrystalFile( ) currently "
         << "periodicFaceIDs to be a 6 element numpy array consisting of "
         << "integers enumerating periodic faces."
         << "regeneratePolycrystalFile() does not currently accept "
         << "a mixture of periodic and aperiodic boundaries."
         << std::endl;
      return;
   }
   if ( TT < 0)
   {
      std::cout << "error, regeneratePolycrystalFile() unacceptable "
         << "value for temperature " << TT << std::endl;
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
      if ( materialIn == std::string(mat)) materialIsAcceptable = true;
   }
   if ( ! ( latticeIsAcceptable && materialIsAcceptable))
   {
      std::cout << "error: unacceptable lattice or material: "
         << lattice << ", " << materialIn << std::endl;
      std::cout << "  acceptable lattice or materials are: ";
      for ( const auto& mat : acceptableMaterials) std::cout << mat << ", ";
      for ( const auto& lat : acceptableLattices) std::cout << lat << ", ";
      return;
   }

   lattice = latticeIn; // assign member of AtomDisplacementGenerator

   std::cout << "lattice: " << lattice << ", material: " << materialIn << std::endl; // debug

   VectorDim x3Crossx1( grain1globalX3.cross( grain1globalX1));
   x3Crossx1 /= sqrt( x3Crossx1.dot( x3Crossx1));

   // normalize vectors
   grain1globalX1 /= sqrt( grain1globalX1.dot( grain1globalX1));
   grain1globalX3 /= sqrt( grain1globalX3.dot( grain1globalX3));

   // recalculate C2G, an orthonormal rotation matrix
   //MatrixDim C2G;
   C2G << grain1globalX1(0), grain1globalX1(1), grain1globalX1(2),
         x3Crossx1(0), x3Crossx1(1), x3Crossx1(2),
         grain1globalX3(0), grain1globalX3(1), grain1globalX3(2);

   std::cout << "C2G:\n" << C2G << std::endl; // debug

   C2Ginv = C2G.inverse();
   std::cout << "C2Ginv:\n" << C2Ginv << std::endl; // debug

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

   // instantiate box scaling in multiples of burgers vector magnitude
   // deltaX is in \AA if using lammps metal units
   //std::vector<double> boxScaling(
   //      { deltaX / (burgersMagnitude/1e-10),
   //      deltaY / (burgersMagnitude/1e-10),
   //      deltaZ / (burgersMagnitude/1e-10)});

   std::cout << "boxScaling: " << boxScalingNp[0] << ", " << boxScalingNp[1] << ", " << boxScalingNp[2] << std::endl; // debug

   MatrixDim lammpsBoxEdgesCrystalVectors; // vectors in columns
   lammpsBoxEdgesCrystalVectors
      << boxEdges1Np[0], boxEdges2Np[0], boxEdges3Np[0],
         boxEdges1Np[1], boxEdges2Np[1], boxEdges3Np[1],
         boxEdges1Np[2], boxEdges2Np[2], boxEdges3Np[2];
   //// ii'th row is the direction of the ii'th box edge in global coordinates
   //lammpsBoxEdges
   //   << deltaX, 0, 0,
   //      lammpsTiltFactors[0], deltaY, 0,
   //      lammpsTiltFactors[1], lammpsTiltFactors[2], deltaZ;
   //MatrixDim lammpsBoxEdgeUnitVectors( lammpsBoxEdges);
   //// renomalize box edges (cols of lammpsBoxEdgeUnitVectors) in global coordinates
   //double tmpMag = 0;
   //for ( size_t ii=0; ii < 3; ++ii)
   //{
   //   tmpMag = sqrt( lammpsBoxEdgeUnitVectors( ii, 0) * lammpsBoxEdgeUnitVectors( ii, 0)
   //                + lammpsBoxEdgeUnitVectors( ii, 1) * lammpsBoxEdgeUnitVectors( ii, 1)
   //                + lammpsBoxEdgeUnitVectors( ii, 2) * lammpsBoxEdgeUnitVectors( ii, 2));
   //   if ( ( tmpMag > 0) && !( isinf( tmpMag)))
   //   {
   //      lammpsBoxEdgeUnitVectors( ii, 0) /= tmpMag;
   //      lammpsBoxEdgeUnitVectors( ii, 1) /= tmpMag;
   //      lammpsBoxEdgeUnitVectors( ii, 2) /= tmpMag;
   //   }
   //   else
   //   {
   //      std::cout << "error, regeneratePolycrystalFile(): "
   //         << "column " << ii << " of lammpsBoxEdgeUnitVectors has magnitude "
   //         << tmpMag << " <= 0" << std::endl;
   //      return;
   //   }
   //}
   // ii'th column is the direction of the ii'th box edge in global coordinates
   MatrixDim lammpsBoxEdges;
   lammpsBoxEdges // in global coordinates
      << deltaX, lammpsTiltFactors[0], lammpsTiltFactors[1],
         0,                    deltaY, lammpsTiltFactors[2],
         0,                         0,               deltaZ;
   std::cout << "lammps boundaries (regeneratePolycrystalFile2): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug

   std::cout << "lammpsBoxEdges (regeneratePolycrystalFile2):\n" << lammpsBoxEdges <<  std::endl;// debug

   //MatrixDim lammpsBoxEdgeUnitVectors( lammpsBoxEdges);
   // renomalize box edges (cols of lammpsBoxEdgeUnitVectors) in global coordinates
   //double tmpMag = 0;
   //for ( size_t ii=0; ii < 3; ++ii)
   //{
   //   tmpMag = sqrt( lammpsBoxEdgeUnitVectors( 0, ii) * lammpsBoxEdgeUnitVectors( 0, ii)
   //                + lammpsBoxEdgeUnitVectors( 1, ii) * lammpsBoxEdgeUnitVectors( 1, ii)
   //                + lammpsBoxEdgeUnitVectors( 2, ii) * lammpsBoxEdgeUnitVectors( 2, ii));
   //   if ( ( tmpMag > 0) && !( isinf( tmpMag)))
   //   {
   //      lammpsBoxEdgeUnitVectors( 0, ii) /= tmpMag;
   //      lammpsBoxEdgeUnitVectors( 1, ii) /= tmpMag;
   //      lammpsBoxEdgeUnitVectors( 2, ii) /= tmpMag;
   //   }
   //   else
   //   {
   //      std::cout << "error, regeneratePolycrystalFile(): "
   //         << "column " << ii << " of lammpsBoxEdgeUnitVectors has magnitude "
   //         << tmpMag << " <= 0" << std::endl;
   //      return;
   //   }
   //}
   //std::cout << "lammpsBoxEdges:\n" << lammpsBoxEdges << std::endl; // debug
   //std::cout << "lammpsBoxEdgeUnitVectors:\n" << lammpsBoxEdgeUnitVectors << std::endl; // debug

   MatrixDim AA;
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

   // Find lattice vectors (columns of modelibBoxEdges) aligned to lammpsBoxEdgesCrystalVectors to be scaled and become columns of FF
   MatrixDim AAinv;
   AAinv = AA.inverse();
   //MatrixDim C2Ginv;
   //C2Ginv = C2G.inverse();

   // interpret lammps box edge vectors
   //MatrixDim lammpsBoxEdgesCrystalUnits( C2Ginv*(lammpsBoxEdgeUnitVectors));
   // try to transform lammpsBoxEdgesCrystalUnits cols into integer vectors
   // This assumes the user has prepared the lammps box edges to be aligned with crystallographic directions which can be expressed as integers.
   //std::cout << "lammpsBoxEdgesCrystalUnits:\n" // debug
   //   << lammpsBoxEdgesCrystalUnits<< std::endl; // debug

   //// TODO: fix
   //// TODO: fix
   //// TODO: fix
   //// TODO: fix
   //for ( size_t jj=0; jj < 3; ++jj)
   //{
   //   double tmpDblMax( 0.0);
   //   tmpDblMax = lammpsBoxEdgesCrystalUnits( Eigen::seq(0,Eigen::last), jj).array().abs().maxCoeff();
   //   for ( size_t ii=0; ii < 3; ++ii)
   //   {
   //      if ( abs( lammpsBoxEdgesCrystalUnits( ii, jj))
   //            > std::numeric_limits<double>::epsilon())
   //      {
   //         lammpsBoxEdgesCrystalUnits( ii, jj)
   //            = lammpsBoxEdgesCrystalUnits( ii, jj)/ tmpDblMax
   //      }
   //      else
   //      {
   //         lammpsBoxEdgesCrystalUnits( ii, jj) = 0;
   //      }
   //   }
   //}
   std::cout << "lammpsBoxEdgesCrystalVectors:\n" // debug
      << lammpsBoxEdgesCrystalVectors << std::endl; // debug
   MatrixDim BB( AAinv * lammpsBoxEdgesCrystalVectors);
   //MatrixDim BB( AAinv * C2Ginv*(lammpsBoxEdgesCrystalUnits.transpose())); // transform boxEdge vectors
   //MatrixDim modelibBoxEdges;//( MatrixDim::Zero());
   MatrixDim FF; // FF scales the mesh: y=F(x-x0), where x is the input mesh
   //std::cout << "AA:\n" << AA << std::endl; // debug
   //std::cout << "BB:\n" << BB << std::endl; // debug
   for ( size_t jj=0; jj < 3; ++jj)
   {
      VectorDim Bcol( BB( Eigen::seq(0,Eigen::last), jj));
      //std::cout << "Bcol:\n" << Bcol << std::endl; // debug
      Bcol /= Bcol.array().abs().maxCoeff(); // rescale box edge vector components to be within [-1,1] interval
      //std::cout << "Bcol rescaled:\n" << Bcol << std::endl; // debug
      Eigen::Matrix<int,3,1> nn; nn << 0,0,0;
      Eigen::Matrix<int,3,1> dd; dd << 1,1,1;
      std::pair<int,int> ff;
      for ( size_t ii=0; ii < 3; ++ii)
      {
         ff = limit_denominator( Bcol[ii], 100);
         nn(ii) = ff.first;
         dd(ii) = ff.second;
      }
      //std::cout << "nn:\n" << nn << std::endl; // debug
      //std::cout << "dd:\n" << dd << std::endl; // debug
      ////dp=np.prod(dd);
      int dp( 1); // product of denominators of box edge vector components
      for ( size_t ii=0; ii < 3; ++ii)
      {
         dp *= dd(ii);
      }
      //std::cout << "dp: " << dp << std::endl; // debug
      //nr=np.array([1,1,1], dtype=int)
      Eigen::Matrix<int,3,1> nr; nr << 1, 1, 1;
      //for i in range(0, 3):
      //    nr[ii]=nn[ii]*dp/dd[ii]
      for ( size_t ii=0; ii < 3; ++ii)
      {
         nr(ii) = nn(ii) * dp/dd(ii); // numerators multiplied by the denominators missing from their rational expression
      }
      //std::cout << "nr:\n" << nr << std::endl; // debug
      //nr=nr/np.gcd.reduce(nr)
      // divide out the greatest common divisor of nr's components
      int nrGcd( 1);
      nrGcd = std::gcd( nr(0), nr(1));
      nrGcd = std::gcd( nrGcd, nr(2));
      nr /= nrGcd;
      //L[:,jj]=self.A@nr.transpose()
      VectorDim nrDubl;
      nrDubl << nr(0), nr(1), nr(2);
      //std::cout << "nrDubl:\n" << nrDubl << std::endl; // debug
      VectorDim AAnr( AA*nrDubl); // L[:,jj]
      //std::cout << "AAnr:\n" << AAnr << std::endl; // debug
      //self.F[:,jj]=self.C2G@modelibBoxEdges[:,jj]*self.boxScaling[jj]
      //self.F[:,jj]=self.C2G@L[:,jj]*self.boxScaling[jj]
      // columns of modelibBoxEdges should be lattice vectors (in global coordinates) aligned to box edges
      VectorDim modelibBoxEdgesCol( (C2G*AAnr)* boxScalingNp[jj]); // transform lattice vectors aligned to box edges from crystal coordinates to global coordinates
      //std::cout << "modelibBoxEdgesCol:\n" << modelibBoxEdgesCol << std::endl; // debug
      for ( size_t ii=0; ii < 3; ++ii)
      {
         FF(ii,jj) = modelibBoxEdgesCol(ii);
      }
   }
   std::cout << "FF:\n" << FF << std::endl; // debug

   // FF is now determined. FF scales the mesh to a size and shape almost equal to that of the lammps box.
   // The lammps box now needs to be deformed a little to ensure the  modelib and lammps box sizes and shapes are identical.

   // Create a transformation to be applied to the atoms to align
   //  their strained atomic planes to the perfect slip systems of DDD.
   // Atomic planes will still need to be shifted to make atomic slip
   //  planes coincide with DDD slip planes.

   lmpDeformationMatrix = FF * (lammpsBoxEdges.inverse()) ;
   lmpDeformationMatrixInverse =  lmpDeformationMatrix.inverse();
   std::cout << "lmpDeformationMatrix :\n" << lmpDeformationMatrix // debug
      << std::endl; // debug
   std::cout << "lmpDeformationMatrixInverse :\n" // debug
      << lmpDeformationMatrixInverse // debug
      << std::endl; // debug
   std::cout << "lammpsBoxEdges:\n" << lammpsBoxEdges << std::endl;// debug
   std::cout << "lammps boundaries (regeneratePolycrystalFile3): " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug

   VectorDim edgeVector1, edgeVector2, edgeVector3;
   edgeVector1 << lammpsBoxBounds[1] - lammpsBoxBounds[0],
                  0,
                  0;
   edgeVector2 << lammpsTiltFactors[0], // xy
                  lammpsBoxBounds[3] - lammpsBoxBounds[2], // yhi - ylo
                  0;
   edgeVector3 << lammpsTiltFactors[1], // xz,
                  lammpsTiltFactors[2], // yz,
                  lammpsBoxBounds[5] - lammpsBoxBounds[4];// zhi - zlo

   lammpsBoxPBCVectors.clear();
   lammpsBoxPBCVectors.push_back( VectorDim());
   lammpsBoxPBCVectors.push_back( VectorDim());
   lammpsBoxPBCVectors.push_back( VectorDim());
   lammpsBoxPBCVectors[0] << edgeVector1(0),
                               edgeVector1(1),
                               edgeVector1(2);
   lammpsBoxPBCVectors[1] << edgeVector2(0),
                               edgeVector2(1),
                               edgeVector2(2);
   lammpsBoxPBCVectors[2] << edgeVector3(0),
                               edgeVector3(1),
                               edgeVector3(2);
   for (size_t ii=0; ii < lammpsBoxPBCVectors.size(); ++ii)
   {
      lammpsBoxPBCVectors[ii](0) *= burgersMagnitude/1e-10;
      lammpsBoxPBCVectors[ii](1) *= burgersMagnitude/1e-10;
      lammpsBoxPBCVectors[ii](2) *= burgersMagnitude/1e-10;
      lammpsBoxPBCVectors[ii] = lmpDeformationMatrix * lammpsBoxPBCVectors[ii];
   }

   lammpsBoxVertices.clear();
   lammpsBoxVertices.push_back( VectorDim());
   lammpsBoxVertices.push_back( VectorDim());
   lammpsBoxVertices.push_back( VectorDim());
   lammpsBoxVertices.push_back( VectorDim());
   lammpsBoxVertices.push_back( VectorDim());
   lammpsBoxVertices.push_back( VectorDim());
   lammpsBoxVertices.push_back( VectorDim());
   lammpsBoxVertices.push_back( VectorDim());
   lammpsBoxVertices[0] <<
            lammpsBoxBounds[0],
            lammpsBoxBounds[2],
            lammpsBoxBounds[4];
   lammpsBoxVertices[1] <<
            lammpsBoxBounds[0] + edgeVector1(0),
            lammpsBoxBounds[2] + edgeVector1(1),
            lammpsBoxBounds[4] + edgeVector1(2);
   lammpsBoxVertices[2] <<
            lammpsBoxBounds[0] + edgeVector1(0) + edgeVector2(0),
            lammpsBoxBounds[2] + edgeVector1(1) + edgeVector2(1),
            lammpsBoxBounds[4] + edgeVector1(2) + edgeVector2(2);
   lammpsBoxVertices[3] <<
            lammpsBoxBounds[0] + edgeVector2(0),
            lammpsBoxBounds[2] + edgeVector2(1),
            lammpsBoxBounds[4] + edgeVector2(2);
   lammpsBoxVertices[4] <<
            lammpsBoxBounds[0] + edgeVector2(0) + edgeVector3(0),
            lammpsBoxBounds[2] + edgeVector2(1) + edgeVector3(1),
            lammpsBoxBounds[4] + edgeVector2(2) + edgeVector3(2);
   lammpsBoxVertices[5] <<
            lammpsBoxBounds[0] + edgeVector3(0),
            lammpsBoxBounds[2] + edgeVector3(1),
            lammpsBoxBounds[4] + edgeVector3(2);
   lammpsBoxVertices[6] <<
            lammpsBoxBounds[0] + edgeVector1(0) + edgeVector3(0),
            lammpsBoxBounds[2] + edgeVector1(1) + edgeVector3(1),
            lammpsBoxBounds[4] + edgeVector1(2) + edgeVector3(2);
   lammpsBoxVertices[7] <<
            lammpsBoxBounds[0] + edgeVector1(0) + edgeVector2(0) + edgeVector3(0),
            lammpsBoxBounds[2] + edgeVector1(1) + edgeVector2(1) + edgeVector3(1),
            lammpsBoxBounds[4] + edgeVector1(2) + edgeVector2(2) + edgeVector3(2);
   for ( size_t ii=0; ii < lammpsBoxVertices.size(); ++ii)
   {
      lammpsBoxVertices[ii](0) *= burgersMagnitude/1e-10;
      lammpsBoxVertices[ii](1) *= burgersMagnitude/1e-10;
      lammpsBoxVertices[ii](2) *= burgersMagnitude/1e-10;
      lammpsBoxVertices[ii] = lmpDeformationMatrix * lammpsBoxVertices[ii];
   }

   ////////////////////////////////////////////////
   /// use the local lmpDeformationMatrix to calculate lammpsDeformedBox
   // the following should be unecessary. lammpsBoxBounds were alread deformed by the LmpReader which was instantiated using lmpDeformationMatrix, which should have been applied to all positions and bounds when the file was read.
   //// transform lammps box boundary coordinates
   //lammpsDeformedBoxBounds.clear();
   //lammpsDeformedBoxBounds.resize(6);
   //lammpsDeformedBoxBounds[0] // xlo
   //   =  lmpDeformationMatrix(0,0) * lammpsBoxBounds[0];
   //   //= (lmpDeformationMatrix * lmpAxis00)[0];

   //lammpsDeformedBoxBounds[1] // xhi
   //   =  lmpDeformationMatrix(0,0) * lammpsBoxBounds[1];
   //   //= (lmpDeformationMatrix * lmpAxis01)[0];

   //lammpsDeformedBoxBounds[2] // ylo
   //   =  lmpDeformationMatrix(1,1) * lammpsBoxBounds[2];
   //   //= (lmpDeformationMatrix * lmpAxis10)[1];

   //lammpsDeformedBoxBounds[3] // yhi
   //   =  lmpDeformationMatrix(1,1) * lammpsBoxBounds[3];
   //   //= (lmpDeformationMatrix * lmpAxis11)[1];

   //lammpsDeformedBoxBounds[4] // zlo
   //   =  lmpDeformationMatrix(2,2) * lammpsBoxBounds[4];
   //   //= (lmpDeformationMatrix * lmpAxis20)[2];

   //lammpsDeformedBoxBounds[5] // zhi
   //   =  lmpDeformationMatrix(2,2) * lammpsBoxBounds[5];
   //   //= (lmpDeformationMatrix * lmpAxis21)[2];

   //std::cout << "lammpsDeformedBoxBounds:\n"; // debug
   //for ( const auto& tmpDbl : lammpsDeformedBoxBounds) // debug
   //{ // debug
   //   std::cout << tmpDbl << ", "; // debug
   //} // debug
   //std::cout << std::endl; // debug

   //// lammpsDeformedTiltFactors
   //lammpsDeformedTiltFactors.clear();
   //lammpsDeformedTiltFactors.resize(3);
   //lammpsDeformedTiltFactors[0] // xy
   //   = lmpDeformationMatrix(0,0) * lammpsTiltFactors[0];
   //lammpsDeformedTiltFactors[1] // xz
   //   = lmpDeformationMatrix(0,0) * lammpsTiltFactors[1];
   //lammpsDeformedTiltFactors[2] // yz
   //   = lmpDeformationMatrix(1,1) * lammpsTiltFactors[2];

   //std::cout << "lammpsDeformedTiltFactors:\n"; // debug
   //for ( const auto& tmpDbl : lammpsDeformedTiltFactors) // debug
   //{ // debug
   //   std::cout << tmpDbl << ", "; // debug
   //} // debug
   //std::cout << std::endl; // debug


   //lammpsDeformedBoxDimensions.clear();
   //lammpsDeformedBoxDimensions.resize(3);
   //lammpsDeformedBoxDimensions[0] = lammpsDeformedBoxBounds[1]
   //                                  - lammpsDeformedBoxBounds[0];
   //lammpsDeformedBoxDimensions[1] = lammpsDeformedBoxBounds[3]
   //                                  - lammpsDeformedBoxBounds[2];
   //lammpsDeformedBoxDimensions[2] = lammpsDeformedBoxBounds[5]
   //                                  - lammpsDeformedBoxBounds[4];
   //std::cout << "lammpsDeformedBoxBounds: ("
   //   << lammpsDeformedBoxBounds[0] << ", "
   //   << lammpsDeformedBoxBounds[1] << ", "
   //   << lammpsDeformedBoxBounds[2] << ", "
   //   << lammpsDeformedBoxBounds[3] << ", "
   //   << lammpsDeformedBoxBounds[4] << ", "
   //   << lammpsDeformedBoxBounds[5] << ") " << std::endl; // debug
   //if ( lammpsDeformedTiltFactors.size() == 3)
   //{
   //   std::cout << "lammpsDeformedTiltFactors: ("
   //      << lammpsDeformedTiltFactors[0] << ", "
   //      << lammpsDeformedTiltFactors[1] << ", "
   //      << lammpsDeformedTiltFactors[2] << ") " << std::endl; // debug
   //}
   ////std::cout << "lammpsDeformedBoxDimensions: ("
   ////   << lammpsDeformedBoxDimensions[0] << ", "
   ////   << lammpsDeformedBoxDimensions[1] << ", "
   ////   << lammpsDeformedBoxDimensions[2] << ")" << std::endl; // debug

   double x0x, x0y, x0z;
   x0x = lammpsBoxBounds[0]/deltaX; // lammpsBoxBounds[0] is xlo, deltaX = xhi -xlo
   x0y = lammpsBoxBounds[2]/deltaY; // lammpsBoxBounds[2] is ylo
   x0z = lammpsBoxBounds[4]/deltaZ; // lammpsBoxBounds[4] is zlo
   VectorDim x0;
   x0 << x0x, x0y, x0z; // scaling of the mesh is y=A(x-x0)

   std::cout << "X0: " << x0x << ", " << x0y << ", " << x0z << std::endl; // debug

   //////////////////////////////////////

   outputFile << "materialFile=" << materialIn << ".txt;" << std::endl;
   outputFile << "enablePartials=" << enablePartials << ";"
      << std::endl;
   outputFile << "absoluteTemperature = "
      << TT << "; # [K] simulation temperature"
      << std::endl;
   outputFile << "meshFile=" << meshFilePathIn << ";"
      << std::endl;
   outputFile << "C2G1="
        << std::setw(23) << std::setprecision(16) << C2G(0,0)
        << std::setw(23) << std::setprecision(16) << C2G(0,1)
        << std::setw(23) << std::setprecision(16) << C2G(0,2)
        << std::endl
        << std::setw(23) << std::setprecision(16) << C2G(1,0)
        << std::setw(23) << std::setprecision(16) << C2G(1,1)
        << std::setw(23) << std::setprecision(16) << C2G(1,2)
        << std::endl
        << std::setw(23) << std::setprecision(16) << C2G(2,0)
        << std::setw(23) << std::setprecision(16) << C2G(2,1)
        << std::setw(23) << std::setprecision(16) << C2G(2,2)
        << ";" << std::endl;
   outputFile << std::endl;
   outputFile << "F="
        << std::setw(23) << std::setprecision(16) << FF(0,0)
        << std::setw(23) << std::setprecision(16) << FF(0,1)
        << std::setw(23) << std::setprecision(16) << FF(0,2)
        << std::endl
        << std::setw(23) << std::setprecision(16) << FF(1,0)
        << std::setw(23) << std::setprecision(16) << FF(1,1)
        << std::setw(23) << std::setprecision(16) << FF(1,2)
        << std::endl
        << std::setw(23) << std::setprecision(16) << FF(2,0)
        << std::setw(23) << std::setprecision(16) << FF(2,1)
        << std::setw(23) << std::setprecision(16) << FF(2,2)
        << ";" << std::endl;
   outputFile << std::endl << std::endl;

   outputFile << "X0="
      << std::setw(23) << std::setprecision(16) << x0(0)
      << std::setw(23) << std::setprecision(16) << x0(1)
      << std::setw(23) << std::setprecision(16) << x0(2)
      << ";" << std::endl;
   outputFile << "periodicFaceIDs= "
      << periodicFaceIDsNp(0) << " "
      << periodicFaceIDsNp(1) << " "
      << periodicFaceIDsNp(2) << " "
      << periodicFaceIDsNp(3) << " "
      << periodicFaceIDsNp(4) << " "
      << periodicFaceIDsNp(5) << ";" << std::endl;
   outputFile << std::endl;

   outputFile << "solidSolutionNoiseMode=" << solidSolutionNoiseMode
      << "; # 0=no noise, 1= read noise, 2=compute noise" << std::endl;
   //outputFile << "solidSolutionNoiseFile_xz=" << ";" << std::endl;
   //outputFile << "solidSolutionNoiseFile_yz=" << ";" << std::endl;
   outputFile << "stackingFaultNoiseMode=" << stackingFaultNoiseMode
      << ";" << std::endl;
   if (! lattice.compare("bcc")) // .compare() returns 0 if they're equal
   {
      //outputFile << "dislocationMobilityType='BCC';" << std::endl;
      outputFile << "dislocationMobilityType=default;" << std::endl;
   }
   else if (! lattice.compare("fcc"))
   {
      //outputFile << "dislocationMobilityType='FCC';" << std::endl;
      outputFile << "dislocationMobilityType=default;" << std::endl;
   }
   outputFile << "gridSize=256 256; # size of grid on the glide plane"
      << std::endl;
   outputFile << "gridSpacing_SI=1e-10 1e-10; # [m] spacing of grid on the glide plane"
      << std::endl;
   //outputFile << "spreadLstress_A=1; # add comment" << std::endl;
   //outputFile << "seed=0; # add comment" << std::endl;
   outputFile << std::endl;

   std::cout << "lammpsBoxBounds 1530: " <<  std::endl;// debug
   for ( const auto& bd : lammpsBoxBounds) std::cout << bd << ", "; // debug
   std::cout << std::endl; // debug
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

std::pair<int,int> model::AtomDisplacementGenerator::limit_denominator(
      const double& xxIn,
      const int& max_denominator
      )
{
   // adapted from https://github.com/python/cpython/blob/66c0d0ac8c55be9c973be1189b0e9ffcfdfb35a4/Lib/fractions.py#L234
   //  and https://en.wikipedia.org/wiki/Continued_fraction

   //std::cout << "limit_denominator(" << xxIn << ", " << max_denominator << ")" << std::endl; // debug
   if ( max_denominator < 1)
   {
      throw std::runtime_error("limit_denominator(), error: argument max_denominator must be greater than 1");
   }
   if ( std::isinf( xxIn))
   {
      throw std::runtime_error("limit_denominator(), error: given inf as argument");
   }
   if ( std::isnan( xxIn))
   {
      throw std::runtime_error("limit_denominator(), error: given NaN as argument");
   }
   if ( abs( xxIn) < 1.0/max_denominator)
   {
      std::cout << "limit_denominator(" << xxIn << ", " << max_denominator << ") returning " << 0 << "/" << 1 << std::endl; // debug
      return std::make_pair<int,int>( 0, 1);
   }
   double remainder( xxIn - floor( xxIn));
   if ( remainder == 0 ) // xxIn is an integer
   {
      std::cout << "limit_denominator(" << xxIn << ", " << max_denominator << ") returning " << static_cast<int>(xxIn) << "/" << 1 << std::endl; // debug
      return std::make_pair<int,int>( static_cast<int>( xxIn), 1);
   }
   int numerator0, denominator0, numerator1, denominator1, numerator2, denominator2;
   double invReciprocal;
   numerator0 = static_cast<int>( floor( xxIn));
   denominator0 = 1;
   numerator1 = 1;
   denominator1 = 0;

   int coefficient;
   std::cout << "limit_denominator("
      << xxIn << ", " << max_denominator << ") entering while loop" << std::endl; // debug
   while ( true)
   {
      invReciprocal = 1.0/remainder;
      coefficient = static_cast<int>( floor( invReciprocal));
      remainder = invReciprocal - coefficient;

       // numerator2 is older than numerator1
      numerator2 = numerator1;
      numerator1 = numerator0;
      denominator2 = denominator1;
      denominator1 = denominator0;

      numerator0 = coefficient * numerator1 + numerator2;
      denominator0 = coefficient * denominator1 + denominator2;

      //std::cout // debug
      //   << "coefficient " << coefficient << ", " // debug
      //   << "remainder " << remainder << ", " // debug
      //   << "numerator0 " << numerator0 << ", " // debug
      //   << "numerator1 " << numerator1 << ", " // debug
      //   << "numerator2 " << numerator2 << ", " // debug
      //   << "denominator0 " << denominator0 << ", " // debug
      //   << "denominator1 " << denominator1 << ", " // debug
      //   << "denominator2 " << denominator2 << std::endl; // debug

      //if (( xxIn - invReciprocal == 0)
      if (( remainder <= std::numeric_limits<double>::epsilon())
            || denominator0 > max_denominator)
      {
         //std::cout << "numerator/denominator: " // debug
         //   << numerator1 << "/" << denominator1 << std::endl; // debug
         break;
      }
   }
   if ( remainder <= std::numeric_limits<double>::epsilon())
   {
      std::pair<int,int> pp{ numerator0, denominator0};
      std::cout << "limit_denominator(" << xxIn << ", " << max_denominator << ") returning " << pp.first << "/" << pp.second << std::endl; // debug
      return pp;
   }
   int kk( floor( (max_denominator - denominator0)/denominator1));
   int bound1_n, bound1_d;
   bound1_n = numerator0 + kk*numerator1;
   bound1_d = denominator0 + kk*denominator1;
   if (
         abs( static_cast<double>(numerator1) / static_cast<double>(denominator1) - xxIn)
         <= abs((static_cast<double>(bound1_n)/static_cast<double>(bound1_d)) - xxIn)
         )
   {
      std::cout << "(n1/d1)";
      //return std::make_pair<int,int>( numerator1, denominator1);
      std::pair<int,int> pp{ numerator1, denominator1};
      std::cout << "limit_denominator(" << xxIn << ", " << max_denominator << ") returning " << pp.first << "/" << pp.second << std::endl; // debug
      return pp;
   }
   else
   {
      std::cout << "(bound1)";
      //else return std::make_pair<int,int>( bound1_n, bound1_d);
      std::pair<int,int> pp{ bound1_n, bound1_d};
      std::cout << "limit_denominator(" << xxIn << ", " << max_denominator << ") returning " << pp.first << "/" << pp.second << std::endl; // debug
      return pp;
   }
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
      .def("writeDisplacementsToLammpsFile",
            &model::AtomDisplacementGenerator::writeDisplacementsToLammpsFile,
            py::arg("outputFile").none(false)
            )
      .def("getDisplacementsNumpy",
            &model::AtomDisplacementGenerator::getDisplacementsNumpy
            )
      .def("getDisplacementsMap",
            &model::AtomDisplacementGenerator::getDisplacementsMap
            )
      .def("getLammpsBoxVertices",
            &model::AtomDisplacementGenerator::getLammpsBoxVertices
            )
      .def("getLammpsBoxPBCVectors",
            &model::AtomDisplacementGenerator::getLammpsBoxPBCVectors
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
            py::arg( "grain1globalX1").none(false), // in crystal coord
            py::arg( "grain1globalX3").none(false), // in crystal coord
            py::arg( "boxScaling").none(false),
            py::arg( "boxEdges1").none(false), // in rows of a matrix, in crystal coord
            py::arg( "boxEdges2").none(false), // in rows of a matrix, in crystal coord
            py::arg( "boxEdges3").none(false), // in rows of a matrix, in crystal coord
            py::arg( "T").none(false), // absolute temperature [K]
            py::arg( "enablePartials"), // True or False. Default: False
            py::arg( "lattice").none(false),
            py::arg( "material").none(false),
            py::arg( "meshFilePath").none(false),
            //py::arg( "X0").none(false),
            py::arg( "periodicFaceIDs").none(false)
          )
      .def("specifyLoops",
            &model::AtomDisplacementGenerator::specifyLoops,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "slipSystemIDs").none(false),
            py::arg( "loopRadii").none(false),
            py::arg( "loopSegmentCounts").none(false),
            py::arg( "loopCenters").none(false)
          )
      .def("specifyDipoles",
            &model::AtomDisplacementGenerator::specifyDipoles,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "slipSystemIDs").none(false),
            py::arg( "exitFaceIDs").none(false),
            py::arg( "points").none(false),
            py::arg( "heights").none(false),
            py::arg( "dipoleNodes").none(false),
            py::arg( "dipoleGlideSteps").none(false)
          )
      .def("specifyLoopDensity",
            &model::AtomDisplacementGenerator::specifyLoopDensity,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "loopDensity").none(false),
            py::arg( "loopSegmentCount").none(false),
            py::arg( "loopRadiusMean").none(false),
            py::arg( "loopRadiusStd").none(false)
          )
      .def("specifyLoopDensitiesPerSlipSystem",
            &model::AtomDisplacementGenerator::specifyLoopDensitiesPerSlipSystem,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "loopDensitiesPerSlipSystem").none(false),
            py::arg( "loopSegmentCount").none(false), // keep
            py::arg( "loopRadiusMean").none(false), // keep
            py::arg( "loopRadiusStd").none(false) // keep
          )
      .def("specifyPrismaticLoops",
            &model::AtomDisplacementGenerator::specifyPrismaticLoops,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "slipSystemIDs").none(false),
            py::arg( "prismaticLoopRadii").none(false),
            py::arg( "prismaticLoopCenters").none(false),
            py::arg( "prismaticLoopSteps").none(false)
          )
      .def("specifyPrismaticLoopDensity",
            &model::AtomDisplacementGenerator::specifyPrismaticLoopDensity,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "prismaticLoopDensity").none(false),
            py::arg( "prismaticLoopRadiusMean").none(false),
            py::arg( "prismaticLoopRadiusStd").none(false),
            py::arg( "prismaticLoopStepMean").none(false),
            py::arg( "prismaticLoopStepStd").none(false)
          )
      .def("specifyPrismaticLoopDensitiesPerSlipSystem",
            &model::AtomDisplacementGenerator::specifyPrismaticLoopDensitiesPerSlipSystem,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "prismaticLoopDensitiesPerSlipSystem").none(false),
            py::arg( "prismaticLoopRadiusMean").none(false), // keep
            py::arg( "prismaticLoopRadiusStd").none(false), // keep
            py::arg( "prismaticLoopStepMean").none(false), // keep
            py::arg( "prismaticLoopStepStd").none(false) // keep
          )
      .def("specifyDipoleDensity",
            &model::AtomDisplacementGenerator::specifyDipoleDensity,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "dipoleDensity").none(false)
          )
      .def("generateMicrostructure", // uses microstructureSpecifications
         &model::AtomDisplacementGenerator::generateMicrostructure
         )
      .def("clearMicrostructureSpecifications", // uses microstructureSpecifications
         &model::AtomDisplacementGenerator::clearMicrostructureSpecifications
         )
      .def("regenerateMicrostructure", // reads inputFiles/initialMicrostructure.txt
            &model::AtomDisplacementGenerator::regenerateMicrostructure
          )
      .def("getCurrentStep",
            &model::AtomDisplacementGenerator::getCurrentStep
          )
      .def("setCurrentStep",
            &model::AtomDisplacementGenerator::setCurrentStep,
            py::arg("currentStep").none(false)
          )
      .def("readDefectiveCrystal",
            &model::AtomDisplacementGenerator::readDefectiveCrystal
          )
      .def("setOutputPath",
            &model::AtomDisplacementGenerator::setOutputPath,
            py::arg( "outputPath").none(false)
          )
      .def("readLammpsBox",
            &model::AtomDisplacementGenerator::readLammpsBox,
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
