/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ddpy_cpp_
#define model_ddpy_cpp_

#include <ddpy.h>

std::map<std::tuple< size_t, size_t>, double>
ddpy::DDInterface::getResolvedShearStresses()
{
   //Eigen::Matrix< double, 3, 3> stress;
   MatrixDim stress;
   std::map<std::tuple<size_t,size_t>,double> rss;
   if ( DC == nullptr)
   {
      std::cout << "error: getResolvedShearStresses(), "
        << " DefectiveCrystal not yet initialized" << std::endl;
      return rss;
   }
   if ( DC->externalLoadController != nullptr)
   {
      stress = DC->externalLoadController->stress( VectorDim::Zero());
   }
   else
   {
      std::cout << "error: getResolvedShearStresses(), "
        << " DC->externalLoadController isn't instantiated" << std::endl;
      return rss;
   }

   for ( const auto& grain : DC->poly.grains)
   {
      rss.clear();
      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
      { // loop over slip system
         std::tuple<size_t,size_t> key( grain.second.grainID, ss->sID);
         rss[ key] = (stress * (ss->unitNormal)).dot(ss->unitSlip);
      } // loop over slip system
   }
   return rss;
}

std::map<std::tuple< size_t, size_t>, double>
ddpy::DDInterface::getResolvedShearStrains()
{
   MatrixDim strain;
   std::map<std::tuple< size_t, size_t>, double> rss;
   if ( DC == nullptr)
   {
      std::cout << "error: getResolvedShearStrains(), "
        << " DefectiveCrystal not yet initialized" << std::endl;
      return rss;
   }
   if ( DC->externalLoadController != nullptr)
   {
      strain = DC->externalLoadController->strain( VectorDim::Zero());
   }
   else
   {
      std::cout << "error: getResolvedShearStrains(), "
        << " DC->externalLoadController == nullptr" << std::endl;
      return rss;
   }
   for ( const auto& grain : DC->poly.grains)
   {
      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
      { // loop over slip system
         std::tuple<size_t,size_t> key( grain.second.grainID, ss->sID);
         rss[ key] = (strain * (ss->unitNormal)).dot(ss->unitSlip);
      } // loop over slip system
   }
   return rss;
}

std::map<std::tuple<size_t,size_t>, double>
ddpy::DDInterface::getDensityPerSlipSystem()
{
   std::map<std::tuple<size_t,size_t>, double> densityPerSlipSystem;
   //std::map<size_t, double> densityPerSlipSystem;
   //densityPerSlipSystem[1] = 2.0; // debug
   //densityPerSlipSystem[3] = 0.1; // debug

   size_t ssID;
   size_t grainID;
   for ( const auto& loop : DC->DN->loops())
   {
      if ( loop.second.lock()->slipSystem() != nullptr)
      {
         grainID = loop.second.lock()->grain.grainID;
         ssID = loop.second.lock()->slipSystem()->sID;
         std::tuple<size_t,size_t> key( grainID, ssID);
         for ( const auto& loopLink : loop.second.lock()->loopLinks())
         {
            if ( loopLink->networkLink())
            {
               if ( ! loopLink->networkLink()->hasZeroBurgers())
               {
                  if (( !loopLink->networkLink()->isBoundarySegment())
                        &&(!loopLink->networkLink()->isGrainBoundarySegment())
                     )
                  {
                     // if ssID isn't in the map yet, instantiate w/value 0
                     densityPerSlipSystem.try_emplace( key, 0.0);
                     densityPerSlipSystem[ key ]
                        += loopLink->networkLink()->chord().norm()
                           /(
                            loopLink->networkLink()->loopLinks().size()
                            * ddBase->mesh.volume()
                            * std::pow( ddBase->poly.b_SI, 2)
                            // mesh.volume lacks |b|^3
                            // chord().norm() probably lacks |b| factor
                            // leaving two |b| factors in denominator
                            );
                  }
               }
            }
         }
      }
      else
      {
         //std::cout << "error: loop.second.lock()->slipSystem(): "
         //   << loop.second.lock()->slipSystem() << std::endl;
      }
   }
   return densityPerSlipSystem;
}

//std::list<std::tuple< size_t, size_t, double>>
std::map<std::tuple< size_t, size_t>, double>
ddpy::DDInterface::getPlasticStrains()
{
   //std::list<std::tuple< size_t, size_t, double>> plasticStrains;
   std::map<std::tuple< size_t, size_t>, double> plasticStrains;
   if ( DC == nullptr)
   {
      std::cout << "error: getPlasticStrains(), "
        << " DefectiveCrystal not yet initialized" << std::endl;
      return plasticStrains;
   }
   std::map<std::pair<int,int>,double> sspd(
         DC->DN->slipSystemPlasticDistortion());
   // [[grain ID, slip system ID], strain value in that slip system]

   //std::cout << "sspd.size: " << sspd.size() << std::endl;

   // copy into list of tuples to be returned
   for (const auto& itr : sspd)
   {
       plasticStrains[
             std::tuple( itr.first.first, itr.first.second)
         ] = itr.second;
             //.push_back(std::tuple( itr.first.first, itr.first.second, itr.second));
   }

   //for ( const auto& itr : plasticStrains)
   //{
   //   std::cout << "(grain,slipSystem,strain): " << "("
   //      << std::get<0>(itr) << ","
   //      << std::get<1>(itr) << ","
   //      << std::get<2>(itr) << ")" << std::endl;
   //}

   return plasticStrains;
}

//std::map<std::pair< size_t, size_t>, ddpy::DDInterface::VectorDim>
//   ddpy::DDInterface::getSlipSystemNormals() const
//{
//   size_t grainCount = 0;
//   size_t slipSystemCount = 0;
//   //std::map< std::pair<size_t,size_t>, model::ReciprocalLatticeVector<3>> normals;
//   std::map< std::pair<size_t,size_t>, VectorDim> normals;
//   //std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
//      normals;
//   if ( DC == NULL)
//   {
//      std::cout << "error: getSlipSystemNormals(), "
//        << " DefectiveCrystal not yet initialized" << std::endl;
//      return normals;
//   }
//   //std::map< std::pair<size_t,size_t>, std::string> normals;
//   // ((grain number, slip system number), plane normal)
//   for ( const auto& grain : DC->poly.grains)
//   {
//      std::cout << "grain " << grainCount << std::endl;
//      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
//      { // loop over slip system
//         //std::cout << "slip system " << slipSystemCount
//         //   << " plane normal:" << std::endl
//         //   << " (" << std::endl << ss->unitNormal
//         //   << ")" << std::endl;
//         normals.push_back(
//               std::tuple( grainCount, slipSystemCount, ss->unitNormal)
//               );
//         ++slipSystemCount;
//      }
//      ++grainCount;
//   }
//   for ( const auto& nn : normals)
//   {
//      std::cout << "grain " << std::get<0>( nn)
//         << " slip system " << std::get<1>( nn)
//         << std::endl << std::get<2>( nn) << std::endl;
//   }
//   return normals;
//}

//std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
//std::map< std::pair<size_t,size_t>, Eigen::Matrix<double,3,1> >
//   ddpy::DDInterface::getSlipSystemBurgersVectors() const
//{
//   size_t grainCount = 0;
//   size_t slipSystemCount = 0;
//   //std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
//   std::map< std::pair<size_t,size_t>, Eigen::Matrix<double,3,1>>
//      burgersVectors;
//   if ( DC == NULL)
//   {
//      std::cout << "error: getSlipSystemBurgersVectors(), "
//        << " DefectiveCrystal not yet initialized" << std::endl;
//      return burgersVectors;
//   }
//   for ( const auto& grain : DC->poly.grains)
//   {
//      ++grainCount;
//      std::cout << "grain " << grainCount << std::endl;
//      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
//      { // loop over slip system
//         ++slipSystemCount;
//         burgersVectors.emplace(
//               std::make_pair(
//                  std::make_pair( grainCount, slipSystemCount),
//                  ss->unitSlip)
//               );
//         //std::cout << "slip system " << slipSystemCount
//         //   << " burgers vector:" << std::endl
//         //   << " (" << std::endl << ss->unitSlip
//         //   << ")" << std::endl;
//      }
//   }
//   return burgersVectors;
//}

double ddpy::DDInterface::getBurgersMagnitude()
{
   if ( DC == nullptr)
   {
      std::cout << "error, getBurgersMagnitude: "
         << " DefectiveCrystal not yet instantiated." << std::endl;
      return -1;
   }
   if ( DC->DN == nullptr)
   {
      std::cout << "error, getBurgersMagnitude: "
         << " DefectiveCrystal::DefectiveNetwork not yet instantiated."
         << std::endl;
      return -1;
   }
   return DC->DN->poly.b_SI/1e-10;
}

void ddpy::DDInterface::readBurgersMagnitude(
      const std::string& materialPath)
{
   burgersMagnitude = model::TextFileParser( materialPath).readScalar<double>("b_SI",true);
   std::cout << "burgersMagnitude: " << burgersMagnitude << std::endl; // debug
   return;
}


size_t ddpy::DDInterface::getCurrentStep()
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot getCurrentStep until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return 0;
   }
   return DC->simulationParameters.runID;
}

void ddpy::DDInterface::setCurrentStep( const long int& step)
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot setCurrentStep until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }
   std::cout << "calling ddpy::DDInterface::setCurrentStep(" // debug
      << step << ")" << std::endl; // debug

   DC->simulationParameters.runID = step;
   model::DDconfigIO<3> evl(DC->simulationParameters.traitsIO.evlFolder);
   evl.read( DC->simulationParameters.runID);
   DC->DN->setConfiguration(evl);
   DC->externalLoadController = DC->getExternalLoadController(
            DC->simulationParameters, *DC, step
            );
   DC->simulationParameters.manageRestart();
   return;
}

//template <>
//void ddpy::DDInterface::setExternalLoad<model::UniformExternalLoadController<typename DefectiveCrystalType>>(
void ddpy::DDInterface::setExternalLoad(
      std::optional< const pybind11::array_t<double,
         pybind11::array::c_style | pybind11::array::forcecast>>&
            ExternalStress0In, // 3x3 matrix [Pa]
      std::optional< const pybind11::array_t<double,
         pybind11::array::c_style | pybind11::array::forcecast>>&
            ExternalStressRateIn, // 3x3 matrix [Pa/s]
      std::optional< const pybind11::array_t<double,
         pybind11::array::c_style | pybind11::array::forcecast>>&
            ExternalStrain0In, // 3x3 matrix [%]
      std::optional< const pybind11::array_t<double,
         pybind11::array::c_style | pybind11::array::forcecast>>&
            ExternalStrainRateIn, // 3x3 matrix [s^-1]
      std::optional< const pybind11::array_t<double,
         pybind11::array::c_style | pybind11::array::forcecast>>&
            MachineStiffnessRatioIn // Voigt format 11 22 33 12 23 13
      )
{
   if ( ddBase == nullptr) readddBase();
   if ( DC == nullptr)
   {
      std::cout << "error: cannot assign external load until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }

   if ( DC->externalLoadController == nullptr)
   {
      // instantiate a load controller by reading a file
      //  (e.g. uniformExternalLoadController.txt)
      DC->externalLoadController = DC->getExternalLoadController(
              DC->simulationParameters,
              *DC,
              DC->simulationParameters.runID
              );
   }

   if ( DC->simulationParameters.externalLoadControllerName
         != "UniformExternalLoadController")
   {
      std::cout << "error ddpy::DDInterface::setExternalLoad(): "
        << "DC->simulationParameters.externalLoadControllerName: "
        << DC->simulationParameters.externalLoadControllerName
        << " != UniformExternalLoadController ."
        << "Only UniformExternalLoadController is supported by "
        << "ddpy::DDInterface" << std::endl;
      return;
   }

   // Adapting calculations from UniformExternalLoadController constructor

   //TODO: convert values from SI to modelib's unitless units
   // read and assign ExternalStressRate if it was specified
   if ( ExternalStressRateIn.has_value())
   {
      auto ExternalStressRateInBuf = ExternalStressRateIn.value().unchecked<2>();
      for ( ssize_t ii=0; ii < ExternalStressRateIn.value().shape()[0]; ++ii)
         for ( ssize_t jj=0; jj < ExternalStressRateIn.value().shape()[1]; ++jj)
         {
            DC->externalLoadController->ExternalStressRate( ii, jj)
               = ExternalStressRateInBuf( ii, jj) // [Pa/s]
                  * DC->poly.b_SI/DC->poly.cs_SI  // [m/(m/s)]
                  / DC->poly.mu_SI; // [Pa]
         }
      assert((
         DC->externalLoadController->ExternalStressRate
            - DC->externalLoadController->ExternalStressRate.transpose()
             ).norm()<DBL_EPSILON && "ExternalStressRate is not symmetric."
            );
   }

   // read and assign ExternalStrainRate if it was specified
   if ( ExternalStrainRateIn.has_value())
   {
      auto ExternalStrainRateInBuf = ExternalStrainRateIn.value().unchecked<2>();
      for ( ssize_t ii=0; ii < ExternalStrainRateIn.value().shape()[0]; ++ii)
         for ( ssize_t jj=0; jj < ExternalStrainRateIn.value().shape()[1]; ++jj)
         {
            DC->externalLoadController->ExternalStrainRate( ii, jj)
               = ExternalStrainRateInBuf( ii, jj) // [s^-1]
                  * DC->poly.b_SI/DC->poly.cs_SI;  // [m/(m/s)]
         }
      assert((
               DC->externalLoadController->ExternalStrainRate
                - DC->externalLoadController->ExternalStrainRate.transpose()
              ).norm()
              < DBL_EPSILON && "ExternalStrainRate is not symmetric."
            );
   }

   // read and assign MachineStiffnessRatio if it was specified
   if ( MachineStiffnessRatioIn.has_value())
   {
      auto MachineStiffnessRatioInBuf = MachineStiffnessRatioIn.value().unchecked<1>();
      for ( ssize_t ii=0; ii < MachineStiffnessRatioIn.value().shape()[0]; ++ii)
      {
         DC->externalLoadController->MachineStiffnessRatio( ii)
               = MachineStiffnessRatioInBuf( ii);
      }
      DC->externalLoadController->MachineStiffnessRatio.block(0,0,1,DefectiveCrystalType::dim)
         =DC->externalLoadController->MachineStiffnessRatio.block(0,0,1,DefectiveCrystalType::dim)
            * (2+DC->externalLoadController->lambda);

      static constexpr int voigtSize = DefectiveCrystalType::dim*(DefectiveCrystalType::dim+1)/2;
      Eigen::Matrix<double,voigtSize,voigtSize>  Cinv=Eigen::Matrix<double,voigtSize,voigtSize>::Identity();
      double nu( DC->DN->poly.nu);
      double nu_use( DC->externalLoadController->nu_use);
      Cinv.block(0,0,DefectiveCrystalType::dim,DefectiveCrystalType::dim)
         <<
             (nu_use)/nu, -(nu_use),       -(nu_use),
            -(nu_use),   (nu_use)/nu,     -(nu_use),
            -(nu_use),   -(nu_use),       (nu_use)/nu;

      Eigen::Matrix<double,voigtSize,voigtSize>  machinestiffness
         =
         DC->externalLoadController->MachineStiffnessRatio.asDiagonal();
      DC->externalLoadController->stressmultimachinestiffness
         = (Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse();
      DC->externalLoadController->strainmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse()*machinestiffness;
   }

   // ExternalStrain is either assigned the value of ExternalStrain0 or
   //  left untouched.
   if ( ExternalStrain0In.has_value())
   {
      auto ExternalStrain0InBuf = ExternalStrain0In.value().unchecked<2>();
      for ( ssize_t ii=0; ii < ExternalStrain0In.value().shape()[0]; ++ii)
         for ( ssize_t jj=0; jj < ExternalStrain0In.value().shape()[1]; ++jj)
         {
            DC->externalLoadController->ExternalStrain( ii, jj)
               = ExternalStrain0InBuf( ii, jj); // [%]
         }
      assert(
               ( DC->externalLoadController->ExternalStrain
                  - DC->externalLoadController->ExternalStrain.transpose()
               ).norm()
               < DBL_EPSILON && "ExternalStrain0 is not symmetric."
            );
   }

   // If ExternalStress0In is specified, then assign it to ExternalStress
   //  and calculate plasticStrain as the difference between ExternalStrain
   //  and elasticstrain(ExternalStress,nu_use).
   // The following mimicks how ExternalStress and plasticStrain are
   //  assigned when reading the state of the external load from F/F_0.txt.
   if ( ExternalStress0In.has_value())
   {
      auto ExternalStress0InBuf = ExternalStress0In.value().unchecked<2>();
      for ( ssize_t ii=0; ii < ExternalStress0In.value().shape()[0]; ++ii)
         for ( ssize_t jj=0; jj < ExternalStress0In.value().shape()[1]; ++jj)
         {
            DC->externalLoadController->ExternalStress( ii, jj) // unitless
               = ExternalStress0InBuf( ii, jj) / DC->poly.mu_SI;
         }
      assert(
               (
                  DC->externalLoadController->ExternalStress
                   - DC->externalLoadController->ExternalStress.transpose()
               ).norm()
               < DBL_EPSILON && "ExternalStress0 is not symmetric."
            );
      DC->externalLoadController->plasticStrain
         = DC->externalLoadController->ExternalStrain
          - DC->externalLoadController->elasticstrain(
                DC->externalLoadController->ExternalStress,
                DC->externalLoadController->nu_use
                ); // mimicking assignment when reading from F/F_0.txt
   }
   else
   {  // If ExternalStress0In is not specified
      //  then mimick assignments used by the load constructor when
      //  F/F_0.txt cannot be read.

      if (
            ( ExternalStrain0In.has_value())
            || ( MachineStiffnessRatioIn.has_value())
         )
      {  // If ExternalStrain or stiffness ratio were passed to this
         //  function, then the ExternalStress and plasticStrain would
         //  need to be updated.
         MatrixDim pdr( DC->DN->plasticDistortion());
         DC->externalLoadController->plasticStrain
            = ( pdr +pdr.transpose())*0.5;

         MatrixDim dstrain(
            DC->externalLoadController->ExternalStrain // possibly new
            - DC->externalLoadController->plasticStrain
            );

         MatrixDim S_stress( DC->externalLoadController->ExternalStress);
         DC->externalLoadController->ExternalStress
            = DC->externalLoadController->stressconsidermachinestiffness(
                  dstrain,
                  S_stress
                  );
      }
      // If none of ExternalStrain, stiffness ratio, or ExternalStress
      //  were passed to this function, the values previously assigned by
      //  the constructor should remain.
   }
   return;
}

void ddpy::DDInterface::setOutputFrequency( const long int& outputFrequency)
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot assign output frequency until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }
   if ( DC->DN == nullptr)
   {
      std::cout << "error: cannot assign output frequency until "
         << "DefectiveNetwork is instantiated" << std::endl;
      return;
   }

   int outputFrequencyInt = static_cast<int>( outputFrequency);

   if ( outputFrequencyInt <= 0)
   {
      std::cout << "error: static_cast<int>( outputFrequency) "
         << outputFrequencyInt << " <= 0" << std::endl;
      return;
   }

   DC->DN->outputFrequency = outputFrequencyInt;
   return;
}

void ddpy::DDInterface::setEndingStep( const long int& endingStep)
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot assign ending step until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }
   if ( endingStep <= DC->simulationParameters.Nsteps)
   {
      std::cout << "warning: assigning endingStep=" << endingStep
         << " <= currentStep:" << DC->simulationParameters.Nsteps
         << std::endl;
   }
   DC->simulationParameters.Nsteps = endingStep;
   return;
}

void ddpy::DDInterface::readddBase()
{
   //resetStaticIDs();
   std::cout << "resetting staticIDs" << std::endl; // debug
   if ( ddBase != nullptr) resetStaticIDs(); // TODO: is this necessary?
   std::cout << "reading " << dddFolderPath << std::endl; // debug
   ddBase =  std::unique_ptr<DislocationDynamicsBaseType>(
         new DislocationDynamicsBaseType( dddFolderPath)
         );
   return;
}

void ddpy::DDInterface::setBoxBounds(
      const double& xlo,
      const double& xhi,
      const double& ylo,
      const double& yhi,
      const double& zlo,
      const double& zhi
      )
{
   if ( burgersMagnitude <= 0)
   {
      std::cout << "error, setBoxBounds: burgersMagnitude <= 0\n"
         << "try calling readBurgersMagnitude( <material file path>) first"
         << std::endl;
      return;
   }
   boxBounds[0] = xlo * 1e-10 / burgersMagnitude;
   boxBounds[1] = xhi * 1e-10 / burgersMagnitude;
   boxBounds[2] = ylo * 1e-10 / burgersMagnitude;
   boxBounds[3] = yhi * 1e-10 / burgersMagnitude;
   boxBounds[4] = zlo * 1e-10 / burgersMagnitude;
   boxBounds[5] = zhi * 1e-10 / burgersMagnitude;
   boxVolume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo); // [\AA^{2}]

   return;
}

void ddpy::DDInterface::readDefectiveCrystal()
{
   //if ( DC != NULL) delete DC;
   if ( ddBase == nullptr) readddBase();
   DC = std::unique_ptr<DefectiveCrystalType>(
         new DefectiveCrystalType( *ddBase)
         );
   //DC->simulationParameters.manageRestart();
   return;
}

void ddpy::DDInterface::generateMicrostructure()
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


void ddpy::DDInterface::writeConfigToTxt()
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot write configuration to file until"
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }
   if ( DC->DN == nullptr)
   {
      std::cout << "error: cannot write configuration to file until "
         << "DefectiveNetwork is instantiated" << std::endl;
      return;
   }
   std::cout << "attempting to write a configuration" << std::endl; // debug
   DC->DN->io().output( DC->simulationParameters.runID);
   return;
}

void ddpy::DDInterface::runGlideSteps( const size_t& stepsToRun)
{
   if ( DC == nullptr) readDefectiveCrystal();
   if ( DC->DN == nullptr)
   {
      std::cout << "error: ddpy::DDInterface::runGlideSteps()"
         << " DefectiveCrystal::DefectiveNetwork not yet instantiated."
         << std::endl;
      return;
   }
   if ( DC->externalLoadController == nullptr)
   {
      std::cout << "error: ddpy::DDInterface::runGlideSteps()"
        << " DC->externalLoadController isn't instantiated" << std::endl;
      return;
   }

   if ( collectMechanicalMeasurements
         //collectPlasticStrainRates || collectResolvedStrainRates
         //|| collectDensities
      )
   {
      ssize_t initialStep = DC->simulationParameters.runID;
      ssize_t endStep = initialStep + stepsToRun;
      if ( DC->poly.grains.size() != 1)
      {
         std::cout << "error: ddpy::DDInterface::runGlideSteps()"
            << " is not yet able to collect measurements when"
            << " the number of grains is not 1."
            << " DC->poly.grains.size() "
            << DC->poly.grains.size()
            //<< ", collectPlasticStrain: " << collectPlasticStrain
            //<< ", collectResolvedStrainRates: " << collectResolvedStrainRates
            //<< ", collectDensities: " << collectDensities
            << std::endl;
         return;
      }
      const auto& grain = *(DC->poly.grains.begin());
      const int grainID = grain.second.grainID;

      size_t measurementNumber;
      if ( times == nullptr)
      {
         //times = new pybind11::array_t<double, pybind11::array::c_style>;
         //runIDs = new pybind11::array_t<ssize_t, pybind11::array::c_style>;
         times = std::make_shared< std::vector<double> >();
         runIDs = std::make_shared< std::vector<ssize_t> >();
      }
      measurementNumber = times->size();
      if ( runIDs->size() != measurementNumber)
      {
         std::cout << "error: ddpy::DDInterface::runGlideSteps()"
            << " pre-existing measurement data has differing sizes"
            << std::endl;
      }

      //times.emplace_back( prevTime);
      //runIDs.emplace_back( initialStep);
      //// DC->simulationParameters.totalTime - prevTime

      // allocate and initiate temporary storage for measured quantities
      std::map<std::pair<int,int>,double> currentSSPD;//key:grainID,ssID
      currentSSPD = DC->DN->slipSystemPlasticDistortion(); // unitless
      //MatrixDim strain
      //   = DC->externalLoadController->strain( VectorDim::Zero());
      MatrixDim stress // units: [DC->poly.mu_SI [Pa]]
         = DC->externalLoadController->stress( VectorDim::Zero())
            * DC->poly.mu_SI; // [Pa]

      //// ensure times, runIDs and measured quantity vectors have same size
      //if ( times.size() != runIDs.size())
      //{
      //   std::cout << "error, ddpy::DDInterface::runGlideSteps()"
      //      << ": times.size() " << times.size()
      //     << " != runIDs.size() " << runIDs.size() << std::endl;
      //   return;
      //}
      //for ( const auto& ss : grain.second.singleCrystal->slipSystems())
      //{
      //   //if ( times.size() != plasticStrainRates[ss->sID].size())
      //   //{
      //   //   std::cout << "error, ddpy::DDInterface::runGlideSteps()"
      //   //      << ": times.size() " << times.size()
      //   //      << " != plasticStrainRates[" << ss->sID << "].size() "
      //   //      << plasticStrainRates[ss->sID].size() << std::endl;
      //   //   return;
      //   //}
      //   //if ( times.size() != resolvedStrainRates[ss->sID].size())
      //   //{
      //   //   std::cout << "error, ddpy::DDInterface::runGlideSteps()"
      //   //      << ": times.size() " << times.size()
      //   //      << " != resolvedStrainRates[" << ss->sID << "].size() "
      //   //      << resolvedStrainRates[ss->sID].size() << std::endl;
      //   //   return;
      //   //}
      //   if ( times.size() != densityPerSlipSystem[ss->sID].size())
      //   {
      //      std::cout << "error, ddpy::DDInterface::runGlideSteps()"
      //         << ": times.size() " << times.size()
      //         << " != densityPerSlipSystem[" << ss->sID << "].size() "
      //         << densityPerSlipSystem[ss->sID].size() << std::endl;
      //      return;
      //   }
      //} // ensured times, runIDs and measured quantity vectors have same size

      // determine the required size of contiguous measurement vectors

      size_t sizeOfMeasurementVector
         = times->size()
            + static_cast<size_t>(
                  std::ceil(
                     static_cast<double>( stepsToRun)
                     / static_cast<double>( stepsBetweenMeasurements)
                  )+1 // allow for storing initial state
               );

      try
      {
         times->reserve( sizeOfMeasurementVector);
         runIDs->reserve( sizeOfMeasurementVector);
      }
      catch(const std::exception& e)
      {
         std::cout << "error: ddpy::DDInterface::runGlideSteps()"
            << " could not reserve vector<> storage for "
            << "times and runIDs. Maybe try decreasing the total number"
            << " of measurements by increasing"
            << " DDInterface::stepsBetweenMeasurements"
            << " or disabling measurements for a greater portion of steps."
            << std::endl;
         return;
      }

      //times->push_back( DC->simulationParameters.totalTime
      //      * (( DC->poly.b_SI)/( DC->poly.cs_SI))); // [sec]
      //runIDs->push_back( DC->simulationParameters.runID);
      ////times[ measurementNumber] = DC->simulationParameters.totalTime;
      ////runIDs[ measurementNumber] = DC->simulationParameters.runID;

      // allocate and initialize storage for densities per slip system
      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
      {
         // if ssID isn't in the map yet, instantiate first step w/0.0
         //  else, append a new element to the slip system's vector
         // instantiate and reserve densityPerSlipSystem
         if  (! densityPerSlipSystem.contains( ss->sID))
         {
            try
            {
               densityPerSlipSystem.try_emplace(
                  ss->sID, // key
                  std::make_shared< std::vector<double> >()//( sizeOfMeasurementVector,
                     //0.0) // skip initial configuration
                  );
               densityPerSlipSystem[ ss->sID]->reserve(
                     sizeOfMeasurementVector
                     );
            }
            catch( const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "densityPerSlipSystem[" << ss->sID << "]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }
         else
         {// densityPerSlipSystem[ ssID] already contains a vector<double>
            // append an appropriate number of steps to it
            try
            {
               densityPerSlipSystem[ ss->sID]->reserve(
                     sizeOfMeasurementVector
                     );
            }
            catch(const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "densityPerSlipSystem[" << ss->sID << "]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }
         //densityPerSlipSystem[ ss->sID]->push_back( 0.0);

         // instantiate and reserve glissileDensityPerSlipSystem
         if  (! glissileDensityPerSlipSystem.contains( ss->sID))
         {
            try
            {
               glissileDensityPerSlipSystem.try_emplace(
                  ss->sID, // key
                  std::make_shared< std::vector<double> >()//( sizeOfMeasurementVector,
                     //0.0) // skip initial configuration
                  );
               glissileDensityPerSlipSystem[ ss->sID]->reserve(
                     sizeOfMeasurementVector
                     );
            }
            catch( const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "glissileDensityPerSlipSystem[" << ss->sID << "]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }
         else
         {// glissileDensityPerSlipSystem[ ssID] already contains a vector<double>
            // append an appropriate number of steps to it
            try
            {
               glissileDensityPerSlipSystem[ ss->sID]->reserve(
                     sizeOfMeasurementVector
                     );
            }
            catch(const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "glissileDensityPerSlipSystem[" << ss->sID << "]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }
         //glissileDensityPerSlipSystem[ ss->sID]->push_back( 0.0);

         // instantiate and reserve sessileDensityPerSlipSystem
         if  (! sessileDensityPerSlipSystem.contains( ss->sID))
         {
            try
            {
               sessileDensityPerSlipSystem.try_emplace(
                  ss->sID, // key
                  std::make_shared< std::vector<double> >()//( sizeOfMeasurementVector,
                     //0.0) // skip initial configuration
                  );
               sessileDensityPerSlipSystem[ ss->sID]->reserve(
                     sizeOfMeasurementVector
                     );
            }
            catch( const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "sessileDensityPerSlipSystem[" << ss->sID << "]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }
         else
         {// sessileDensityPerSlipSystem[ ssID] already contains a vector<double>
            // append an appropriate number of steps to it
            try
            {
               sessileDensityPerSlipSystem[ ss->sID]->reserve(
                     sizeOfMeasurementVector
                     );
            }
            catch(const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "sessileDensityPerSlipSystem[" << ss->sID << "]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }
         //sessileDensityPerSlipSystem[ ss->sID]->push_back( 0.0);
         // density variables have had storage reserved

         // if ssID isn't in the map yet, instantiate first step w/0.0
         //  else, append a new element to the slip system's vector
         if  (! slipSystemPlasticDistortion.contains( ss->sID))
         {
            try
            {
               slipSystemPlasticDistortion.try_emplace(
                  ss->sID, // key
                  std::make_shared< std::vector<double> >()
                  );
               slipSystemPlasticDistortion[ ss->sID]->reserve(
                     sizeOfMeasurementVector
                     );
            }
            catch( const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "slipSystemPlasticDistortion[" << ss->sID << "]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }
         else
         {// slipSystemPlasticDistortion[ ssID] already contains a vector<double>
            // append an appropriate number of steps to it
            try
            {
               slipSystemPlasticDistortion[ ss->sID]->reserve(
                     sizeOfMeasurementVector);
            }
            catch(const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "slipSystemPlasticDistortion[" << ss->sID << "]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }

         //slipSystemPlasticDistortion[ ss->sID]->push_back(
         //   currentSSPD[ std::pair<int,int>( grainID, ss->sID)] //unitless
         //   );

      } // for ( const auto& ss : grain.second.singleCrystal->slipSystems())

      //// initialize densityPerSlipSystem
      //for ( const auto& loop : DC->DN->loops())
      //{
      //   if ( loop.second.lock()->slipSystem() != nullptr)
      //   {
      //      for ( const auto& loopLink : loop.second.lock()->loopLinks())
      //      {
      //         if ( loopLink->networkLink())
      //         {
      //            if ( ! loopLink->networkLink()->hasZeroBurgers())
      //            {
      //               if (( !loopLink->networkLink()->isBoundarySegment())
      //                     &&(!loopLink->networkLink()->isGrainBoundarySegment())
      //                  )
      //               {
      //                  double densityContribution(
      //                        loopLink->networkLink()->chord().norm()
      //                        /(
      //                         loopLink->networkLink()->loopLinks().size()
      //                         * ddBase->mesh.volume() // [b^3]
      //                         * std::pow( ddBase->poly.b_SI, 2) //[m^2/b^2]
      //                         ) // [m^-2];
      //                        );

      //                  densityPerSlipSystem[
      //                     loop.second.lock()->slipSystem()->sID
      //                  ]->back() += densityContribution;

      //                  if ( loopLink->networkLink()->isGlissile())
      //                  {
      //                     glissileDensityPerSlipSystem[
      //                        loop.second.lock()->slipSystem()->sID
      //                     ]->back() += densityContribution;
      //                  }

      //                  if ( loopLink->networkLink()->isSessile())
      //                  {
      //                     sessileDensityPerSlipSystem[
      //                        loop.second.lock()->slipSystem()->sID
      //                     ]->back() += densityContribution;
      //                  }
      //               }
      //            }
      //         }
      //      }
      //   }
      //   else
      //   {
      //         std::cout
      //            << "warning: loop.second.lock()->slipSystem() == nullptr"
      //            << "loop.second.lock()->tag(): "
      //            << loop.second.lock()->tag()<< std::endl
      //            << "loop.second.lock()->glidePlane: "
      //            << loop.second.lock()->glidePlane << std::endl
      //            << "loop.second.lock()->slippedArea(): "
      //            << loop.second.lock()->slippedArea() << std::endl
      //            << std::endl;
      //   }
      //} // accumulated density per slip system, initial density skipped

      // allocate stressTensorComponents // and strainTensorComponents
      for ( const auto& key : tensorComponentKeys)
      {
         if ( ! stressTensorComponents.contains( key))
         {
            try
            {
               stressTensorComponents.try_emplace(
                     key,
                     std::make_shared< std::vector<double> >()// sizeOfMeasurementVector)
                     );
               stressTensorComponents[ key]->reserve(
                     sizeOfMeasurementVector
                     );
            }
            catch( const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "stressTensorComponents[("
                  << key.first << "," << key.second << ")]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }
         else
         {// stressTensorComponents[key] already contains a vector<double>
            // append an appropriate number of steps to it
            try
            {
               stressTensorComponents[ key]->reserve(
                     sizeOfMeasurementVector
                     );
            }
            catch(const std::exception& e)
            {
               std::cout << "error: ddpy::DDInterface::runGlideSteps()"
                  << " could not reserve vector<> storage for "
                  << "stressTensorComponents[("
                  << key.first << "," << key.second << ")]"
                  << ". Maybe try increasing "
                  << "DDInterface::stepsBetweenMeasurements."
                  << std::endl;
               return;
            }
         }

         //if ( ! strainTensorComponents.contains( key))
         //{
         //   try
         //   {
         //      strainTensorComponents.try_emplace(
         //            key,
         //            std::vector<double>()// sizeOfMeasurementVector)
         //            );
         //      strainTensorComponents[ key].reserve(
         //            sizeOfMeasurementVector
         //            );
         //   }
         //   catch( const std::exception& e)
         //   {
         //      std::cout << "error: ddpy::DDInterface::runGlideSteps()"
         //         << " could not reserve vector<> storage for "
         //         << "strainTensorComponents[("
         //         << key.first << "," << key.second << ")]"
         //         << ". Maybe try increasing "
         //         << "DDInterface::stepsBetweenMeasurements."
         //         << std::endl;
         //      return;
         //   }
         //}
         //else
         //{// strainTensorComponents[key] already contains a vector<double>
         //   // append an appropriate number of steps to it
         //   try
         //   {
         //      strainTensorComponents[ key].reserve(
         //            sizeOfMeasurementVector
         //            );
         //   }
         //   catch(const std::exception& e)
         //   {
         //      std::cout << "error: ddpy::DDInterface::runGlideSteps()"
         //         << " could not reserve vector<> storage for "
         //         << "strainTensorComponents[("
         //         << key.first << "," << key.second << ")]"
         //         << ". Maybe try increasing "
         //         << "DDInterface::stepsBetweenMeasurements."
         //         << std::endl;
         //      return;
         //   }
         //}

         //// collect initial stress and strain components
         //stressTensorComponents[ key]->push_back( //[ measurementNumber]
         //   stress( key.first-1, key.second-1) // index from 0, not 1
         //   );
         ////strainTensorComponents[ key].push_back( //[ measurementNumber]
         ////   strain( key.first, key.second)
         ////   );
      } // for ( const auto& key : tensorComponentKeys)

      while ( DC->simulationParameters.runID < endStep)
      {
         //measurementNumber += 1;
         // set end step
         DC->simulationParameters.Nsteps
            = DC->simulationParameters.runID + stepsBetweenMeasurements;

         // run DDD steps until reaching the specified end step
         DC->runGlideSteps();
         DC->DN->updateGeometry();
         // ??updateLoadControllers(simulationParamerters.runID, false);
         // measurements are collected after each DDD time step, but
         //  runGlideSteps() leaves the data in need of updating, hence
         //  the need to run updateGeometry()

         // append measurements to:
         //  std::vector<double> times
         //  std::vector<double> runIDs
         //  std::map< size_t, // slip system,
         //     std::vector<double> // time series of data
         //      slipSystemPlasticDistortion
         //  std::map<size_t, std::vector< double> > densityPerSlipSystem;
         //  std::map<size_t, std::vector< double> > stressTensorComponents; //keys:11,22,33,12,13,23
         //  std::map<size_t, std::vector< double> > strainTensorComponents; //keys:11,22,33,12,13,23
         // std::vector<std::pair<size_t,size_t> > tensorComponentKeys({11,22,33,12,13,23});

         times->push_back( DC->simulationParameters.totalTime
                  * ( DC->poly.b_SI)/(DC->poly.cs_SI)
               );
         runIDs->push_back( DC->simulationParameters.runID);
         //times[ measurementNumber] = DC->simulationParameters.totalTime;
         //runIDs[ measurementNumber] = DC->simulationParameters.runID;

         // retrieve stress and strain tensors and store their components
         stress = DC->externalLoadController->stress( VectorDim::Zero())
                  * DC->poly.mu_SI; // [Pa]

         //strain = DC->externalLoadController->strain( VectorDim::Zero());
         // store stressTensorComponents
         for ( const auto& key : tensorComponentKeys)//11,22,33,12,13,23
         {
            //// check vector bounds
            //if ( measurementNumber >= stressTensorComponents[ key].size())
            //{
            //   std::cout << "error, ddpy::DDInterface::runGlideSteps()"
            //      << " measurementNumber " << measurementNumber
            //      << " >= stressTensorComponents[("
            //      << key.first << "," << key.second << ")].size() "
            //      << stressTensorComponents[ ss->sID].size()
            //   return;
            //}
            //if ( measurementNumber >= strainTensorComponents[ key].size())
            //{
            //   std::cout << "error, ddpy::DDInterface::runGlideSteps()"
            //      << " measurementNumber " << measurementNumber
            //      << " >= strainTensorComponents[("
            //      << key.first << "," << key.second << ")].size() "
            //      << strainTensorComponents[ ss->sID].size()
            //   return;
            //}

            // TODO: consider evaluating and storing resolved shear
            //        stress and strain, since it will be done anyway
            //        by getMechanicalMeasurements(), and it would be
            //        more efficient to return those results as pointers

            // assign values
            stressTensorComponents[ key]->push_back(
               stress( key.first-1, key.second-1) // [Pa]
               );
            //strainTensorComponents[ key].push_back(
            //   strain( key.first, key.second)
            //   );
            //stressTensorComponents[ key][ measurementNumber]
            //   = stress( key.first, key.second);
            //strainTensorComponents[ key][ measurementNumber]
            //   = strain( key.first, key.second);
         }

         // std::map<std::pair<int,int>,double> currentSSPD;//key:grainID,ssID
         currentSSPD = DC->DN->slipSystemPlasticDistortion(); // unitless

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
         for ( const auto& ss : grain.second.singleCrystal->slipSystems())
         {
            //// check vector bounds
            //if ( measurementNumber >= densityPerSlipSystem[ ss->sID].size())
            //{
            //   std::cout << "error, ddpy::DDInterface::runGlideSteps()"
            //      << " measurementNumber " << measurementNumber
            //      << " >= densityPerSlipSystem["
            //      << ss->sID << "].size() "
            //      << densityPerSlipSystem[ ss->sID].size()
            //   return;
            //}

            //if ( measurementNumber >= slipSystemPlasticDistortion[ ss->sID].size())
            //{
            //   std::cout << "error, ddpy::DDInterface::runGlideSteps()"
            //      << " measurementNumber " << measurementNumber
            //      << " >= slipSystemPlasticDistortion["
            //      << ss->sID << "].size() "
            //      << slipSystemPlasticDistortion[ ss->sID].size()
            //   return;
            //}

            // assign initial values to slipSystemPlasticDistortion
            slipSystemPlasticDistortion[ ss->sID]->push_back(
               currentSSPD[ std::pair<int,int>( grainID, ss->sID)]
               ); // unitless

            //slipSystemPlasticDistortion[ ss->sID][ measurementNumber]
            //   = currentSSPD[ std::pair<int,int>( grainID, ss->sID)];

            densityPerSlipSystem[ ss->sID]->push_back( 0.0); // [m^{-2}]
            glissileDensityPerSlipSystem[ ss->sID]->push_back( 0.0); // [m^{-2}]
            sessileDensityPerSlipSystem[ ss->sID]->push_back( 0.0); // [m^{-2}]
         } // for ( const auto& ss : grain.second.singleCrystal->slipSystems())

         // accumulate density per slip system, initial density skipped
         for ( const auto& loop : DC->DN->loops())
         {
            if ( loop.second.lock()->slipSystem() != nullptr)
            {
               for ( const auto& loopLink : loop.second.lock()->loopLinks())
               {
                  if ( loopLink->networkLink())
                  {
                     if ( ! loopLink->networkLink()->hasZeroBurgers())
                     {
                        if (( !loopLink->networkLink()->isBoundarySegment())
                              &&(!loopLink->networkLink()->isGrainBoundarySegment())
                           )
                        {
                           double densityContribution(
                                 loopLink->networkLink()->chord().norm()
                                 /(
                                  loopLink->networkLink()->loopLinks().size()
                                  * ddBase->mesh.volume() // [b^3]
                                  * std::pow( ddBase->poly.b_SI, 2) //[m^2/b^2]
                                  ) // [m^-2];
                                 );

                           densityPerSlipSystem[
                              loop.second.lock()->slipSystem()->sID
                           ]->back() += densityContribution;

                           if ( loopLink->networkLink()->isGlissile())
                           {
                              glissileDensityPerSlipSystem[
                                 loop.second.lock()->slipSystem()->sID
                              ]->back() += densityContribution;
                           }

                           if ( loopLink->networkLink()->isSessile())
                           {
                              sessileDensityPerSlipSystem[
                                 loop.second.lock()->slipSystem()->sID
                              ]->back() += densityContribution;
                           }
                        }
                     }
                  }
               }
            }
         } // accumulated density per slip system, initial density skipped
      }
   }
   else
   { // if no measurements were enabled, run without interruption
      // set step to end on
      DC->simulationParameters.Nsteps
         = DC->simulationParameters.runID + stepsToRun;
      DC->runGlideSteps();
   }
   return;
}

std::tuple<
   //std::vector<double>, // time
   pybind11::array_t<double, pybind11::array::c_style>, // time
   //std::vector<size_t>, // runIDs
   pybind11::array_t<ssize_t, pybind11::array::c_style>, // runIDs
   std::map< // measurementsMap
      std::string, // property name
      std::map<
         ssize_t, // slip system ID or tensor component index
         pybind11::array_t<double, pybind11::array::c_style>
         //std::vector<double>  // time series data
         >
      >
   >
ddpy::DDInterface::getMechanicalMeasurements()
{
   std::map<
      std::string,
      std::map< ssize_t,
         //std::shared_ptr<
            pybind11::array_t<double, pybind11::array::c_style>
            //std::vector<double>
         //   >
         >
      >
      measurements;

   if ( times == nullptr)
   {
      times = std::make_shared< std::vector<double> >();
   }
   if ( runIDs == nullptr)
   {
      runIDs = std::make_shared< std::vector<ssize_t> >();
   }

   pybind11::array_t<double> timesPy( pybind11::cast( *times)); // TODO: this fails
   pybind11::array_t<ssize_t> runIDsPy( pybind11::cast( *runIDs));

   if ( DC == nullptr)
   {
      std::cout << "error: "
         << "ddpy::DDInterface::getMechanicalMeasurements()"
         << " called while DC == nullptr" << std::endl;
      return std::make_tuple(
            timesPy,
            runIDsPy,
            measurements // map<string,map<size_t,shared_ptr<vector<double>>>>
            );
   }

   measurements.try_emplace(
         "density",
         std::map<
               ssize_t,
               pybind11::array_t< double, pybind11::array::c_style>
            >()
         );
   for ( const auto& [key, value] : densityPerSlipSystem)
   {
      //std::cout << "densityPerSlipSystem[" << key << "]:" ;
      //for (size_t ii=0; ii< value->size(); ++ii)
      //   std::cout << (*value)[ii] << ",";
      //std::cout << "moving density of key " << key
      //   << " with [0]: "
      //   << (*value)[0] << std::endl; // debug

      // TODO: fix. This is still copying the data, rather than passing its ownership to Python.
      measurements["density"][ key] = pybind11::cast( *value);
   }

   for ( const auto& [key, value] : glissileDensityPerSlipSystem)
   {
      // TODO: fix. This is still copying the data, rather than passing its ownership to Python.
      measurements["glissileDensity"][ key] = pybind11::cast( *value);
   }

   for ( const auto& [key, value] : sessileDensityPerSlipSystem)
   {
      // TODO: fix. This is still copying the data, rather than passing its ownership to Python.
      measurements["sessileDensity"][ key] = pybind11::cast( *value);
   }

   measurements.try_emplace(
         "slipSystemPlasticDistortion",
         std::map<
               ssize_t,
               pybind11::array_t< double, pybind11::array::c_style>
            >()
         );
   for ( const auto& [key, value] : slipSystemPlasticDistortion)
   {
      measurements["slipSystemPlasticDistortion"][ key]
         = pybind11::cast( *value);
   }

   measurements.try_emplace(
         "stressTensorComponent",
         std::map<
               ssize_t,
               pybind11::array_t< double, pybind11::array::c_style>
            >()
         );
   for ( const auto& key : tensorComponentKeys)
   {
      ssize_t outputKey( key.first * 10 + key.second);
      measurements[ "stressTensorComponent"].try_emplace(
            outputKey,
            pybind11::cast( *(stressTensorComponents[ key])) // [Pa]
            );
   }

   const auto& grain = *(DC->poly.grains.begin());
   MatrixDim tmpTensor; // for resolving tensors onto slip systems
   std::shared_ptr< std::vector< double> > rss;
   // determine quantities indexed by slip system
   for ( const auto& ss : grain.second.singleCrystal->slipSystems())
   {
      // resolvedStrain
      // evaluate resolved strain for the current slip system at time[ ii]

      measurements.try_emplace(
            "resolvedShearStress",
            std::map<
                  ssize_t,
                  pybind11::array_t< double, pybind11::array::c_style>
               >()
            );
      //measurements["resolvedShearStress"].try_emplace(
      //      ss->sID,
      //      //std::make_shared< std::vector<double> >( times->size())
      //      );
      ////measurements["resolvedShearStress"][ ss->sID]->reserve( times->size());
      rss = std::make_shared< std::vector< double> >();// times->size());
      rss->reserve( times->size());

      // iterate over time, evaluate resolved shear stress and strain,
      //  and assign them to output variables
      std::pair< size_t, size_t> tmpKey( *(tensorComponentKeys.begin()));
      for ( size_t tt=0; tt < stressTensorComponents[ tmpKey]->size(); ++tt)
      {
         // assign tensor components of current time step to a tensor
         for ( const auto& key : tensorComponentKeys)
         {
            tmpTensor( key.first-1, key.second-1)
               = (*(stressTensorComponents[ key]))[ tt]; // [Pa]
            if ( key.second != key.first)
            {
               tmpTensor( key.second-1, key.first-1)
                  = tmpTensor( key.first-1, key.second-1);
            }
         }
         // evaluate resolved shear stress
         //(*(measurements["resolvedShearStress"][ ss->sID]))[ tt]
         rss->push_back(
               (
                  tmpTensor * ( // [Pa]
                     ss->unitNormal // [b] or unitless ?
                     //   * ddBase->poly.b_SI // [m/b]
                     ) //[Pa m] or [Pa] ?
               ).dot( ss->unitSlip // [b] or unitless ?
                  //* ( ddBase->poly.b_SI // [m/b]
                  //   )
                  )); // [Pa]
      }
      measurements[ "resolvedShearStress"][ ss->sID]
         = pybind11::cast( *rss);
   } // for ( const auto& ss : grain.second.singleCrystal->slipSystems())
   return std::make_tuple(
         timesPy,
         runIDsPy,
         measurements // map<string,map<size_t,vector<double>>>
         );
}

void ddpy::DDInterface::clearMechanicalMeasurements()
{
   times->clear();
   times->shrink_to_fit();
   runIDs->clear();
   runIDs->shrink_to_fit();

   if ( DC == nullptr)
   {
      std::cout << "warning: called "
         << "ddpy::DDInterface::clearMechanicalMeasurements()"
         << " while DC == nullptr" << std::endl;
      return;
   }

   const auto& grain = *(DC->poly.grains.begin());

   // measurements indexed by slip system IDs
   for ( const auto& ss : grain.second.singleCrystal->slipSystems())
   {
      slipSystemPlasticDistortion[ ss->sID]->clear();
      slipSystemPlasticDistortion[ ss->sID]->shrink_to_fit();
      densityPerSlipSystem[ ss->sID]->clear();
      densityPerSlipSystem[ ss->sID]->shrink_to_fit();
      glissileDensityPerSlipSystem[ ss->sID]->clear();
      glissileDensityPerSlipSystem[ ss->sID]->shrink_to_fit();
      sessileDensityPerSlipSystem[ ss->sID]->clear();
      sessileDensityPerSlipSystem[ ss->sID]->shrink_to_fit();
   }

   // measurements indexed by tensor indices
   for ( const auto& key : tensorComponentKeys)//11,22,33,12,13,23
   {
      stressTensorComponents[ key]->clear();
      stressTensorComponents[ key]->shrink_to_fit();
      //strainTensorComponents[ key]->clear();
      //strainTensorComponents[ key]->shrink_to_fit();
   }

   return;
}

//void ddpy::DDInterface::generateLoopDensity()
//void ddpy::DDInterface::generateDipoleDensity()
//void ddpy::DDInterface::generateDipole()
//{
//   // TODO:  emulate MicrostructureGenerator.cpp:57-108
//   // Can a periodicDipoleGenerator be instantiated without reading
//   //   microstructure parameters from a file? No.
//   // maybe create alternatives to the following which don't read from files:
//   //     MicrostructureGeneratorBase, PeriodicDipoleGenerator, ...
//   //
//   // PeriodicDipoleGenerator inherits from MicrostructureGeneratorBase
//   //
//   // MicrostructureGeneratorBase instantiates a TextFileParser and reads
//   //  type
//   //  style
//   //  tag
//   //
//   // PeriodicDipoleGenerator::generateDensity uses parser to read a file
//   //
//   // I don't want to instantiate a MicrostructureGenerator
//   //  because it will attempt to read the files listed in traitsIO.microstructureFile
//   //TODO: empty traitsIO.microstructureFile
//
//   // things to replace from TextFileParser parser found in periodicDipoleGenerator::generateIndividual:
//   const std::vector<int> periodicDipoleExitFaceIDs(this->parser.readArray<int>("periodicDipoleExitFaceIDs",true));
//   const Eigen::Matrix<double,Eigen::Dynamic,dim> periodicDipolePoints(this->parser.readMatrix<double>("periodicDipolePoints",periodicDipoleSlipSystemIDs.size(),dim,true));
//   const std::vector<double> periodicDipoleHeights(this->parser.readArray<double>("periodicDipoleHeights",true));
//   const std::vector<int> periodicDipoleNodes(this->parser.readArray<int>("periodicDipoleNodes",true));
//   const std::vector<double> periodicDipoleGlideSteps(this->parser.readArray<double>("periodicDipoleGlideSteps",true));
//
//   return;
//}

//void ddpy::DDInterface::regenerateMicrostructure()
//{
//   if ( ddBase == nullptr) readddBase();
//   //
//   // instantiate a MicrostructureGenerator
//   model::MicrostructureGeneratorInMem mg( *ddBase); // requires MicrostructureSpecifications
//   if ( debugFlag)
//      std::cout << "finished call to MicrostructureGeneratorInMem" << std::endl;
//
//   setCurrentStep(0); // DefectiveCrystal will use runID to read some things
//
//   // instantiate a DefectiveCrystalType DC( ddBase)
//   DC = std::unique_ptr<DefectiveCrystalType>(
//         new DefectiveCrystalType( *ddBase)
//         );
//   if ( debugFlag)
//      std::cout << "finished call to DefectiveCrystalType" << std::endl;
//
//   return;
//}

void ddpy::DDInterface::resetStaticIDs()
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

void ddpy::DDInterface::specifyLoops(
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

void ddpy::DDInterface::specifyLoopDensitiesPerSlipSystem(
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

void ddpy::DDInterface::specifyLoopDensity(
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
void ddpy::DDInterface::specifyPrismaticLoops(
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

void ddpy::DDInterface::specifyPrismaticLoopDensitiesPerSlipSystem(
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

void ddpy::DDInterface::specifyPrismaticLoopDensity(
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
void ddpy::DDInterface::specifyDipoleDensity(
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


void ddpy::DDInterface::specifyDipoles(
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

   std::cout << "pointsIn.shape(): " << pointsIn.shape() << std::endl; // debug
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

void ddpy::DDInterface::setOutputPath( const std::string& outputPath)
{ // TODO: change return value to allow error checking
   if ( ddBase == nullptr) readddBase();
   // if the path doesn't exist, try to make it
   //if ( ! std::filesystem::is_directory( outputPath))
   //{
   //   // TODO: create the folders: outputPath
   //}
   //if (! std::filesystem::is_directory( outputPath + "/evl"))
   //{
   //   try
   //   {
   //      // TODO: create the folders: {outputPath}/evl
   //   }
   //   catch
   //   {
   //      return false;
   //   }
   //}
   //if (! std::filesystem::is_directory( outputPath + "/F"))
   //{
   //   try
   //   {
   //      // TODO: create the folders: {outputPath}/F
   //   }
   //   catch
   //   {
   //      return false;
   //   }
   //}
   //if (! std::filesystem::is_directory( outputPath + "/evl"))
   //{
   //   std::cout << "not found: " << outputPath + "/evl" << std::endl;
   //}
   //else
   //{
   //   std::cout << "found: " << outputPath + "/evl" << std::endl;
   //}
   //if ( DC == nullptr)
   //{
   //   std::cout << "error: defective crystal not yet read" << std::endl;
   //   return;
   //}
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
   return;
}

void ddpy::DDInterface::regeneratePolycrystalFile(
            const pybind11::array_t<double,
               pybind11::array::c_style | pybind11::array::forcecast>
               c2gIn,
            const std::string& lattice,
            const std::string& material,
            const std::string& meshFilePath
            )
{
   MatrixDim c2g;
   auto c2gNp = c2gIn.unchecked<2>(); // for reading c2g input np.array
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
      if ( lattice == std::string(lat)) latticeIsAcceptable = true;
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
   }
   std::cout << "lattice: " << lattice << ", material: " << material << std::endl; // debug

   c2g << c2gNp(0,0), c2gNp(0,1), c2gNp(0,2),
         c2gNp(1,0), c2gNp(1,1), c2gNp(1,2),
         c2gNp(2,0), c2gNp(2,1), c2gNp(2,2);
   //std::cout << "c2g:\n" << c2g << std::endl; // debug
   std::string outputFilePath = dddFolderPath + "/inputFiles/polycrystal.txt";

   // TODO: detect and move any existing polycrystal.txt file
   std::ofstream outputFile( outputFilePath,
              std::ofstream::out | std::ofstream::trunc);

   //std::cout << "lammps boundaries: " <<  std::endl;// debug
   //for ( const auto& bd : boxBounds) std::cout << bd << ", "; // debug
   //std::cout << std::endl; // debug

   double deltaX, deltaY, deltaZ;
   deltaX = abs( boxBounds[1] - boxBounds[0]);
   deltaY = abs( boxBounds[3] - boxBounds[2]);
   deltaZ = abs( boxBounds[5] - boxBounds[4]);

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

   AA = c2g * AA;

   MatrixDim f12, f31, f23;
   f12 << 1.0, boxSkew, 0.0,
       0.0, 1.0, 0.0,
       0.0, 0.0, 1.0;
   f31 << 1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       boxSkew, 0.0, 1.0;
   f23 << 1.0, 0.0, 0.0,
       0.0, 1.0, boxSkew,
       0.0, 0.0, 1.0;
   MatrixDim deformingMatrix;
   deformingMatrix = f12 * (f23 * f31);

   MatrixDim scalingMatrix;
   scalingMatrix << deltaX, 0, 0,
      0, deltaY, 0,
      0, 0, deltaZ;
   deformingMatrix = ( AA.inverse()) * (deformingMatrix * scalingMatrix);
   deformingMatrix = deformingMatrix.array().round();
   AA = AA * deformingMatrix;

   double x0x, x0y, x0z; // boxBounds required to identify x0
   x0x = boxBounds[0]/deltaX; // boxBounds[0] is xlo, deltaX = xhi -xlo
   x0y = boxBounds[2]/deltaY; // boxBounds[2] is ylo
   x0z = boxBounds[4]/deltaZ; // boxBounds[4] is zlo
   VectorDim x0;
   x0 << x0x, x0y, x0z; // scaling of the mesh is y=A(x-x0)

   outputFile << "materialFile=" << material + ".txt;" << std::endl;
   outputFile << "enablePartials=0;"
      << std::endl;
   outputFile << "absoluteTemperature = 300; # [K] simulation temperature"
      << std::endl;
   outputFile << "meshFile=" << meshFilePath << ";"
      << std::endl;
   outputFile << "C2G1="
        << std::setw(22) << std::setprecision(15) << c2g(0,0)
        << std::setw(22) << std::setprecision(15) << c2g(0,1)
        << std::setw(22) << std::setprecision(15) << c2g(0,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << c2g(1,0)
        << std::setw(22) << std::setprecision(15) << c2g(1,1)
        << std::setw(22) << std::setprecision(15) << c2g(1,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << c2g(2,0)
        << std::setw(22) << std::setprecision(15) << c2g(2,1)
        << std::setw(22) << std::setprecision(15) << c2g(2,2)
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
   outputFile << std::endl;

   if (! lattice.compare("bcc")) // .compare() returns 0 if they're equal
      outputFile << "dislocationMobilityType='BCC';" << std::endl;
   else if (! lattice.compare("fcc")) // .compare() returns 0 if they're equal
      outputFile << "dislocationMobilityType='FCC';" << std::endl;
   else
   {
      std::cout << "error: lattice type not recognized while trying to "
         << "write polycrystal.txt" << std::endl;
      return;
   }
   //outputFile << "solidSolutionNoiseFile_xz=" << ";" << std::endl;
   //outputFile << "solidSolutionNoiseFile_yz=" << ";" << std::endl;

   return;
}

void ddpy::DDInterface::clearMicrostructureSpecifications()
{
   microstructureSpecifications.clear();
   return;
}

pybind11::array_t<double, pybind11::array::c_style>
   ddpy::nonUniformConvolutionWithAGaussian(
      const pybind11::array_t<double, pybind11::array::c_style>& timesIn,
      const pybind11::array_t<double, pybind11::array::c_style>& yyIn,
      const ssize_t& smoothedStartIdx,
      const ssize_t& smoothedEndIdx,
      const double& sigma,
      const ssize_t& sigmaCount
      )
{
//#ifdef _OPENMP
//   const size_t nThreads = omp_get_max_threads();
//#else
//   const size_t nThreads = 1;
//#endif
   auto times = timesIn.unchecked<1>();
   auto yy = yyIn.unchecked<1>();
   // check input validity
   if ( times.shape(0) != yy.shape(0))
   {
      std::cout << "error: nonUniformConvolutionWithAGaussian() was given domain and range arrays of differing sizes, "
         << "times.shape(0) " << times.shape(0)
         << ", yy.shape(0) " << yy.shape(0)
         << std::endl; 
      pybind11::array_t<double, pybind11::array::c_style> tmp;
      return tmp;
   }
   if ( smoothedStartIdx >= smoothedEndIdx)
   {
      std::cout << "error: nonUniformConvolutionWithAGaussian() was given"
        << " smoothedStartIdx " << smoothedStartIdx
        << " <= smoothedEndIdx " << smoothedEndIdx;
      pybind11::array_t<double, pybind11::array::c_style> tmp;
      return tmp;
   }
   if ( ( smoothedStartIdx < 0) || ( smoothedStartIdx > yy.shape(0) -1))
   {
      std::cout << "error: nonUniformConvolutionWithAGaussian() was given"
         << " an out of bounds smoothedStartIdx of " << smoothedStartIdx
         << ", while yy.shape(0) is " << yy.shape(0)
         << std::endl; 
      pybind11::array_t<double, pybind11::array::c_style> tmp;
      return tmp;
   }
   if ( ( smoothedEndIdx < 0) || ( smoothedEndIdx > yy.shape(0) -1))
   {
      std::cout << "error: nonUniformConvolutionWithAGaussian() was given"
         << " an out of bounds smoothedEndIdx of " << smoothedStartIdx
         << ", while yy.shape(0) is " << yy.shape(0)
         << std::endl; 
      pybind11::array_t<double, pybind11::array::c_style> tmp;
      return tmp;
   }

   // allocate vector to return as output and initialize with 0.0 values
   pybind11::array_t<double, pybind11::array::c_style> smoothedyy ( times.shape(0));
   py::buffer_info smoothedyybuf = smoothedyy.request();
   double *smoothedyyptr = static_cast< double *>( smoothedyybuf.ptr);
   for ( ssize_t idx=0; idx < smoothedyybuf.shape[0]; ++idx)
   {
      smoothedyyptr[idx] = 0.0;
   }

   //std::shared_ptr< std::vector< double> > smoothedyy;
   //smoothedyy = std::make_shared< std::vector<double> >( times.shape(0), 0.0);
   //py::buffer_info smoothedyybuf = smoothedyy->request();
   //double* smoothedyyptr = static_cast< double *>( smoothedyybuf.ptr);
   ////smoothedyy->resize( times.size(), 0.0);
   ////smoothedyy( yy.shape(0), 0.0);
   //std::vector< double> smoothedyy( times.shape(0), 0.0);

   double sigmaCountSigmas = sigmaCount * sigma;
   //std::cout << "nThreads: " << nThreads << std::endl; // debug
   // iterate yIdx from yy[ smoothedStartIdx] to yy[ smoothedEndIdx]
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for ( ssize_t yIdx = smoothedStartIdx; yIdx < smoothedEndIdx; ++yIdx)
   {
      // iterate tIdx over all but the final element of times
      for ( ssize_t tIdx=0; tIdx < times.size(); ++tIdx)
      {
         if ( abs( times[ tIdx] - times[ yIdx]) < sigmaCountSigmas)
         {
            smoothedyyptr[ yIdx]
               += yy[ tIdx]
                  * gaussian( times[ tIdx], sigma, times[ yIdx])
                  * ( times[ tIdx+1] - times[ tIdx]);
         }
      }
   }
   return smoothedyy;
}

pybind11::array_t<double, pybind11::array::c_style>
   ddpy::nonUniformWeightByAGaussian(
      const pybind11::array_t<double, pybind11::array::c_style>& timesIn,
      const pybind11::array_t<double, pybind11::array::c_style>& yyIn,
      const pybind11::array_t<ssize_t, pybind11::array::c_style>& IdxsIn,
      const double& sigma,
      const ssize_t& sigmaCount
      )
{
   auto times = timesIn.unchecked<1>();
   auto yy = yyIn.unchecked<1>();
   auto Idxs = IdxsIn.unchecked<1>();
   double sigmaCountSigmas = sigmaCount * sigma; // in same units as times
   // check input validity
   if ( times.shape(0) != yy.shape(0))
   {
      std::cout << "error: nonUniformWeightByAGaussian() was given domain and range arrays of differing sizes, "
         << "times.shape(0) " << times.shape(0)
         << ", yy.shape(0) " << yy.shape(0)
         << "; returning an empty array"
         << std::endl; 
      pybind11::array_t<double, pybind11::array::c_style> tmp;
      return tmp;
   }
   ssize_t tmpIdx;
   for ( ssize_t IdxsIdx=0; IdxsIdx < Idxs.shape(0); ++IdxsIdx)
   {
      tmpIdx = Idxs[ IdxsIdx];
      if (( tmpIdx < 0) || (tmpIdx > yy.shape(0) -1))
      {
         std::cout << "error: nonUniformWeightByAGaussian() was given"
            << " an out of bounds of index" << tmpIdx
            << ", while yy.shape(0) is " << yy.shape(0)
            << "; returning an empty array"
            << std::endl;
         pybind11::array_t<double, pybind11::array::c_style> tmp;
         return tmp;
      }
      if (
            abs( times[tmpIdx] - times[0]) < sigmaCountSigmas
            ||
            abs( times[tmpIdx] - times[times.shape(0)-1]) < sigmaCountSigmas
         )
      {
         std::cout << "error: nonUniformWeightByAGaussian() was given an index "
            << tmpIdx
            << " to a time that is within sigmaCount * sigma"
            << sigmaCount << "*" << sigma
            << " from the first or last time value "
            << " times[0] " << times[0] << ", times[times.shape(0)-1] " 
            << times[times.shape(0)-1]
            << "; returning an empty array"
            << std::endl;
         pybind11::array_t<double, pybind11::array::c_style> tmp;
         return tmp;
      }
   }
   pybind11::array_t<double, pybind11::array::c_style> weightedValues( Idxs.shape(0));
   py::buffer_info weightedValuesbuf = weightedValues.request();
   double *weightedValuesptr = static_cast< double *>( weightedValuesbuf.ptr);
   for ( ssize_t idx=0; idx < weightedValuesbuf.shape[0]; ++idx)
   {
      weightedValuesptr[idx] = 0.0;
   }
   // TODO: adapt lines 2526-2540 to only weight around indices in Idxs
   ssize_t wvIdx; wvIdx=0; // weightedValues has the same length as Idxs
   ssize_t centerIdx;
   for ( ssize_t IdxsIdx=0; IdxsIdx < Idxs.shape(0); ++IdxsIdx)
   {
      centerIdx = Idxs[ IdxsIdx];
      // iterate neighIdx over all but the final element of times
      for ( ssize_t neighIdx=0; neighIdx < times.size()-1; ++neighIdx)
      {
         if ( abs( times[ neighIdx] - times[ centerIdx]) < sigmaCountSigmas)
         {
            weightedValuesptr[ wvIdx]
               += yy[ neighIdx]
                  * gaussian( times[ neighIdx], sigma, times[ centerIdx])
                  * ( times[ neighIdx+1] - times[ neighIdx]);
         }
      }
      ++wvIdx;
   }
   return weightedValues;
}


PYBIND11_MODULE( ddpy, m) {
   namespace py = pybind11;
   m.doc() = "TODO: revise m.doc() in src/modelib2py.cpp";
   py::class_<ddpy::DDInterface>( m, "DDInterface")
      .def( py::init(
               // using a lambda function that returns an instantiation
               [](const std::string& modelibFolderPath)
               {
                  return std::unique_ptr< ddpy::DDInterface>(
                        new ddpy::DDInterface( modelibFolderPath)
                        );
               }), py::arg("modelibFolderPath").none(false)
            )
      .def("getCurrentStep",
            &ddpy::DDInterface::getCurrentStep
          )
      .def("setCurrentStep",
            &ddpy::DDInterface::setCurrentStep,
            py::arg("currentStep").none(false)
          )
      .def("setOutputFrequency",
            &ddpy::DDInterface::setOutputFrequency,
            py::arg("outputFrequency").none(false)
          )
      .def("runGlideSteps",
            &ddpy::DDInterface::runGlideSteps,
            py::arg("Nsteps").none(false)
          )
      .def("getResolvedShearStresses",
            &ddpy::DDInterface::getResolvedShearStresses
          )
      .def("getResolvedShearStrains",
            &ddpy::DDInterface::getResolvedShearStrains
          )
      .def("getPlasticStrains",
            &ddpy::DDInterface::getPlasticStrains
          )
      .def("getBurgersMagnitude",
            &ddpy::DDInterface::getBurgersMagnitude
          )
      .def("getDensityPerSlipSystem",
            &ddpy::DDInterface::getDensityPerSlipSystem
          )
      .def("readBurgersMagnitude",
            &ddpy::DDInterface::readBurgersMagnitude,
            py::arg("materialFilePath").none(false)
          )
      .def("regeneratePolycrystalFile",
            &ddpy::DDInterface::regeneratePolycrystalFile,
            py::kw_only(),
            py::arg( "c2g").none(false),
            py::arg( "lattice").none(false),
            py::arg( "material").none(false),
            py::arg( "meshFilePath").none(false)
          )
      .def("specifyLoops",
            &ddpy::DDInterface::specifyLoops,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "slipSystemIDs").none(false),
            py::arg( "loopRadii").none(false),
            py::arg( "loopSegmentCounts").none(false),
            py::arg( "loopCenters").none(false)
          )
      .def("specifyDipoles",
            &ddpy::DDInterface::specifyDipoles,
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
            &ddpy::DDInterface::specifyLoopDensity,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "loopDensity").none(false),
            py::arg( "loopSegmentCount").none(false),
            py::arg( "loopRadiusMean").none(false),
            py::arg( "loopRadiusStd").none(false)
          )
      .def("specifyLoopDensitiesPerSlipSystem",
            &ddpy::DDInterface::specifyLoopDensitiesPerSlipSystem,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "loopDensitiesPerSlipSystem").none(false),
            py::arg( "loopSegmentCount").none(false), // keep
            py::arg( "loopRadiusMean").none(false), // keep
            py::arg( "loopRadiusStd").none(false) // keep
          )
      .def("specifyPrismaticLoops", // TODO: implememnt
            &ddpy::DDInterface::specifyPrismaticLoops,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "slipSystemIDs").none(false),
            py::arg( "prismaticLoopRadii").none(false),
            py::arg( "prismaticLoopCenters").none(false),
            py::arg( "prismaticLoopSteps").none(false)
          )
      .def("specifyPrismaticLoopDensity", // TODO: implememnt
            &ddpy::DDInterface::specifyPrismaticLoopDensity,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "prismaticLoopDensity").none(false),
            py::arg( "prismaticLoopRadiusMean").none(false),
            py::arg( "prismaticLoopRadiusStd").none(false),
            py::arg( "prismaticLoopStepMean").none(false),
            py::arg( "prismaticLoopStepStd").none(false)
          )
      .def("specifyPrismaticLoopDensitiesPerSlipSystem", // TODO: implememnt
            &ddpy::DDInterface::specifyPrismaticLoopDensitiesPerSlipSystem,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "prismaticLoopDensitiesPerSlipSystem").none(false),
            py::arg( "prismaticLoopRadiusMean").none(false), // keep
            py::arg( "prismaticLoopRadiusStd").none(false), // keep
            py::arg( "prismaticLoopStepMean").none(false), // keep
            py::arg( "prismaticLoopStepStd").none(false) // keep
          )
      .def("specifyDipoleDensity",
            &ddpy::DDInterface::specifyDipoleDensity,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "dipoleDensity").none(false)
          )
      .def("generateMicrostructure",
         &ddpy::DDInterface::generateMicrostructure
         )
      .def("writeConfigToTxt",
            &ddpy::DDInterface::writeConfigToTxt
          )
      .def("readDefectiveCrystal",
            &ddpy::DDInterface::readDefectiveCrystal
          )
      .def("setOutputPath",
            &ddpy::DDInterface::setOutputPath,
            py::arg( "outputPath").none(false)
          )
      .def("setBoxBounds", // prerequisite: readBurgersMagnitude()
            &ddpy::DDInterface::setBoxBounds,
            py::arg( "xlo").none(false), // [m]
            py::arg( "xhi").none(false),
            py::arg( "ylo").none(false),
            py::arg( "yhi").none(false),
            py::arg( "zlo").none(false),
            py::arg( "zhi").none(false)
      )
      .def("clearMicrostructureSpecifications",
            &ddpy::DDInterface::clearMicrostructureSpecifications
      )
      .def("setExternalLoad", // prerequisite: readBurgersMagnitude()
            &ddpy::DDInterface::setExternalLoad,
            py::kw_only(),
            py::arg( "stress") = py::none(), //.none(true), // 3x3 matrix
            py::arg( "stressRate") = py::none(), // 3x3 matrix
            py::arg( "strain") = py::none(), // 3x3 matrix
            py::arg( "strainRate") = py::none(), // 3x3 matrix
            py::arg( "stiffnessRatio") = py::none() // Voigt format 11 22 33 12 23 13
          )
      .def("enableMechanicalMeasurements",
            &ddpy::DDInterface::enableMechanicalMeasurements,
            py::arg( "stepsBetweenMeasurements") = 1
          )
      .def("disableMechanicalMeasurements",
            &ddpy::DDInterface::disableMechanicalMeasurements
          )
      .def("getMechanicalMeasurements", // TODO
            &ddpy::DDInterface::getMechanicalMeasurements // TODO
          )
      .def("clearMechanicalMeasurements", // TODO
            &ddpy::DDInterface::clearMechanicalMeasurements // TODO
          )

      //.def("getSlipSystemNormals"
      //      &ddpy::DDInterface::getSlipSystemNormals
      //    )
      //.def("getSlipSystemBurgersVectors"
      //      &ddpy::DDInterface::getSlipSystemBurgersVectors
      //    )
      ;
   m.def("nonUniformConvolutionWithAGaussian",
         &ddpy::nonUniformConvolutionWithAGaussian,
         //py::kw_only(),
         py::arg( "x").none(false),
         py::arg( "y").none(false),
         py::arg( "startIdx").none(false),
         py::arg( "endIdx").none(false),
         py::arg( "sigma").none(false),
         py::arg( "sigmaCount").none(false)
        );
         //"smoothes the input via Gaussian convolution");
   m.def("nonUniformWeightByAGaussian",
         &ddpy::nonUniformWeightByAGaussian,
         //py::kw_only(),
         py::arg( "x").none(false),
         py::arg( "y").none(false),
         py::arg( "indices").none(false),
         py::arg( "sigma").none(false),
         py::arg( "sigmaCount").none(false)
        );
   //py::class_<ddpy::DDInterface>( m, "DDInterface")
   //   .def(
   //         py::init(
   //            // using a lambda function that returns an instantiation
   //            [](const std::string& folderName)
   //            {
   //               model::DislocationDynamicsBase ddBase( folderName);
   //               return std::unique_ptr< model::DefectiveCrystal>(
   //                     new model::DefectiveCrystal( ddBase)
   //                     );
   //            }
   //            ), py::arg("dddFolderPath")
   //         )
   //   .def("getCurrentStep",
   //         &ddpy::DDInterface::getCurrentStep
   //         )
   //   .def("setCurrentStep",
   //         &ddpy::DDInterface::setCurrentStep,
   //         py::arg("currentStep")
   //       )
   //   .def("runGlideSteps",
   //         &ddpy::DDInterface::runGlideSteps
   //       );
   //   // TODO: expose any other class members to pybind here
   //   //.def("printDisplacements", &ddpy::DDInterface::printDisplacements)
   //   ;
} //  PYBIND11_MODULE

#endif
