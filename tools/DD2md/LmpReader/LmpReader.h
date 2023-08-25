
#ifndef dd2md_LmpReader_H_
#define dd2md_LmpReader_H_

#include <iostream>
#include <fstream>
#include <ostream> // necessary?
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h> // EXIT_SUCCESS, EXIT_FAILURE

#include <Polycrystal.h>
#include <DefectiveCrystal.h>
#include <UniformExternalLoadController.h>
#include <IDreader.h>
//#include <VTKsegments.h>

#include <Eigen/Core>
#include <Eigen/LU>

namespace model
{
class LmpReader
{
   public:
   typedef Eigen::Matrix<double,3,1> VectorDim;
   typedef Eigen::Matrix<double,3,3> MatrixDim;
   //typedef std::vector<model::FEMnodeEvaluation<model::LagrangeElement<3,2>,3,1>> atomPositionType;
   typedef std::vector<VectorDim> atomPositionType;
   typedef std::vector<size_t> pointIDsType;

   private:
   int skipLinesOfWhitespace(
         std::ifstream& inFile,
         std::string& line
         );
   void outputCommonErrorMessage(
         std::ostream& outStream,
         const std::string& lammpsFilePath,
         const size_t& lineNumber,
         const std::string& line
         );
   size_t lineNumber;
   std::string lammpsFilePath;
   //const model::DefectiveCrystal<3,0>* const DC;
   double scaleFactor;
   MatrixDim deformationMatrix;
   bool debugFlag;


   public:
   LmpReader(
         const std::string& lmpFilePath,
         const double& scaleFactorIn,
         const MatrixDim& deformationMatrixIn,
         const bool& dbgFlag
         ) :
        lammpsFilePath( lmpFilePath)
      , scaleFactor( scaleFactorIn)
      , deformationMatrix( deformationMatrixIn)
      , debugFlag( dbgFlag)
   {};

   int readLmpStream(
         std::vector<double>& bounds,
         std::map<size_t, double>& masses,
         std::vector<size_t>& atomIDs,
         std::vector<size_t>& atomTypes,
         atomPositionType& atomPositions,
         pointIDsType& pointIDs
         );
   int readLmpStreamBounds( std::vector<double>& bounds);
   void trimWhitespace( std::string& str);
};
} // namespace model
#endif
