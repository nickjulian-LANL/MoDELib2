
#ifndef dd2md_LmpReader_CPP_
#define dd2md_LmpReader_CPP_

#include <LmpReader.h>

int model::LmpReader::skipLinesOfWhitespace(
      std::ifstream& inFile,
      std::string& line)
{
   size_t first = line.find_first_not_of(" \t");
   //std::cout << "attempting to skip empty lines" << std::endl; // debug
   while ( first == std::string::npos) // npos : max size of a string
   {
      ++lineNumber;
      if (!( std::getline( inFile, line) && inFile.good()))
      {
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
         }
         inFile.close();
         return EXIT_FAILURE;
      }
      first = line.find_first_not_of(" \t");
   }
   //std::cout << "skipped to line " << lineNumber << std::endl; // debug
   return EXIT_SUCCESS;
}

int model::LmpReader::readLmpStream(
         std::vector<double>& bounds,
         std::map<size_t, double>& masses,
         std::vector<size_t>& atomIDs,
         std::vector<size_t>& atomTypes,
         atomPositionType& atomPositions
      )
{
   // Following specification found at:
   //   https://docs.lammps.org/2001/data_format.html
   // required sections: header, Masses, Atoms

   // empty the vectors to write be written to
   if ( bounds.size() != 6) bounds.clear();
   if ( masses.size() != 0) masses.clear();
   if ( atomIDs.size() != 0) atomIDs.clear();
   if ( atomTypes.size() != 0) atomTypes.clear();
   if ( atomPositions.size() != 0) atomPositions.clear();

   // open the file
   std::ifstream inFile(lammpsFilePath.c_str(),std::ifstream::in);
   if (!( inFile.good()))
   {
      std::cout << "Error: failed to open the file "
        << lammpsFilePath << std::endl;
      inFile.close();
      return EXIT_FAILURE;
   }

   //std::cout<<"Reading file "<< lammpsFilePath <<"..."<<std::endl; // debug
   std::string line;
   std::stringstream ss;
   std::string sectionName;
   lineNumber = 0;
   size_t atomPopulationCount; atomPopulationCount = 0;
   size_t atomTypeCount; atomTypeCount = 0;
   size_t tmpAtomType; tmpAtomType = 0;
   size_t tmpAtomID; tmpAtomID = 0;
   double tmpMass; tmpMass = 0;
   double tmpDbl1, tmpDbl2, tmpDbl3; tmpDbl1 = 0; tmpDbl2 = 0; tmpDbl3 = 0;
   std::string term1, term2, term3, term4, term5, term6;

   // throw away the first two lines
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      std::cout << "Error: failed to read the first line of file "
         << lammpsFilePath << std::endl;
      inFile.close();
      return EXIT_FAILURE;
   }
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      std::cout << "Error: failed to read the first line of file "
         << lammpsFilePath << std::endl;
      inFile.close();
      return EXIT_FAILURE;
   }
   
   // read the header section
   sectionName = "header";
   // header: 
   //    500 atoms
   //    0 bonds     #skip, only required if there are > 0 of them
   //    0 angles    #skip, only required if there are > 0 of them
   //    0 dihedrals #skip, only required if there are > 0 of them
   //    0 impropers #skip, only required if there are > 0 of them
   //    # force field coefficients?
   //
   //    1 atom types
   //    0 bond types      #skip, only required if there are > 0 of them
   //    0 angle types     #skip, only required if there are > 0 of them
   //    0 dihedral types  #skip, only required if there are > 0 of them
   //    0 improper types  #skip, only required if there are > 0 of them
   //
   //    -1 1 xlo xhi
   //    -1 1 ylo yhi
   //    -1 1 zlo zhi

   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   lineStream >> atomPopulationCount;
   lineStream >> term2;
   }
   //std::cout << "line: " << line << std::endl; // debug
   //std::cout << "atomPopulationCount: " << atomPopulationCount << std::endl; // debug
   // match to atomPopulationCount ' atoms'
   if ( (! atomPopulationCount)
         || (term2.compare("atoms")) // 0 if they're equal
         )
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                     lineNumber, line);
         std::cout << "current line should be <size_t> atoms" << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   //if ( debugFlag )
   //{
   //   std::cout << atomPopulationCount << " atoms declared"
   //   << std::endl;
   //}

   // skip '<size_t> {bonds,angle,dihedrals,impropers}' lines++lineNumber;
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   term1.clear();
   term2.clear();
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   lineStream >> term1;
   lineStream >> term2;
   }
   while (!( (term2.compare("bonds"))
         || (term2.compare("angles"))
         || (term2.compare("dihedrals"))
         || (term2.compare("impropers"))
         ))
   {
      //if ( debugFlag) std::cout << "skipping line " << line << std::endl;
      ++lineNumber;
      if (!( std::getline(inFile, line) && inFile.good()))
      {
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
         }
         inFile.close();
         return EXIT_FAILURE;
      }
      std::istringstream lineStream( line);
      term1.clear();
      term2.clear();
      lineStream >> term1;
      lineStream >> term2;
   }
   if ( skipLinesOfWhitespace( inFile, line) != EXIT_SUCCESS)
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                               lineNumber, line);
         std::cout << "LmpReader failed during skipLinesOfWhitespace"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   term1.clear();
   term2.clear();
   term3.clear();
   lineStream >> term1;
   lineStream >> term2;
   lineStream >> term3;
   }
   // term1,2,3 should be "<size_t> atom types"
   if ( term2.compare("atom") || term3.compare("types"))
   {
      // failure
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
         std::cout << "LmpReader, line should be '<size_t> atom types'"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( term1);
   lineStream >> atomTypeCount;
   }

   // read system x bounds
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   term1.clear();
   term2.clear();
   lineStream >> term1;
   lineStream >> term2;
   }
   while (!( (term2.compare("bonds"))
         || (term2.compare("angles"))
         || (term2.compare("dihedrals"))
         || (term2.compare("impropers"))
         ))
   { // skip '<size_t> {bonds,angle,dihedrals,impropers} types' lines
      ++lineNumber;
      if (!( std::getline(inFile, line) && inFile.good()))
      {
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
         }
         inFile.close();
         return EXIT_FAILURE;
      }
      //std::istringstream lineStream2( line);
      //term1.clear();
      //term2.clear();
      //lineStream2 >> term1;
      //lineStream2 >> term2;
   }
   if ( skipLinesOfWhitespace( inFile, line) != EXIT_SUCCESS)
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                               lineNumber, line);
         std::cout << "LmpReader failed during skipLinesOfWhitespace"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   //term3.clear();
   //term4.clear();
   lineStream >> tmpDbl1;
   lineStream >> tmpDbl2;
   lineStream >> term3;
   lineStream >> term4;
   }
   // term1,2,3,4 should be "<double> <double> xlo xhi"
   if ( term3.compare("xlo") || term4.compare("xhi"))
   {
      // failure
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
         std::cout << "LmpReader, line should be "
            << "'<double> <double> xlo xhi" << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   bounds.emplace_back( tmpDbl1 * scaleFactor);
   bounds.emplace_back( tmpDbl2 * scaleFactor);

   // read system y bounds
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   lineStream >> tmpDbl1;
   lineStream >> tmpDbl2;
   lineStream >> term3;
   lineStream >> term4;
   }
   // term1,2,3,4 should be "<double> <double> ylo yhi"
   if ( term3.compare("ylo") || term4.compare("yhi"))
   {
      // failure
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
         std::cout << "LmpReader, line should be "
            << "'<double> <double> ylo yhi" << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   bounds.emplace_back( tmpDbl1 * scaleFactor);
   bounds.emplace_back( tmpDbl2 * scaleFactor);

   // read system z bounds
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   lineStream >> tmpDbl1;
   lineStream >> tmpDbl2;
   lineStream >> term3;
   lineStream >> term4;
   }
   // term1,2,3,4 should be "<double> <double> zlo zhi"
   if ( term3.compare("zlo") || term4.compare("zhi"))
   {
      // failure
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
         std::cout << "LmpReader, line should be "
            << "'<double> <double> zlo zhi" << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   bounds.emplace_back( tmpDbl1 * scaleFactor);
   bounds.emplace_back( tmpDbl2 * scaleFactor);
   //if ( debugFlag)
   //{
   //   std::cout << "bounds: " << std::endl
   //      << "  " << bounds[0] << " " << bounds[1] << std::endl
   //      << "  " << bounds[2] << " " << bounds[3] << std::endl
   //      << "  " << bounds[4] << " " << bounds[5] << std::endl;
   //}
   // see if box is sheared
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   lineStream >> tmpDbl1;
   lineStream >> tmpDbl2;
   lineStream >> tmpDbl3;
   lineStream >> term4;
   lineStream >> term5;
   lineStream >> term6;
   }
   // term1,2,3,4,5 might be "<double> <double> <double> xy xz yz"
   if ( !term4.compare("xy"))
   {
      if ( !term5.compare("xz"))
      {
         if ( !term6.compare("yz"))
         {
            if ((tmpDbl1 != 0.0) || (tmpDbl2 != 0.0) || (tmpDbl3 != 0.0))
            {
               // failure
               if ( debugFlag)
               {
                  outputCommonErrorMessage( std::cout, lammpsFilePath,
                                             lineNumber, line);
                  std::cout << "LmpReader, box tilt factors must be 0"
                     << " line: " << line
                     << std::endl;
               }
               inFile.close();
               return EXIT_FAILURE;
            }
         }
         else
         {
            // failure
            if ( debugFlag)
            {
               outputCommonErrorMessage( std::cout, lammpsFilePath,
                                             lineNumber, line);
               std::cout << "LmpReader, unkown line format" << std::endl;
            }
            inFile.close();
            return EXIT_FAILURE;
         }
      }
      else
      {
         // failure
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                          lineNumber, line);
            std::cout << "LmpReader, unkown line format" << std::endl;
         }
         inFile.close();
         return EXIT_FAILURE;
      }
   } // if ( term4.compare("xy")) // if not, then continue


   /********************************************************************/
   /* read sections: Masses, Atoms, and discard other sections *********/
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   if ( skipLinesOfWhitespace( inFile, line) != EXIT_SUCCESS)
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                               lineNumber, line);
         std::cout << "LmpReader failed during skipLinesOfWhitespace"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   term1.clear();
   lineStream >> term1; // trims whitespace
   }
   if ( term1.compare("Masses"))
   {
      // failure
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
         std::cout << "LmpReader, line should be "
            << "Masses" << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   sectionName = "Masses";
   // skip blank lines
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   if ( skipLinesOfWhitespace( inFile, line) != EXIT_SUCCESS)
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                               lineNumber, line);
         std::cout << "LmpReader failed during skipLinesOfWhitespace"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   // read atomTypeCount number of masses
   for ( size_t ii=0; ii < atomTypeCount; ++ii)
   {
      { // scope so that lineStream can be redeclared later
      std::istringstream lineStream( line);
      lineStream >> tmpAtomType;
      //if ( (tmpAtomType <= 0) || (tmpAtomType > atomTypeCount))
      if ( (tmpAtomType <= 0) || (tmpAtomType > 118))
      {
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
            std::cout << "((" << tmpAtomType << " <= 0) || ("
               << tmpAtomType << " > 118 ))"
               << std::endl;
            std::cout <<  "Error: atom type integer must be among "
               << "1 through 118" << std::endl;
         }
         inFile.close();
         return EXIT_FAILURE;
      }
      lineStream >> tmpMass;
      }
      masses[ tmpAtomType] = tmpMass;
      if ( masses[ tmpAtomType] <=0)
      {
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
            std::cout <<  "Error: atom mass "
               << masses[ tmpAtomType]
               << " of type " << tmpAtomType << " must be > 0"
               << std::endl;
         }
         inFile.close();
         return EXIT_FAILURE;
      }
      ++lineNumber;
      if (!( std::getline(inFile, line) && inFile.good()))
      {
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
         }
         inFile.close();
         return EXIT_FAILURE;
      }
   } // read atomTypeCount number of masses

   //if ( debugFlag)
   //{
   //   std::cout << "masses: ";
   //   for ( const auto& [key, value] : masses)
   //   {
   //      std::cout << "  " << key << " " << value << std::endl;
   //   }
   //}

   // read next section 
   // if line is not blank, exit 
   size_t first = line.find_first_not_of(" \t");
   if ( first != std::string::npos) // npos : max size of a string
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                               lineNumber, line);
         std::cout << "Either too many lines found in Masses section, or "
            << " the Masses section was not followed by a blank line." 
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   if ( skipLinesOfWhitespace( inFile, line) != EXIT_SUCCESS)
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath, 
                               lineNumber, line);
         std::cout << "LmpReader failed during skipLinesOfWhitespace"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   // ignore all sections that are not named 'Atoms'
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   term1.clear();
   lineStream >> term1; // trims whitespace
   while ( term1.compare("Atoms") && inFile.good())
   {
      ++lineNumber;
      if (!( std::getline(inFile, line) && inFile.good()))
      {
         // failure
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
            std::cout << "LmpReader, could not find section named Atoms"
               << std::endl;
         }
         inFile.close();
         return EXIT_FAILURE;
      }
   }
   term2.clear();
   lineStream >> term2; // style of atom output.
   }
   sectionName = "Atoms";
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   if ( skipLinesOfWhitespace( inFile, line) != EXIT_SUCCESS)
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                               lineNumber, line);
         std::cout << "LmpReader failed during skipLinesOfWhitespace"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   VectorDim tmpPosition;
   // size_t moleculeTag; 

   // read atomPopulationCount number of atomPositionType
   for ( size_t ii=0; ii < atomPopulationCount; ++ii)
   {
      // parse contents of the line
      { // scope so that lineStream can be redeclared later
      std::istringstream lineStream( line);

      lineStream >> tmpAtomID;
      lineStream >> tmpAtomType;
      // 'molecule-tag' is ignored.
      // lineStream >> moleculeTag;
      lineStream >> tmpPosition[ 0];
      lineStream >> tmpPosition[ 1];
      lineStream >> tmpPosition[ 2];
      }

      atomIDs.emplace_back( tmpAtomID);
      atomTypes.emplace_back( tmpAtomType);

      atomPositions.emplace_back(
            tmpAtomID,
            deformationMatrix * tmpPosition * scaleFactor // 1e-10/(DC->DN->poly.b_SI)
            );

      // If the whole population was read, then don't read any more lines.
      if ( ii == atomPopulationCount -1) break;
      // read another line
      ++lineNumber;
      if (!( std::getline(inFile, line) && inFile.good()))
      {
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
         }
         inFile.close();
         return EXIT_FAILURE;
      }
   }

   //if ( debugFlag)
   //{
   //   std::cout << "finished reading " << atomPopulationCount
   //      <<  " atom positions after " << lineNumber << " lines"
   //      << std::endl;
   //}
   
   return EXIT_SUCCESS;
}

int model::LmpReader::readLmpStreamBounds( std::vector<double>& bounds)
{
   // Read only the header section of the lammps file.

   // Following specification found at:
   //   https://docs.lammps.org/2001/data_format.html
   // required sections: header, Masses, Atoms

   // empty the vectors to write be written to
   if ( bounds.size() != 6) bounds.clear();

   // open the file
   std::ifstream inFile(lammpsFilePath.c_str(),std::ifstream::in);
   if (!( inFile.good()))
   {
      std::cout << "Error: failed to open the file "
        << lammpsFilePath << std::endl;
      inFile.close();
      return EXIT_FAILURE;
   }

   //std::cout<<"Reading file "<< lammpsFilePath <<"..."<<std::endl; // debug
   std::string line;
   std::stringstream ss;
   std::string sectionName;
   lineNumber = 0;
   size_t atomPopulationCount; atomPopulationCount = 0;
   size_t atomTypeCount; atomTypeCount = 0;
   double tmpDbl1, tmpDbl2; tmpDbl1 = 0; tmpDbl2 = 0;
   std::string term1, term2, term3, term4;

   // throw away the first two lines
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      std::cout << "Error: failed to read the first line of file "
         << lammpsFilePath << std::endl;
      inFile.close();
      return EXIT_FAILURE;
   }
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      std::cout << "Error: failed to read the first line of file "
         << lammpsFilePath << std::endl;
      inFile.close();
      return EXIT_FAILURE;
   }
   
   // read the header section
   sectionName = "header";
   // header: 
   //    500 atoms
   //    0 bonds     #skip, only required if there are > 0 of them
   //    0 angles    #skip, only required if there are > 0 of them
   //    0 dihedrals #skip, only required if there are > 0 of them
   //    0 impropers #skip, only required if there are > 0 of them
   //    # force field coefficients?
   //
   //    1 atom types
   //    0 bond types      #skip, only required if there are > 0 of them
   //    0 angle types     #skip, only required if there are > 0 of them
   //    0 dihedral types  #skip, only required if there are > 0 of them
   //    0 improper types  #skip, only required if there are > 0 of them
   //
   //    -1 1 xlo xhi
   //    -1 1 ylo yhi
   //    -1 1 zlo zhi

   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   lineStream >> atomPopulationCount;
   lineStream >> term2;
   }
   // match to atomPopulationCount ' atoms'
   if ( (! atomPopulationCount)
         || (term2.compare("atoms")) // 0 if they're equal
         )
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                     lineNumber, line);
         std::cout << "current line should be <size_t> atoms" << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   // skip '<size_t> {bonds,angle,dihedrals,impropers}' lines++lineNumber;
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   term1.clear();
   term2.clear();
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   lineStream >> term1;
   lineStream >> term2;
   }
   while (!( (term2.compare("bonds"))
         || (term2.compare("angles"))
         || (term2.compare("dihedrals"))
         || (term2.compare("impropers"))
         ))
   {
      //if ( debugFlag) std::cout << "skipping line " << line << std::endl;
      ++lineNumber;
      if (!( std::getline(inFile, line) && inFile.good()))
      {
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
         }
         inFile.close();
         return EXIT_FAILURE;
      }
      std::istringstream lineStream( line);
      term1.clear();
      term2.clear();
      lineStream >> term1;
      lineStream >> term2;
   }
   if ( skipLinesOfWhitespace( inFile, line) != EXIT_SUCCESS)
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                               lineNumber, line);
         std::cout << "LmpReader failed during skipLinesOfWhitespace"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   term1.clear();
   term2.clear();
   term3.clear();
   lineStream >> term1;
   lineStream >> term2;
   lineStream >> term3;
   }
   // term1,2,3 should be "<size_t> atom types"
   if ( term2.compare("atom") || term3.compare("types"))
   {
      // failure
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
         std::cout << "LmpReader, line should be '<size_t> atom types'"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( term1);
   lineStream >> atomTypeCount;
   }

   // read system x bounds
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }

   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   term1.clear();
   term2.clear();
   lineStream >> term1;
   lineStream >> term2;
   }
   while (!( (term2.compare("bonds"))
         || (term2.compare("angles"))
         || (term2.compare("dihedrals"))
         || (term2.compare("impropers"))
         ))
   { // skip '<size_t> {bonds,angle,dihedrals,impropers} types' lines
      ++lineNumber;
      if (!( std::getline(inFile, line) && inFile.good()))
      {
         if ( debugFlag)
         {
            outputCommonErrorMessage( std::cout, lammpsFilePath,
                                       lineNumber, line);
         }
         inFile.close();
         return EXIT_FAILURE;
      }
      //std::istringstream lineStream2( line);
      //term1.clear();
      //term2.clear();
      //lineStream2 >> term1;
      //lineStream2 >> term2;
   }
   if ( skipLinesOfWhitespace( inFile, line) != EXIT_SUCCESS)
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                               lineNumber, line);
         std::cout << "LmpReader failed during skipLinesOfWhitespace"
            << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   //term3.clear();
   //term4.clear();
   lineStream >> tmpDbl1;
   lineStream >> tmpDbl2;
   lineStream >> term3;
   lineStream >> term4;
   }
   // term1,2,3,4 should be "<double> <double> xlo xhi"
   if ( term3.compare("xlo") || term4.compare("xhi"))
   {
      // failure
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
         std::cout << "LmpReader, line should be "
            << "'<double> <double> xlo xhi" << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   bounds.emplace_back( tmpDbl1 * scaleFactor);
   bounds.emplace_back( tmpDbl2 * scaleFactor);

   // read system y bounds
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   lineStream >> tmpDbl1;
   lineStream >> tmpDbl2;
   lineStream >> term3;
   lineStream >> term4;
   }
   // term1,2,3,4 should be "<double> <double> ylo yhi"
   if ( term3.compare("ylo") || term4.compare("yhi"))
   {
      // failure
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
         std::cout << "LmpReader, line should be "
            << "'<double> <double> ylo yhi" << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   bounds.emplace_back( tmpDbl1 * scaleFactor);
   bounds.emplace_back( tmpDbl2 * scaleFactor);

   // read system z bounds
   ++lineNumber;
   if (!( std::getline(inFile, line) && inFile.good()))
   {
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   { // scope so that lineStream can be redeclared later
   std::istringstream lineStream( line);
   lineStream >> tmpDbl1;
   lineStream >> tmpDbl2;
   lineStream >> term3;
   lineStream >> term4;
   }
   // term1,2,3,4 should be "<double> <double> zlo zhi"
   if ( term3.compare("zlo") || term4.compare("zhi"))
   {
      // failure
      if ( debugFlag)
      {
         outputCommonErrorMessage( std::cout, lammpsFilePath,
                                    lineNumber, line);
         std::cout << "LmpReader, line should be "
            << "'<double> <double> zlo zhi" << std::endl;
      }
      inFile.close();
      return EXIT_FAILURE;
   }
   bounds.emplace_back( tmpDbl1 * scaleFactor);
   bounds.emplace_back( tmpDbl2 * scaleFactor);

   return EXIT_SUCCESS;
}

void model::LmpReader::outputCommonErrorMessage(
      std::ostream& outStream,
      const std::string& lammpsFilePath,
      const size_t& lineNumber,
      const std::string& line)
{
   outStream << "Error: failed to read "
      << lammpsFilePath << ", line " << lineNumber << ": "
      << line << std::endl;
   return;
}
#endif
