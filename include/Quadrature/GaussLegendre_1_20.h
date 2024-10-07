 /* This file is part of MoDELib, the Mechanics Of Defects Evolution Library. 
 * 
 * 
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
#ifndef model_GAUSSLEGENDRE_1_20_H_ 
#define model_GAUSSLEGENDRE_1_20_H_ 

namespace model
{

   template<>
   struct GaussLegendre<1,20>
   {
       static Eigen::Matrix<double,2,20> abcsissasAndWeights()
       {
           Eigen::Matrix<double,20,2> aw;
           aw<<3.435700407452502e-03, 8.807003569579235e-03, 
               1.801403636104343e-02, 2.030071490016726e-02, 
               4.388278587433764e-02, 3.133602416765558e-02, 
               8.044151408889044e-02, 4.163837078783605e-02, 
               1.268340467699243e-01, 5.096505990862069e-02, 
               1.819731596367427e-01, 5.909726598075927e-02, 
               2.445664990245863e-01, 6.584431922458842e-02, 
               3.131469556422902e-01, 7.104805465919081e-02, 
               3.861070744291775e-01, 7.458649323623891e-02, 
               4.617367394332512e-01, 7.637669356536389e-02, 
               5.382632605667486e-01, 7.637669356536278e-02, 
               6.138929255708221e-01, 7.458649323630152e-02, 
               6.868530443577097e-01, 7.104805465919073e-02, 
               7.554335009754136e-01, 6.584431922458885e-02, 
               8.180268403632573e-01, 5.909726598075893e-02, 
               8.731659532300756e-01, 5.096505990861991e-02, 
               9.195584859111090e-01, 4.163837078835188e-02, 
               9.561172141256624e-01, 3.133602416705428e-02, 
               9.819859636389568e-01, 2.030071490019512e-02, 
               9.965642995925476e-01, 8.807003569575599e-03; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 

