 /* This file is part of MoDELib, the Mechanics Of Defects Evolution Library. 
 * 
 * 
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
#ifndef model_GAUSSLEGENDRE_1_6_H_ 
#define model_GAUSSLEGENDRE_1_6_H_ 

namespace model
{

   template<>
   struct GaussLegendre<1,6>
   {
       static Eigen::Matrix<double,2,6> abcsissasAndWeights()
       {
           Eigen::Matrix<double,6,2> aw;
           aw<<3.376524289842381e-02, 8.566224618958604e-02, 
               1.693953067668680e-01, 1.803807865240662e-01, 
               3.806904069584016e-01, 2.339569672863475e-01, 
               6.193095930415985e-01, 2.339569672863454e-01, 
               8.306046932331321e-01, 1.803807865240689e-01, 
               9.662347571015757e-01, 8.566224618958555e-02; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 

