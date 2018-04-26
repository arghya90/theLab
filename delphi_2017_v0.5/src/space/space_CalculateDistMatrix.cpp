/* 
 * Author : aRgo
 * Date : Aug 2017
 */

#include "space.h"
#include <fstream>

using namespace std;
bool debug_space = true;
SGrid <delphi_real> xyz_i, xyz_j;
delphi_real rad_i, rad_j, Rij2;


delphi_real delR = 0.5;	//Ang

void CDelphiSpace::CalculateDistMatrix() {

	//brute-force all-by-all distance calculation
	if (debug_space) ofstream fdist2("dist2.dat");

	for (int iv = 1; iv < iNatom; iv++) {
		xyz_i = sDelPhiPDB[iv].xyz;
		rad_i = sDelPhiPDB[iv].radius;

		for ( int jv = iv + 1; jv <= iNatom; jv++ ) {
			
			xyz_j = sDelPhiPDB[jv].xyz;
			rad_j = sDelPhiPDB[jv].radius;

			Rij2 = (rad_i + rad_j + delR) * (rad_i + rad_j + delR);
			distMatrix2[iv-1][jv-1] = calculateDistance2(xyz_i, xyz_j);

			if ( distMatrix2[iv-1][jv-1] < Rij2 ) adjMatrix[iv-1][jv-1] = true;		// they are adjecant

			if ( debug_space ) fdist2 << distMatrix2[iv-1][jv-1] << " " ;

			
		}//jv

		if ( debug_space ) fdist2 << endl;

	}// iv

	if ( debug_space ) fdist2.close();

}


delphi_real CDelphiSpace::calculateDistance2 ( SGrid <delphi_real> xyz1, SGrid < delphi_real> xyz2 ) {
	delphi_real dist2 = 0;

	dist2 +=  (xyz1.nX - xyz2.nX) * (xyz1.nX - xyz2.nX);
	dist2 +=  (xyz1.nY - xyz2.nY) * (xyz1.nY - xyz2.nY);
	dist2 +=  (xyz1.nZ - xyz2.nZ) * (xyz1.nZ - xyz2.nZ);

	return(dist2);
}

