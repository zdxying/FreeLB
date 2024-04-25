/*CA method manager class, implementations */
#include <atomic>
#include "CAManager.h"
#include "nuc.hh"
#include "grow.hh"


template <typename T>
void CAManager<T>::init_Mushzone(int Grain_Id_)
{
	/*init mushzone of newly nucleated cell with cell centre overlap with grain centre*/
	int Id = Grains[Grain_Id_].Cell_Id;
	T Orien = Grains[Grain_Id_].Orien;
	T x_Grain = Grains[Grain_Id_].x;
	T y_Grain = Grains[Grain_Id_].y;
	T trans[2][2] = {{cos(Orien), sin(Orien)},
						 {-sin(Orien), cos(Orien)}};
	/*Find the most distant Cell centre to determine Amax*/
	T Length_max = 0;
	for (int k = 0; k < 8; k++)
	{
		int is = trav_in[k];
		int js = trav_jn[k];
		T x_Cell = Geo.getVoxel(Id)[0] + is;
		T y_Cell = Geo.getVoxel(Id)[1] + js;
		T x = trans[0][0] * (x_Cell - x_Grain) + trans[0][1] * (y_Cell - y_Grain);
		T y = trans[1][0] * (x_Cell - x_Grain) + trans[1][1] * (y_Cell - y_Grain);
		T Dist = fabs(x) + fabs(y);
		if (Dist > Length_max)
		{
			Length_max = Dist;
		}
	}
	CA[Id].A_max = Length_max * Length_max;
	CA[Id].A_min = 0;
	CA[Id].Mush_Frac = 0;
}


template <typename T>
void CAManager<T>::Simple_Phase_Transition()
{
	NewSolidCells.clear();
	NewMushyCells.clear();
	//traverse fluid cells
	int Id = 0;
	std::vector<int>::iterator iter = FluidCells.begin();
	while (iter != FluidCells.end())
	{
		Id = *iter;

		if (lbmT.pop[Id].rho < get_Temp_Liquidus(lbmC.pop[Id].rho))
		{
			Geo.getVoxel(Id).getFlag() = 1;
			CA[Id].State = 1;
			CA[Id].Mush = 1;
			NewMushyCells.emplace_back(Id);
			iter = FluidCells.erase(iter);
			Nuc_Count++;
		}
		else
		{
			iter++;
		}
	}
	std::vector<int>::iterator iter1 = NewMushyCells.begin();
	while (iter1 != NewMushyCells.end())
	{
		Id = *iter1;
		if (get_CA_Neighbor_State(Id) != 0) // at least one neighbor cell is liquid
		{
			CA[Id].State = -1; // all neighbor cells are NO LONGER LIQUID
			CA[Id].Mush = -1;
			Geo.getVoxel(Id).getFlag() = -1;
			NewSolidCells.emplace_back(Id);
			iter1 = NewMushyCells.erase(iter1);
		}
		else
		{
			iter1++;
		}
	}
	ActiveCells.insert(ActiveCells.end(), NewMushyCells.begin(), NewMushyCells.end());
	std::vector<int>::iterator iter2 = ActiveCells.begin();
	while (iter2 != ActiveCells.end())
	{
		Id = *iter2;
		if (get_CA_Neighbor_State(Id) != 0) // at least one neighbor cell is liquid
		{
			CA[Id].State = -1; // all neighbor cells are NO LONGER LIQUID
			CA[Id].Mush = -1;
			Geo.getVoxel(Id).getFlag() = -1;
			NewSolidCells.emplace_back(Id);
			iter2 = ActiveCells.erase(iter2);
		}
		else
		{
			iter2++;
		}
	}
}

/* ------------------
	nucleation
------------------ */
/*create nucleation sites with critial undercooling*/
	// 23.03.08:
	//  Actually, Sites_Num denotes max nucleation sites within the computation area,
	//	we use normal distribution to create undercooling for each cell,
	//	MAYBE this mothod simulateS the thermo fluctuation in real conditions.
	//  Then we could use the [max nucleation sites] ,namely,Sites_Num to mutiply a
	//	normal distribution coefficient, thus we get nucleation sites' number.
	//  This is what we 're going to do next: give each cell a undercooling
	//	and all the undercoolings make up a normal distribution.
	// 23.03.09:
	//  A kind of classic view:
	//   We calculate each cell's undercooling and use normal distribution to
	//   caculate the possibility of nucleation within the cell.
	//  Actually, this method is equal to the following statement:
	//   We give each cell a undercooling [DTcritical] which follows normal distribution,
	//   Then we compare the undercooling [DT] interpolated from FD/FE mesh,
	//   Assuming DT>0, if [DT]>[DTcritical], nucleation takes place.

/*nucleate at nucleation sites*/
template <typename T>
void CAManager<T>::Nucleation()
{
	NewMushyCells.clear();
	NewSolidCells.clear();
	/*traverse NucSite (pre-nucleation sites)*/
	int Id = 0; // Cell Id
				// #pragma omp parallel for private(Id)
	for (int i = 0; i < Pre_NucSites.size(); i++)
	{
		Id = Pre_NucSites.at(i);
		if (CA[Id].State == 0)
		{
			if (lbmT.pop[Id].rho < get_Temp_Liquidus(lbmC.pop[Id].rho - CA[Id].DT_Critic))
			{
				// #pragma omp critical
				{
					Nucleated(Id);
				}
			}
		}
	}
	/*update Nucleated cells*/
	std::vector<int>::iterator iter = NewMushyCells.begin();
	while (iter != NewMushyCells.end())
	{
		Id = *iter;
		if (get_CA_Neighbor_State(Id) != 0) // at least one neighbor cell is liquid
		{
			CA[Id].State = -1; // all neighbor cells are NO LONGER LIQUID
			CA[Id].Mush = -1;
			Geo.getVoxel(Id).getFlag() = -1;
			NewSolidCells.emplace_back(Id);
			iter = NewMushyCells.erase(iter);
		}
		else
		{
			iter++;
		}
	}
	ActiveCells.insert(ActiveCells.end(), NewMushyCells.begin(), NewMushyCells.end());
}

template <typename T>
T CAManager<T>::get_estimated_Max_TimeStep(T CA_TimeStep_Coeff_)
{
	T max_DT = Lat_DT_Mean_Bulk + 3 * Lat_DT_Std_Bulk;
	T max_Velo = getVelo(max_DT);
	return 1.0f / max_Velo * CA_TimeStep_Coeff_;
}

/* ------------------
	capture
------------------ */

template <typename T>
void CAManager<T>::Grain_capture()
{
	CapturedCells.clear();
	std::vector<int>::iterator iter = ActiveCells.begin();
	int Id = 0;
	int NeighbId = 0;
	int GrainId = 0;
	while (iter != ActiveCells.end())
	{
		Id = *iter; // cell id
		GrainId = CA[Id].Grain_Id;
		/*traverse neighbor cells and capture*/
		// #pragma omp parallel for private(NeighbId)
		for (int k = 1; k < 9; k++)
		{
			// NeighbId = Id + Lat.Neighbor9[k];
			// if (CA[NeighbId].State == 0)
			// {
			// 	CaptureOrNot(&Grains[GrainId], NeighbId);
			// }
		}
		// if solified then remove the grain from Grains, i.e. stop growing
		if (all_captured(Grains[GrainId]))
		{
			NewSolidCells.emplace_back(Id);
			iter = ActiveCells.erase(iter);
			Geo.getVoxel(Id).getFlag() = -1;
			CA[Id].Mush = -1;
			CA[Id].State = -1;
		}
		else
		{
			iter++;
		}
	}
	/*merge */
	ActiveCells.insert(ActiveCells.end(), CapturedCells.begin(), CapturedCells.end());
	NewMushyCells.insert(NewMushyCells.end(), CapturedCells.begin(), CapturedCells.end());
	// std::cout << "Grain = " << Grain_List.size() << std::endl;
}

/*
// void get_LocGrain_Cord(T Grain, T Cell, T Orien, T *x, T *y)
// {
// 	T trans[2][2] = {{cos(Orien), sin(Orien)},
// 						 {-sin(Orien), cos(Orien)}};
// }

// template <typename T>
// void CAManager<T>::CaptureOrNot(Grain<T> *grain, int Cell_Id)
// {
// 	bool isCaptured = false;

// 	T DENDRITE_ID = grain->Dendrite_Id;
// 	T ARMLEN[4] = {grain->ArmLen[0], grain->ArmLen[1], grain->ArmLen[2], grain->ArmLen[3]};

// 	/*get Loc cord of [cell centre] regarding [grain cetre]*/
// 	T x_Grain = grain->x;
// 	T y_Grain = grain->y;
// 	T x_Cell = Geo.getVoxel(Cell_Id)[0];
// 	T y_Cell = Geo.getVoxel(Cell_Id)[1]; // cords of cell to be captured
// 	T Orien = grain->Orien;	   // 0-90 degree
// 	T trans[2][2] = {{cos(Orien), sin(Orien)},
// 						 {-sin(Orien), cos(Orien)}}; // rotate the coordinate system counter-clockwise by an angle is equivalent to rotating the vector clockwise
// 	T x_Loc = trans[0][0] * (x_Cell - x_Grain) + trans[0][1] * (y_Cell - y_Grain);
// 	T y_Loc = trans[1][0] * (x_Cell - x_Grain) + trans[1][1] * (y_Cell - y_Grain);

// 	/*calc quadrant and choose associated ArmLen*/
// 	int quad = get_quad(x_Loc, y_Loc);

// 	T ArmLen0 = grain->ArmLen[quad - 1]; // 0 1 2 3
// 	T ArmLen1 = grain->ArmLen[quad % 4]; // 1 2 3 0

// 	/*line equation: x/a0 + y/a1 = 1 -> a1*x + a0*y = a0*a1*/
// 	if ((fabs(x_Loc * ArmLen1) + fabs(y_Loc * ArmLen0)) <= ArmLen0 * ArmLen1)
// 		isCaptured = true; // if the cell is inside the grain

// 	/*2D DECENTRED SQUARE CA GROWTH ALGORITHM: Gandin 1997*/
// 	if (isCaptured)
// 	{
// 		// Find endpoints of capture line         //  1   2   3   4
// 		T Angle0 = Orien + (quad - 1) * Pi / 2; //  0  1/2  1  3/2  Pi
// 		T Angle1 = Angle0 + Pi / 2;				// 1/2  1  3/2  0	Pi
// 		/*Arm0*/
// 		// absolute position
// 		T x_Arm0 = x_Grain + ArmLen0 * cos(Angle0);
// 		T y_Arm0 = y_Grain + ArmLen0 * sin(Angle0);
// 		/*Arm1*/
// 		// absolute position
// 		T x_Arm1 = x_Grain + ArmLen1 * cos(Angle1);
// 		T y_Arm1 = y_Grain + ArmLen1 * sin(Angle1);

// 		/*Project cell center onto capture line*/
// 		T x_Vec_Edge = x_Arm1 - x_Arm0;
// 		T y_Vec_Edge = y_Arm1 - y_Arm0; // capture line vector 0->1
// 		T x_Vec_Arm0_to_Cell = x_Cell - x_Arm0;
// 		T y_Vec_Arm0_to_Cell = y_Cell - y_Arm0;												 // vector 0->Cell
// 		T NormSquare_Vec_Edge = x_Vec_Edge * x_Vec_Edge + y_Vec_Edge * y_Vec_Edge;			 // norm square of Vec_Edge
// 		T Dot_Edge_PtoC = x_Vec_Edge * x_Vec_Arm0_to_Cell + y_Vec_Edge * y_Vec_Arm0_to_Cell; // Vec_Edge . Vec_Arm0_to_Cell
// 		// Proj: v1/|v1|*|v2|costheta = v1/|v1|*|v2|*(v1.v2)/|v1|^2 = v1*(v1.v2)/|v1|^2
// 		T x_Vec_Proj_Cell_to_Edge = x_Arm0 + x_Vec_Edge * Dot_Edge_PtoC / NormSquare_Vec_Edge; // Proj. to (1,1) face(Edge) thus devide the (1,1) into L1 and L2
// 		T y_Vec_Proj_Cell_to_Edge = y_Arm0 + y_Vec_Edge * Dot_Edge_PtoC / NormSquare_Vec_Edge; // Vec_Arm0 + proj. Vec_Arm0_to_Cell to Vec_Edge = Proj.Point

// 		// Compute line segments: Proj.Point on the Edge and divides the line into 2 parts: L1 and L2
// 		T Edge_L1 = sqrt((x_Vec_Proj_Cell_to_Edge - x_Arm0) * (x_Vec_Proj_Cell_to_Edge - x_Arm0) + (y_Vec_Proj_Cell_to_Edge - y_Arm0) * (y_Vec_Proj_Cell_to_Edge - y_Arm0));
// 		T Edge_L2 = sqrt((x_Vec_Proj_Cell_to_Edge - x_Arm1) * (x_Vec_Proj_Cell_to_Edge - x_Arm1) + (y_Vec_Proj_Cell_to_Edge - y_Arm1) * (y_Vec_Proj_Cell_to_Edge - y_Arm1));

// 		// Compute new grain size: Gandin's thesis
// 		T ArmLength_New = std::min(static_cast<T>(Edge_L1 / sqrt(2.0)), static_cast<T>(1)) + std::min(static_cast<T>(Edge_L2 / sqrt(2.0)), static_cast<T>(1));

// 		// Compute closest corner
// 		T x_Point = x_Arm0;
// 		T y_Point = y_Arm0;
// 		T ArmLen = ArmLen0;
// 		if (Edge_L2 < Edge_L1)
// 		{
// 			x_Point = x_Arm1;
// 			y_Point = y_Arm1;
// 			ArmLen = ArmLen1;
// 		}
// 		// Compute new grain center
// 		T Norm_Vec_Point = sqrt((x_Point - x_Grain) * (x_Point - x_Grain) + (y_Point - y_Grain) * (y_Point - y_Grain)); // = armLength
// 		T x_Grain_New = x_Grain + (ArmLen - ArmLength_New) * (x_Point - x_Grain) / Norm_Vec_Point;
// 		T y_Grain_New = y_Grain + (ArmLen - ArmLength_New) * (y_Point - y_Grain) / Norm_Vec_Point; // + *sin grow along Point direction

// 		// #pragma omp critical
// 		{
// 			CA[Cell_Id].State = 1;
// 			CA[Cell_Id].Mush = 1;
// 			Geo.Flag[Cell_Id] = 1;
// 			// Create new grain with same orientation and assign it to this cell
// 			int Grain_Id = NewGrain(x_Grain_New, y_Grain_New, Cell_Id, Orien);
// 			/*add to captured grain list*/
// 			CapturedCells.emplace_back(Cell_Id);
// 			/*add to dendrite*/
// 			Addto_Dendrite(DENDRITE_ID, Grain_Id, Cell_Id);

// 			/*set ArmLen[]*/
// 			T ArmLen_Coeff = ArmLength_New / ArmLen;
// 			for (int k = 0; k < 4; k++)
// 			{
// 				Grains[Grain_Id].ArmLen[k] = ArmLen_Coeff * ARMLEN[k];
// 			}

// 			set_Mushzone(Grain_Id);
// 		}
// 	}
// }


/*initilize the mushzone of a grain: calc Mush_Frac of a initially growing mushzone*/
template <typename T>
void CAManager<T>::set_Mushzone(int Grain_Id_)
{
	int Id = Grains[Grain_Id_].Cell_Id;
	T Orien = Grains[Grain_Id_].Orien;
	T x_Grain = Grains[Grain_Id_].x;
	T y_Grain = Grains[Grain_Id_].y;
	T trans[2][2] = {{cos(Orien), sin(Orien)},
						 {-sin(Orien), cos(Orien)}}; // rotate coord sys counter-clockwise
	T x_Cell0 = Geo.getVoxel(Id)[0];
	T y_Cell0 = Geo.getVoxel(Id)[1];

	/*calc LEN_MIN*/
	T x_Loc = trans[0][0] * (x_Cell0 - x_Grain) + trans[0][1] * (y_Cell0 - y_Grain);
	T y_Loc = trans[1][0] * (x_Cell0 - x_Grain) + trans[1][1] * (y_Cell0 - y_Grain);
	/*get associated quadrant and corresponding arm*/
	int quad = get_quad(x_Loc, y_Loc);
	T ArmLen0 = Grains[Grain_Id_].ArmLen[quad - 1];
	T ArmLen1 = Grains[Grain_Id_].ArmLen[quad % 4];
	T Reduct = (fabs(x_Loc * ArmLen1) + fabs(y_Loc * ArmLen0)) / (ArmLen0 * ArmLen1); // min half of the diagonal of the square

	/*calc LEN_MAX*/
	T Max_Ampli = 0;
	/*Find the most distant Cell centre to determine Amax*/
	for (int k = 0; k < 8; k++)
	{
		T x_Cell = x_Cell0 + trav_in[k];
		T y_Cell = y_Cell0 + trav_jn[k]; // celllen = 1
		T x_Loc = trans[0][0] * (x_Cell - x_Grain) + trans[0][1] * (y_Cell - y_Grain);
		T y_Loc = trans[1][0] * (x_Cell - x_Grain) + trans[1][1] * (y_Cell - y_Grain);
		/*get associated quadrant and corresponding arm*/
		int quad = get_quad(x_Loc, y_Loc);
		T ArmLen0 = Grains[Grain_Id_].ArmLen[quad - 1];
		T ArmLen1 = Grains[Grain_Id_].ArmLen[quad % 4];
		/*line equation: x/a0 + y/a1 = 1 -> a1*x + a0*y = a0*a1 */
		T Ampli = (fabs(x_Loc * ArmLen1) + fabs(y_Loc * ArmLen0)) / (ArmLen0 * ArmLen1); // Amplification factor of ArmLen
		if (Ampli > Max_Ampli)
			Max_Ampli = Ampli;
	}																																	 // minimum dendrite envelop area																			 // maximum dendrite envelop area
	T A = (Grains[Grain_Id_].ArmLen[0] + Grains[Grain_Id_].ArmLen[2]) * (Grains[Grain_Id_].ArmLen[1] + Grains[Grain_Id_].ArmLen[3]); // dendrite envelop area
	T Amin = A * Reduct * Reduct;																									 // minimum dendrite envelop area
	T Amax = A * Max_Ampli * Max_Ampli;																								 // maximum dendrite envelop area
	T Mush_frac = (A - Amin) / (Amax - Amin);
	if (Mush_frac > 1)
		Mush_frac = 1;
	/* for a given cell: Soli_Vol_Frac = MushZone_Vol_Frac * Soli_Vol_Frac_in_Envelop */
	CA[Id].A_max = Amax;
	CA[Id].A_min = Amin;
	CA[Id].DMush_Frac = Mush_frac - CA[Id].Mush_Frac;
	CA[Id].Mush_Frac = Mush_frac; // will be used in mush zone calculation,
	if (Mush_frac == 1)
		CA[Id].Mush = -1;
}

