#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;
using LatSet = D2Q9<T>;

int Ni;
int Nj;
T Cell_Len;
T RT;
int Thread_Num;
int BlockCellNx;

// physical properties
T rho_ref;    // g/mm^3
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
T Ra;         // Rayleigh number
// init conditions
Vector<T, 2> U_Ini;  // mm/s
T U_Max;

// bcs
Vector<T, 2> U_Wall;  // mm/s

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  iniReader param_reader("mpicomm.ini");
  // Thread_Num = param_reader.getValue<int>("OMP", "Thread_Num");
  // mesh
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellNx = param_reader.getValue<int>("Mesh", "BlockCellNx");
  // physical properties
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Dyna_Visc = param_reader.getValue<T>("Physical_Property", "Dyna_Visc");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
  // init conditions
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  // bcs
  U_Wall[0] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall1");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");
  tol = param_reader.getValue<T>("tolerance", "tol");

  MPI_RANK(0)

  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "Tolerance:         " << tol << "\n"
#ifdef MPI_ENABLED
            << "Running on " << mpi().getSize() << " processors\n"
#endif
            << "----------------------------------------------" << std::endl;
}

int main(int argc, char* argv[]) {
  std::uint8_t VoidFlag = std::uint8_t(1);
  std::uint8_t AABBFlag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  // Printer::Print_BigBanner(std::string("Initializing..."));

  mpi().init(&argc, &argv);

  MPI_DEBUG_WAIT

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  // use 0.5 for refined block
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len / 2, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni - 0.5) * Cell_Len), T(Nj * Cell_Len)));
  // for refined block
  AABB<T, 2> innercavity(
    Vector<T, 2>(T(Ni * Cell_Len) / BlockCellNx, T(Nj * Cell_Len) / BlockCellNx),
    Vector<T, 2>(T(Ni * Cell_Len) * (BlockCellNx - 1) / BlockCellNx,
                 T(Nj * Cell_Len) * (BlockCellNx - 1) / BlockCellNx));

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, cavity, Cell_Len, Ni / BlockCellNx);
  GeoHelper.forEachBlockCell([&](BasicBlock<T, 2>& block) {
    if (!isOverlapped(block, innercavity)) {
      block.refine();
    }
  });
  GeoHelper.CreateBlocks();
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FlagField, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity,
                 [&](FlagField& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  FlagFM.forEach(toplid, [&](FlagField& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, BBMovingWallFlag);
  });

  vtmo::ScalerWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  // velocity field
  BlockFieldManager<VectorFieldAOS<T, 2>, T, 2> VelocityFM(Geo);
  // set initial value of field
  Vector<T, 2> LatU_Wall = BaseConv.getLatticeU(U_Wall);
  VelocityFM.forEach(
    toplid, FlagFM, BBMovingWallFlag,
    [&](VectorFieldAOS<T, 2>& field, std::size_t id) { field.SetField(id, LatU_Wall); });

  vtmo::VectorWriter VecWriter("Velocity", VelocityFM);
  vtmo::vtmWriter<T, LatSet::d> NSWriter("mpicommtest", Geo, 1);
  NSWriter.addWriterSet(VecWriter);

  BlockLatticeManager<T, LatSet> NSLattice(Geo, BaseConv, VelocityFM);

  // ok
  NSLattice.getPopsFM().CommunicateAll();
  // ok
  VelocityFM.CommunicateAll();

  Timer MainLoopTimer;
  NSLattice.Communicate(MainLoopTimer());

  while (MainLoopTimer() < MaxStep){
  NSLattice.UpdateRho(MainLoopTimer(), AABBFlag, FlagFM);
  NSLattice.UpdateU(MainLoopTimer(), AABBFlag, FlagFM);

  NSLattice.Communicate(MainLoopTimer());
  std::cout << "Time: " << MainLoopTimer() << " of Rank: " << mpi().getRank() << std::endl;
  ++MainLoopTimer;
  }

  NSWriter.WriteBinary(MainLoopTimer());

  return 0;
}