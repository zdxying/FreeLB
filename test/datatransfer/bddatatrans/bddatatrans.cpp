
#include "ca/zhu_stefanescu2d.h"
#include "ca/zhu_stefanescu2d.hh"
#include "freelb.h"
#include "freelb.hh"
// Known bugs: Segmentation fault may occur, but running rhe executable again
// may resolve this without re-compile
//  this may be caused by parallel

// int Total_Macro_Step = 0;
using T = FLOAT;
using LatSet = D2Q9<T>;

/*----------------------------------------------
                        Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
T Cell_Len;
T RT;

T TimeStep;
int Thread_Num;
int BlockCellNx;
/*---------------------
        physical param
---------------------*/
/*physical property*/
T rho_ref = T(1);  // g/mm^3
T Kine_Visc;       // mm^2/s kinematic viscosity of the liquid
/*init conditions*/
Vector<T, 2> U_Ini;  // mm/s
T U_Max;

/*bcs*/
T Temp_Wall;          // K
T Conc_Wall;          // wt.%
Vector<T, 2> U_Wall;  // mm/s
/*---------------------
        LB param
---------------------*/
T Th;  // char high temp
T Tl;  // char low temp
T Ch;  // char high conc
T Cl;  // char low conc

// Simulation settings
int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  /*reader*/
  iniReader param_reader("bddatatrans.ini");
  /*mesh*/
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  /*CA mesh*/
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellNx = param_reader.getValue<int>("Mesh", "BlockCellNx");

  Kine_Visc = param_reader.getValue<T>("Phys_Prop", "Kine_Visc");

  /*init conditions*/
  U_Ini[0] = param_reader.getValue<T>("ICs", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("ICs", "U_Ini1");
  U_Max = param_reader.getValue<T>("ICs", "U_Max");
  /*bcs*/
  U_Wall[0] = param_reader.getValue<T>("BCs", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("BCs", "Velo_Wall1");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  /*output to console*/
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------" << std::endl;
}

int main() {
  std::uint8_t voidFlag = std::uint8_t(1);
  std::uint8_t AABBFlag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // ------------------ define converters ------------------
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConvertFromViscosity(Ni, U_Max, Kine_Visc);

  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  // use 0.5 for refined block
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len / 2, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni - 0.5) * Cell_Len), T(Nj * Cell_Len)));
  // use for refined block
  AABB<T, 2> outercavity(
    Vector<T, 2>(T(Ni * Cell_Len) / BlockCellNx, T(Nj * Cell_Len) / BlockCellNx),
    Vector<T, 2>(T(Ni * Cell_Len) * (BlockCellNx - 1) / BlockCellNx,
                 T(Nj * Cell_Len) * (BlockCellNx - 1) / BlockCellNx));
  AABB<T, 2> innercavity(
    Vector<T, 2>(T(Ni * Cell_Len) * (BlockCellNx / 2 - 1) / BlockCellNx,
                 T(Nj * Cell_Len) * (BlockCellNx / 2 - 1) / BlockCellNx),
    Vector<T, 2>(T(Ni * Cell_Len) * (BlockCellNx / 2 + 1) / BlockCellNx,
                 T(Nj * Cell_Len) * (BlockCellNx / 2 + 1) / BlockCellNx));

  // geometry helper
  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, Ni / BlockCellNx, cavity, Cell_Len);
  GeoHelper.CreateBlocks();
  GeoHelper.AdaptiveOptimization(Thread_Num);

  // geometry
  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FlagField, T, 2> FlagFM(Geo, voidFlag);
  FlagFM.forEach(cavity,
                 [&](FlagField& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  FlagFM.forEach(toplid, [&](FlagField& field, std::size_t id) {
    if (field.get(id) == BouncebackFlag) field.SetField(id, BBMovingWallFlag);
  });

  // ------------------ define lattice ------------------
  // velocity field
  BlockFieldManager<VectorFieldAOS<T, 2>, T, 2> VelocityFM(Geo);

  // lbm
  BlockLatticeManager<T, LatSet> NSLattice(Geo, BaseConv, VelocityFM);

  std::vector<T> RefThold;
  std::vector<T> CoaThold;
  DynamicBlockLatticeHelper2D<T, LatSet> NSDynLatHelper(NSLattice, GeoHelper, VelocityFM,
                                                        RefThold, CoaThold, 1);

  VelocityFM.forEach(FlagFM, BBMovingWallFlag,
                     [&](auto& field, std::size_t id) { field.SetField(id, U_Wall); });

  vtmo::ScalerWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag", Geo, 1);
  GeoWriter.addWriterSet(&FlagWriter);
  GeoWriter.WriteBinary();

  // --------------------- BCs ---------------------
  // NS
  BBLikeFixedBlockBdManager<T, LatSet, BounceBackLikeMethod<T, LatSet>::normal_bounceback>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, voidFlag);
  BBLikeFixedBlockBdManager<T, LatSet,
                            BounceBackLikeMethod<T, LatSet>::movingwall_bounceback>
    NS_BBMW("NS_BBMW", NSLattice, FlagFM, BBMovingWallFlag, voidFlag);

  // writer
  // 0
  vtmo::ScalerWriter RhoWriter("Rho", NSLattice.getRhoFM());
  vtmo::VectorWriter VecWriter("Velocity", VelocityFM);
  vtmo::vtmWriter<T, 2> MainWriter("beforetrans", Geo, 1);
  MainWriter.addWriterSet(&RhoWriter, &VecWriter);
  // 1
  vtmo::ScalerWriter RhoWriter1("Rho", NSLattice.getRhoFM());
  vtmo::VectorWriter VecWriter1("Velocity", VelocityFM);
  vtmo::vtmWriter<T, 2> MainWriter1("aftertrans", Geo, 1);
  MainWriter1.addWriterSet(&RhoWriter1, &VecWriter1);

  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  MainWriter.WriteBinary(MainLoopTimer());

  while (MainLoopTimer() < MaxStep) {
    NSLattice.UpdateRho(MainLoopTimer(), AABBFlag, FlagFM);
    NSLattice.UpdateU(MainLoopTimer(), AABBFlag, FlagFM);
    NSLattice.BGK<Equilibrium<T, LatSet>::SecondOrder>(
      MainLoopTimer(), std::uint8_t(AABBFlag | BouncebackFlag | BBMovingWallFlag),
      FlagFM);

    NSLattice.Stream(MainLoopTimer());

    NS_BB.Apply(MainLoopTimer());
    NS_BBMW.Apply(MainLoopTimer());

    NSLattice.Communicate(MainLoopTimer());

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      // Velocity and Conc Field Communication for output
      VelocityFM.CommunicateAll();
      NSLattice.getRhoFM().CommunicateAll();

      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());

      if (MainLoopTimer() == (MaxStep / 2)) {
        GeoHelper.forEachBlockCell([&](BasicBlock<T, 2>& block) {
          // if (isOverlapped(block, innercavity)) {
            if (!isOverlapped(block, outercavity)) {
            block.refine();
          }
        });
        GeoHelper.CreateBlocks();
        GeoHelper.AdaptiveOptimization(Thread_Num);
        // Geo Init
        Geo.Init(GeoHelper);

        // field reconstruction
        FlagFM.Init(voidFlag);
        FlagFM.forEach(cavity, [&](FlagField& field, std::size_t id) {
          field.SetField(id, AABBFlag);
        });
        FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
        FlagFM.forEach(toplid, [&](FlagField& field, std::size_t id) {
          if (field.get(id) == BouncebackFlag) field.SetField(id, BBMovingWallFlag);
        });

        // field data transfer
        VelocityFM.InitAndComm(GeoHelper);
        VelocityFM.forEach(FlagFM, BBMovingWallFlag, [&](auto& field, std::size_t id) {
          field.SetField(id, U_Wall);
        });

        NSDynLatHelper.PopFieldInit();

        NSLattice.Init();

        // Bcs init
        NS_BB.Init();
        NS_BBMW.Init();

        VecWriter.Init(VelocityFM);
        RhoWriter.Init(NSLattice.getRhoFM());
        MainWriter.Init();
        MainWriter.addWriterSet(&RhoWriter, &VecWriter);

        VecWriter1.Init(VelocityFM);
        RhoWriter1.Init(NSLattice.getRhoFM());
        MainWriter1.Init();
        MainWriter1.addWriterSet(&RhoWriter1, &VecWriter1);

        MainWriter1.WriteBinary(MainLoopTimer());

        FlagWriter.Init(FlagFM);
        GeoWriter.Init();
        GeoWriter.addWriterSet(&FlagWriter);
        GeoWriter.WriteBinary(MainLoopTimer());

        //
        // NSLattice.UpdateRho(0, AABBFlag, FlagFM);
        // NSLattice.UpdateU(0, AABBFlag, FlagFM);
        // NSLattice.BGK<Equilibrium<T, LatSet>::SecondOrder>(
        //   0, std::uint8_t(AABBFlag | BouncebackFlag | BBMovingWallFlag), FlagFM);

        // NSLattice.Stream(0);

        // NS_BB.Apply(0);
        // NS_BBMW.Apply(0);

        // NSLattice.Communicate(0);

        //          
        
      }
    }
  }

  MainWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();
  return 0;
}

// attention! DO NOT call std::srand(std::time(0)) in a while loop
// If the while loop doesn't take more than a second to complete,
// then time(0) will return the same value every time we call srand(time(0)),
// which means we are setting the same seed for rand() repeatedly.