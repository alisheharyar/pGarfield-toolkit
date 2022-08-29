#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ namespace Garfield;

#pragma link C++ class Garfield::AvalancheMicroscopic;
#pragma link C++ class Garfield::AvalancheMC;
#pragma link C++ class Garfield::DriftLineRKF;

#pragma link C++ class Garfield::Medium;
#pragma link C++ class Garfield::MediumGas;
#pragma link C++ class Garfield::MediumMagboltz;
#pragma link C++ class Garfield::MediumSilicon;
#pragma link C++ class Garfield::MediumDiamond;
#pragma link C++ class Garfield::MediumGaAs;
#pragma link C++ class Garfield::MediumGaN;
#pragma link C++ class Garfield::MediumCdTe;
#pragma link C++ class Garfield::OpticalData;

#pragma link C++ class Garfield::Solid;
#pragma link C++ class Garfield::SolidBox;
#pragma link C++ class Garfield::SolidHole;
#pragma link C++ class Garfield::SolidExtrusion;
#pragma link C++ class Garfield::SolidRidge;
#pragma link C++ class Garfield::SolidSphere;
#pragma link C++ class Garfield::SolidTube;
#pragma link C++ class Garfield::SolidWire;

#pragma link C++ class Garfield::Geometry;
#pragma link C++ class Garfield::GeometrySimple;
#pragma link C++ class Garfield::GeometryRoot;

#pragma link C++ class Garfield::Component;
#pragma link C++ class Garfield::ComponentAnalyticField;
#pragma link C++ class Garfield::ComponentFieldMap;
#pragma link C++ class Garfield::ComponentAnsys123;
#pragma link C++ class Garfield::ComponentAnsys121;
#pragma link C++ class Garfield::ComponentCST;
#pragma link C++ class Garfield::ComponentComsol;
#pragma link C++ class Garfield::ComponentElmer;
#pragma link C++ class Garfield::ComponentElmer2d;
#pragma link C++ class Garfield::ComponentConstant;
#pragma link C++ class Garfield::ComponentGrid;
#pragma link C++ class Garfield::ComponentNeBem2d;
#pragma link C++ class Garfield::ComponentNeBem3d;
#pragma link C++ class Garfield::ComponentParallelPlate;
#pragma link C++ class Garfield::ComponentTcad2d;
#pragma link C++ class Garfield::ComponentTcad3d;
#pragma link C++ class Garfield::ComponentUser;

#pragma link C++ class Garfield::Sensor;

#pragma link C++ class Garfield::ViewCell;
#pragma link C++ class Garfield::ViewDrift;
#pragma link C++ class Garfield::ViewFEMesh;
#pragma link C++ class Garfield::ViewField;
#pragma link C++ class Garfield::ViewGeometry;
#pragma link C++ class Garfield::ViewIsochrons;
#pragma link C++ class Garfield::ViewMedium;
#pragma link C++ class Garfield::ViewSignal;

#pragma link C++ class Garfield::Track;
#pragma link C++ class Garfield::TrackHeed;
#pragma link C++ class Garfield::TrackSrim;
#pragma link C++ class Garfield::TrackBichsel;

#pragma link C++ function Garfield::RndmUniform();
#pragma link C++ function Garfield::SetDefaultStyle();

#endif
