#include <mpi.h>

#include <iostream>
#include <fstream>
#include <cmath>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>

#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"

#include "RandomEngineMPIServer.hh"

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "read_card.h"

#define DIETAG      1
#define WORKTAG     2
#define RESULT_TAG  3
#define RESULTE_TAG 4
#define RESULTI_TAG 5

using namespace Garfield;
using namespace std;

typedef struct EventParam_s {
  int eventID;
  double x0;
  double y0;
  double z0;
  double t0;
  double e0;
} EventParam;

typedef struct EventResult_s {
  int eventID;

  int ne;
  int ni;

  int tupleSize;
} EventResult;

MPI_Datatype mpi_event_param_type;
MPI_Datatype mpi_event_result_type;

/* Local functions */

static void master(int argc, char** argv, int size, int nEvents);
static void slave(int myrank, int size, int random_server);
static void random_number_server(void);
void create_mpi_event_param_type(void);
void create_mpi_event_result_type(void);

int
main(int argc, char **argv)
{
  int myrank, size;

  int nEvents = 0;

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int res = read_card(argc, argv);

  if(!res) {
    std::cout << "Error while reading the card file.\n" ;
    return 1;
  }

  nEvents = param_num_events;
  if(myrank==0) 
    cout << "NUMBER OF EVENTS: " << nEvents << endl;

  if(size<3) {
    cerr << "Error: Requries at least 3 processes (1 master, 1 random_server, 1 or more slaves).\n";
    MPI_Finalize();
    return 1;
  }

  if(nEvents < size-2) {
    cerr << "Error: Number of events should not be less than number of slaves.\n";
    MPI_Finalize();
    return 1; 
  }
  
  // Create a MPI type for EventResult structure  
  create_mpi_event_param_type();

  // Create a MPI type for EventResult structure
  create_mpi_event_result_type();

  if (myrank == 0) {
    master(argc, argv, size, nEvents);
  } else if (myrank == size-1) {
    random_number_server();
  } else {
    slave(myrank, size, size-1);
  }

  /* Shut down MPI */

  MPI_Finalize();
  return 0;
}


void random_number_server(void) {
  RandomEngineMPIServer server;
  cout << "Random Server: starting the random number generator" << endl;
  server.Run();
  cout << "Random Server: random number generator killed" << endl;
}

// Dimensions of the GEM
//const double pitch = 0.014;
//const double kapton = 50.e-4;
//const double metal = 5.e-4;
//const double outdia = 70.e-4;
//const double middia = 50.e-4;


// Field map
ComponentAnsys123* fm;
// Gas
MediumMagboltz* gas;
// Sensor
Sensor* sensor;
// Microsopic avalanche
AvalancheMicroscopic* aval;
// Drift
AvalancheMC* drift;

void master(int argc, char** argv, int size, int nEvents) {
  TApplication app("app", &argc, argv);
//  plottingEngine.SetDefaultStyle();

  // Histograms
  int nBinsGain = 100;
  double gmin =   0.;
  double gmax = 100.;
  //TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons",
   //                           nBinsGain, gmin, gmax);
  //TH1F* hIons = new TH1F("hIons", "Number of ions",
    //                     nBinsGain, gmin, gmax);

  int nBinsChrg = 100;
  //TH1F* hChrgE = new TH1F("hChrgE", "Electrons on plastic",
   //                       nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);
  //TH1F* hChrgI = new TH1F("hChrgI", "Ions on plastic",
  //                        nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);


  double sumIonsTotal = 0.;
  double sumIonsDrift = 0.;
  double sumIonsPlastic = 0.;

  double sumElectronsTotal = 0.;
  double sumElectronsPlastic = 0.;
  double sumElectronsUpperMetal = 0.;
  double sumElectronsLowerMetal = 0.;
  double sumElectronsTransfer = 0.;
  double sumElectronsOther = 0.;

  vector<float> v_chrgE, v_chrgI;
  vector<double> v_xe1, v_ye1, v_ze1, v_te1, v_e1;
  vector<double> v_xe2, v_ye2, v_ze2, v_te2, v_e2;
  EventResult res;
  MPI_Status status, status2;

  ofstream out;
  out.open("progress.dat");

  // Event parameters structure
  EventParam param;

  // Load the root file containing the random energies
  //TFile *f = new TFile("driftEdep.root");
  //TTree *t = (TTree*)f->Get("drift_energy");
  //TH1F *initialE = new TH1F("initialE","initialE",100,0.0,0.1);

  TFile *f = new TFile("ntuples.root", "RECREATE");
  TNtuple* ntuple  = new TNtuple("ntuple","","xe1:ye1:ze1:te1:e1:xe2:ye2:ze2:te2:e2");
  TNtuple* ntuple1 = new TNtuple("ntuple1", "", "ne:ni");

  int iEvent=0, iCompleted=0;
  for(; iEvent<size-2; iEvent++) {
    int slave_id = iEvent + 1;
    cout << "######## SEND : " << iEvent << endl;
    
    param.eventID = iEvent;
    param.x0 = -0.0025;
    param.y0 = 0.025;
    param.z0 = 0.006;
    param.t0 = 0.;
    param.e0 = 0.1;//initialE->GetRandom();
    MPI_Send(&param, 1, mpi_event_param_type, slave_id, WORKTAG, MPI_COMM_WORLD);
  }

  while(iCompleted<nEvents) {
    
    MPI_Recv(&res, 1, mpi_event_result_type, 
      MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
    
    cout << "######## RECV : " << res.eventID << endl;

    v_xe1.resize(res.tupleSize);
    v_ye1.resize(res.tupleSize);
    v_ze1.resize(res.tupleSize);
    v_te1.resize(res.tupleSize);
    v_e1.resize(res.tupleSize);
    v_xe2.resize(res.tupleSize);
    v_ye2.resize(res.tupleSize);
    v_ze2.resize(res.tupleSize);
    v_te2.resize(res.tupleSize);
    v_e2.resize(res.tupleSize);

    MPI_Recv(&v_xe1[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);
    MPI_Recv(&v_ye1[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);
    MPI_Recv(&v_ze1[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);
    MPI_Recv(&v_te1[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);
    MPI_Recv(&v_e1[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);
    MPI_Recv(&v_xe2[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);
    MPI_Recv(&v_ye2[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);
    MPI_Recv(&v_ze2[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);
    MPI_Recv(&v_te2[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);
    MPI_Recv(&v_e2[0], res.tupleSize, MPI_DOUBLE, status.MPI_SOURCE, 
      RESULTE_TAG, MPI_COMM_WORLD, &status2);

    out << "Event ID: " << res.eventID << ", ne= " <<  res.ne << endl;

    // Accumulate the results
    //hElectrons->Fill(res.ne);
    //hIons->Fill(res.ni);

    ntuple1->Fill(res.ne,res.ni);

    for(int i=0; i<res.tupleSize; i++) {
      ntuple->Fill(v_xe1[i],v_ye1[i],v_ze1[i],v_te1[i],v_e1[i],v_xe2[i],v_ye2[i],v_ze2[i],v_te2[i],v_e2[i]);
    }
    // for(vector<float>::iterator it=v_chrgE.begin();it!=v_chrgE.end(); ++it)
    //   hChrgE->Fill(*it);

    // for(vector<float>::iterator it=v_chrgI.begin();it!=v_chrgI.end(); ++it)
    //   hChrgI->Fill(*it);

    iCompleted++;

    // Send more work (if any left) to the same slave who just finished the work
    if(iEvent<nEvents) {
      cout << "######## SEND : " << iEvent << endl;
      param.eventID = iEvent;
      param.e0 = 0.1;//initialE->GetRandom();
      MPI_Send(&param, 1, mpi_event_param_type, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
      //MPI_Send(&iEvent, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);    
      iEvent++;
    }
  }

  f->Write();

  out.close();

  int dummy=0;

  // Send DIETAG to all slaves
  for(int i=1; i<=size-2; i++)
    MPI_Send(&dummy, 1, MPI_INT, i, DIETAG, MPI_COMM_WORLD);

  // Terminate the random number server
  MPI_Send(&dummy, 1, MPI_INT, size-1, RandomEngineMPI::KILL, MPI_COMM_WORLD );

  double fFeedback = 0.;
  if (sumIonsTotal > 0.) fFeedback = sumIonsDrift / sumIonsTotal;
  std::cout << "Fraction of ions drifting back: " << fFeedback << "\n";

  //const double neMean = hElectrons->GetMean();
  //std::cout << "Mean number of electrons: " << neMean << "\n";
  //const double niMean = hIons->GetMean();
  //std::cout << "Mean number of ions: " << niMean << "\n";


/*TFile f("histos.root","new"); 
TH1F h1("hgaus","histo from a gaussian",100,-3,3); 
h1.FillRandom("gaus",10000); 
h1.Write();*/


  const bool plotHistogram = false;
  if (plotHistogram) {
    TCanvas* cH = new TCanvas("cH", "Histograms", 800, 700);
    cH->Divide(2, 2);
    cH->cd(1);
    //hElectrons->Draw();
    cH->cd(2);
    //hIons->Draw();
    cH->cd(3);
    //hChrgE->Draw();
    cH->cd(4);
    //hChrgI->Draw();
    cH->SaveAs("histogram.eps");

	TFile f("histos.root","new");
	cH->Write();

  }

  //app.Run(kTRUE);
}

void do_work(const EventParam& param, EventResult& res, 
    vector<double>& v_xe1, vector<double>& v_ye1, vector<double>& v_ze1, vector<double>& v_te1, vector<double>& v_e1,
    vector<double>& v_xe2, vector<double>& v_ye2, vector<double>& v_ze2, vector<double>& v_te2, vector<double>& v_e2);

void slave(int myrank, int size, int random_server) {

  const bool debug = true;

  // Load the field map.
  fm = new ComponentAnsys123();
  const std::string efile = param_mesh_dir + "/ELIST.lis";
  const std::string nfile = param_mesh_dir + "/NLIST.lis";
  const std::string mfile = param_mesh_dir + "/MPLIST.lis";
  const std::string sfile = param_mesh_dir + "/PRNSOL.lis";
  fm->Initialise(efile, nfile, mfile, sfile, "mm");
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  fm->PrintRange();
/*
  if(param_opt_cache_boundingboxes)
    fm->EnableCachingOfBoundingBoxes();
  else
	fm->DisableCachingOfBoundingBoxes();

  if(param_opt_search_tetratree)
    fm->EnableTetrahedralTreeForElementSearch();
  else
    fm->DisableTetrahedralTreeForElementSearch();

  if(param_opt_search_neighbors)
    fm->EnableSearchThroughNeighbors();
  else
    fm->DisableSearchThroughNeighbors();
*/
  // Dimensions of the GEM
  const double pitch = param_dim_pitch;
  const double kapton = param_dim_kapton;
  const double metal = param_dim_metal;
  const double outdia = param_dim_outdia;
  const double middia = param_dim_middia;

  // Setup the gas.
  gas = new MediumMagboltz();
  gas->SetComposition(param_gas1_name, param_gas1_comp, param_gas2_name, param_gas2_comp, param_gas3_name, param_gas3_comp); 
  gas->SetTemperature(param_gas_temperature);
  gas->SetPressure(param_gas_pressure);
  gas->EnableDebugging();
  gas->Initialise();
  //gas->DisableDebugging();
  // Set the Penning transfer efficiency.
  const double rPenning = param_gas_penning;
  const double lambdaPenning = 0.;
  gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  gas->LoadIonMobility(param_ion_filepath);

  // Associate the gas with the corresponding field map material.
  const int nMaterials = fm->GetNumberOfMaterials();
  for (int i = 0; i < nMaterials; ++i) {
    const double eps = fm->GetPermittivity(i);
    if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
  }
  fm->PrintMaterials();

  // Create the sensor.
  sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(param_sensor_x1, param_sensor_y1, param_sensor_z1,
                   param_sensor_x2,  param_sensor_y2,  param_sensor_z2);

  aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);

  drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);

  EventParam param = {0};

  EventResult res = {0};

  vector<float> v_chrgE, v_chrgI;
  vector<double> v_xe1, v_ye1, v_ze1, v_te1, v_e1;
  vector<double> v_xe2, v_ye2, v_ze2, v_te2, v_e2;

  int eventID;
  MPI_Status status;

  while(1) {
    /* Receive a message from the master */
    MPI_Recv(&param, 1, mpi_event_param_type, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

    /* Check the tag of the received message. */
    if (status.MPI_TAG == DIETAG) {
      break;
    } else if(status.MPI_TAG == WORKTAG) {
      /* Do the work */
     // v_chrgE.clear();
      //v_chrgI.clear();
      v_xe1.clear(); v_ye1.clear(); v_ze1.clear(); v_te1.clear(); v_e1.clear();
	  v_xe2.clear(); v_ye2.clear(); v_ze2.clear(); v_te2.clear(); v_e2.clear();

      res = {0};
      do_work(param, res, v_xe1, v_ye1, v_ze1, v_te1, v_e1, v_xe2, v_ye2, v_ze2, v_te2, v_e2);
      
      /* Send the result back */
      MPI_Send(&res, 1, mpi_event_result_type, 0, RESULT_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_xe1[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_ye1[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_ze1[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_te1[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_e1[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_xe2[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_ye2[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_ze2[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_te2[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_e2[0], res.tupleSize, MPI_DOUBLE, 0, RESULTE_TAG, MPI_COMM_WORLD);
    } else {
      cerr << "Slave: Unknow MPI tag" << endl;
    }
  }
}

void do_work(const EventParam& param, EventResult& res, 
    vector<double>& v_xe1, vector<double>& v_ye1, vector<double>& v_ze1, vector<double>& v_te1, vector<double>& v_e1,
    vector<double>& v_xe2, vector<double>& v_ye2, vector<double>& v_ze2, vector<double>& v_te2, vector<double>& v_e2) {
  // Randomize the initial position.
  res = {0};
  res.eventID = param.eventID;
  //const double smear = param_dim_pitch / 2.;
  //double x0 = -smear + RndmUniform() * smear;
  //double y0 = -smear + RndmUniform() * smear;
  //double x0 = param.x0;
  //double y0 = param.y0;
  //double z0 = param.z0;
  //double t0 = param.t0;
  //double e0 = param.e0;
  //aval->SetUserHandleElectronDriftPosition(userHandleElectronDriftPosition);

  double x0 = -0.0025;
  double y0 = 0.025;
  double z0 = 0.006; 
  double t0 = 0.;
  double e0 = 0.5;
  aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
  int ne = 0, ni = 0;
  aval->GetAvalancheSize(ne, ni);

  std::cout << "Avalanche Size: ne=" << ne << ", ni=" << ni << std::endl;
  //hElectrons->Fill(ne);
  //hIons->Fill(ni);
  res.ne = ne;
  res.ni = ni;
  const int np = aval->GetNumberOfElectronEndpoints();
  double xe1, ye1, ze1, te1, e1;
  double xe2, ye2, ze2, te2, e2;
  double xi1, yi1, zi1, ti1;
  double xi2, yi2, zi2, ti2;
  int status;
  for (int j = np; j--;) {
    aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1,
                                 xe2, ye2, ze2, te2, e2, status);
    v_xe1.push_back(xe1);
    v_ye1.push_back(ye1);
    v_ze1.push_back(ze1);
    v_te1.push_back(te1);
    v_e1.push_back(e1);
    v_xe2.push_back(xe2);
    v_ye2.push_back(ye2);
    v_ze2.push_back(ze2);
    v_te2.push_back(te2);
    v_e2.push_back(e2);
    // res.sumElectronsTotal += 1.;
    // if (ze2 > -kapton / 2. && ze2 < kapton / 2.) {
    //   //hChrgE->Fill(ze2 * 1.e4);
    //   v_chrgE.push_back(ze2 * 1.e4);
    //   res.sumElectronsPlastic += 1.;
    // } else if (ze2 >= kapton / 2. && ze2 <= kapton  / 2. + metal) {
    //   res.sumElectronsUpperMetal += 1.;
    // } else if (ze2 <= -kapton / 2. && ze2 >= -kapton / 2. - metal) {
    // } else if (ze2 <= -kapton / 2. && ze2 >= -kapton / 2. - metal) {
    //   res.sumElectronsLowerMetal += 1.;
    // } else if (ze2 < -kapton / 2. - metal) {
    //   res.sumElectronsTransfer += 1.;
    // } else {
    //   res.sumElectronsOther += 1.;
    // }
    drift->DriftIon(xe1, ye1, ze1, te1);
    drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1,
                             xi2, yi2, zi2, ti2, status);
    // if (zi1 < 0.01) {
    //   res.sumIonsTotal += 1.;
    //   if (zi2 > 0.01) res.sumIonsDrift += 1.;
    // }
    if (zi2 > -param_dim_kapton / 2. && zi2 < param_dim_kapton / 2.) {
      //hChrgI->Fill(zi2 * 1.e4);
      // v_chrgI.push_back(zi2 * 1.e4);
      // res.sumIonsPlastic += 1.;
    }
  }

  res.tupleSize = v_xe1.size();
  // res.size_chrgE = v_chrgE.size();
  // res.size_chrgI = v_chrgI.size();
}

void create_mpi_event_param_type(void) {
  
  EventParam param;

  const int nitems=2;
  int blocklengths[2] = {1,1};
  MPI_Datatype old_types[2] = {MPI_INT,MPI_DOUBLE};
  
  MPI_Aint displ[2];
  MPI_Get_address(&param.eventID, &displ[0]);
  MPI_Get_address(&param.e0, &displ[1]);

  for(int i=nitems-1; i>=0; i--)
    displ[i] -= displ[0];

  MPI_Type_create_struct(nitems, blocklengths, displ, old_types, &mpi_event_param_type);
  MPI_Type_commit(&mpi_event_param_type);
}

void create_mpi_event_result_type(void) {
  
  EventResult res;

  const int nitems=4;
  int blocklengths[4] = {1,1,1,1};
  MPI_Datatype old_types[4] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT};

  MPI_Aint displ[4];
  MPI_Get_address(&res.eventID, &displ[0]);
  MPI_Get_address(&res.ne, &displ[1]);
  MPI_Get_address(&res.ni, &displ[2]);
  MPI_Get_address(&res.tupleSize, &displ[3]);
  
  for(int i=3; i>=0; i--)
    displ[i] -= displ[0];

  MPI_Type_create_struct(nitems, blocklengths, displ, old_types, &mpi_event_result_type);
  MPI_Type_commit(&mpi_event_result_type);
}
