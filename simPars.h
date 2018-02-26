#ifndef SIMPARS_HH
#define SIMPARS_HH

// Simulation parameters
const char *    cSimFileName = "VcthOnly_WithChargeSharing45_10000events"; // simulation result file name

// Sensor properties

const double    pitch = 90.0;                 // um
const double    csRange = 45.0;                // um, range over which there is charge sharing (assumed to be in between inplants)

// Currently see no reason to modify these once they are set

const int       NEvents = 10000;              // Number of events for each RunTest call
const int       Nck = 4;                      // Number of clock cycles
const int       HipSuppress = 1;              // 0 means don't suppress

// Knobs that are overwritten when performing scans

int             Vcth = 20;                    // DAC units, already pedestal-subtracted
int             DLL = 10;                     // DLL delay. If DLL=0, signal pulse starts at t=0. If DLL=10, it starts at 10

// Knobs that are not overwritten during simulations

bool            chargeSharing = true;        // Include charge sharing or not

// Verbosity etc
int             verbose = 0;                  // 0: none, 1: every event, 2: every clock, 3: every ns
bool            diagnosticHistograms = false; // Make histograms for every event (reduce NEvents to avoid slowing down)
bool            doVcthScan = true;            // Do a 1-D Vcth scan
bool            doDLLScan = false;            // Do a 1-D Delay scan
bool            do2DScan = false;             // Do a 2-D scan over both Vcth and the Delay

#endif
