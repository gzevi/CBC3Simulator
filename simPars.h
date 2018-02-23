#ifndef SIMPARS_HH
#define SIMPARS_HH

// Simulation parameters
const char * cSimFileName = "VcthOnly_10000events";

// Currently see no reason to modify these once they are set
const int       NEvents = 10000; // Number of events for each RunTest call
const int       Nck = 4; // Number of clock cycles
const int       HipSuppress = 1; // 0 means don't suppress

// Knobs that are overwritten when performing scans
int       Vcth = 20; // DAC units, already pedestal-subtracted
int       DLL = 10; // DLL delay. If DLL=0, signal pulse starts at t=0. If DLL=10, it starts at 10

// Verbosity etc
int       verbose = 0; // 0: none, 1: every event, 2: every clock, 3: every ns
bool      diagnosticHistograms = false; // Make histograms for every event (reduce NEvents to avoid slowing down)
bool      doVcthScan = true;
bool      doDLLScan = false;
bool      do2DScan = false;

#endif
