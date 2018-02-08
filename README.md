# CBC3Simulator

Run: root -b -q CBC3Sim.C

Parameters (at the top of CBC3Sim.C):
```
// Currently see no reason to modify these once they are set
const int       NEvents = 100; // Number of events for each RunTest call
const int       Nck = 4; // Number of clock cycles
const int       HipSuppress = 1; // 0 means don't suppress

// Knobs that are overwritten when performing scans
int       Vcth = 20; // DAC units, already pedestal-subtracted
int       DLL = 10; // DLL delay. If DLL=0, signal pulse starts at t=0. If DLL=10, it starts at 10

// Verbosity etc
int       verbose = 0; // 0: none, 1: every event, 2: every clock, 3: every ns
bool      diagnosticHistograms = false; // Make histograms for every event (reduce NEvents to avoid slowing down)
bool      doVcthScan = true;
bool      doDLLScan = true;
bool      do2DScan = true;
```

The main function is CBC3sim(), which can scan parameters and call RunTest() for every configuration.

RunTest() runs over NEvents and produces efficiency plots for each readout mode (Sampled, Latched, OR, SampledHIP, ORHIP)

Currently setup for three scans:
- Vcth scan
- DLL scan
- Vcth/DLL 2D scan




