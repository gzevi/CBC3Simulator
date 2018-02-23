#ifndef SIMSTEPS_HH
#define SIMSTEPS_HH

#include "Utils.h"

float Step0_GetCharge ( TF1* ); // simStep zero, get the charge of an event
float Step0_GetNoise ( TF1* ); // simStep zero, get the noise of an event

float Step1_PreampShaper ( int, TF1*, int ); // First simStep, pre-amplifier shaper funtion
float Step1_PreampShaperDoubleHit ( int, TF1*, TF1*, int ); // Second simStep, pre-amplifier shaper for a double hit 

bool Step2_Comparator ( float, float ); // Second simStep, comparator (binary output)

bool Step3_HitDetectSampled ( bool, int ); // Third simStep, hit detect logic of CBC3 set to sampled
bool Step3_HitDetectLatched ( bool, int, bool ); // Third simStep, hit detect logic of CBC3 set to Latched

bool Step4_HIPSuppress ( bool , int & ); // Fourth simStep, HIP suppress logic of the CBC3

float Step0_GetCharge( TF1* f) {
  // Return fC (based on Landau*Gaussian)
  float charge = 2.5; // 1 MIP ~ 2.5 fC
  charge = f->GetRandom();
  return charge; 
}

float Step0_GetNoise( TF1* f) {
  // Return fC (based on Landau*Gaussian)
  //float noise = 2.5; // 1 MIP ~ 2.5 fC
  float noise = f->GetRandom();
  return noise; 
}

float Step1_PreampShaper( int time, TF1* f, int dll) {
  // Return mV at each ns (based on parametrized shape)
  // Pulse shape is longer than a clock cycle
  
  // DLL delay shifts the shape 
  int timeCorr = time - dll;
  float voltage = 0;

  // External function
  voltage = f->Eval(timeCorr);

  return voltage;
}

float Step1_PreampShaperDoubleHit( int time, TF1* f, TF1* f2, int dll ) {
  // Return mV at each ns (based on parametrized shape)
  // Pulse shape is longer than a clock cycle
  
  // DLL delay shifts the shape 
  int timeCorr = time - dll;
  float voltage = 0;

  // External function
  voltage = TMath::Max( f->Eval(timeCorr), f2->Eval(timeCorr-25));

  return voltage;
}

bool Step2_Comparator( float voltage , float threshold ) {
  // Return 0 or 1
  // Need to implement hysteresis (will requirer historical knowledge)
  return voltage > threshold;
}

bool Step3_HitDetectSampled( bool comp , int time ) {
  // Return 0 or 1 
  // Also called "40 MHz"
  // [CBC3 Manual] The output from the comparator is sampled using the 40MHz clock from the Delay Locked Loop. Only comparator outputs present on the rising edge of the clock will be captured and the output will only return to zero on the first rising clock edge following the comparators return to zero. The minimum width of output pulse is one clock cycle. Hits following immediately one after another in subsequent clock cycles will be captured even if the channel does not return below the comparator threshold for each hit.
  if (time != 0) return false;
  else return comp;
} 

bool Step3_HitDetectLatched( bool comp , int time, bool compMinusOneNs ) {
  // Also called "Fixed Width"
  // [CBC3 Manual] The output from the Hit comparator is latched for a full 25ns clock period. Non-synchronous comparator output transitions are captured. The output is a fixed 25ns pulse, regardless of the width of the comparator output pulse. Hits following immediately one after another in subsequent clock cycles will be captured provided the channel returns below the comparator threshold for each hit.
  
  // [CBC3 Manual] 2.7 ns "Setup & Hold blind spot" at the beginning of clock cycle
  if (time < 3) return false;

  // Capture transition from 0 to 1
  if ( comp && !compMinusOneNs ) return true;
  else return false;
}

bool Step4_HIPSuppress( bool hit , int &countPreviousHits, int hipSuppress ) {
  if ( countPreviousHits >= hipSuppress ) {
    return false;
  }
  else return hit;
}

#endif
