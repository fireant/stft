/*Copyright (c) 2013, Mosalam Ebrahimi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/


#include <iostream>
#include <iomanip> // setprecision
#include <chrono> // timing
#include <string.h> // memcpy
#include <string>

#include "fft.h"

//#define NDEBUG

using namespace std;
using namespace chrono;

void SignalSpectro(double (*funcp)(double), double (*funcp2)(double)) {
  cout<<fixed<<setprecision(10);
  ofstream spectroFile("spectro.csv");
  ofstream phaseFile("phase.csv");
  ofstream signalFile("signal.csv");

  size_t sampling_frq = 256; // Hz
  size_t dft_points = 102; // points

  Fft<double> fft(dft_points, Fft<double>::windowFunc::MULTITAPER, sampling_frq);

  vector<double> powers(dft_points/2+1);
  vector<double> phases(dft_points/2+1);

  float signal_power = 10.0f;
  float signal_phase = 0;
  float signal_frq = 62.0f; // Hz

  // 1 second signal
  double point;
  for (size_t i = 0; i <= sampling_frq; i++ ) {
    if (i<sampling_frq/2)
      point = (*funcp)(2.0*M_PI * float(i)/float(sampling_frq)
                          * signal_frq - signal_phase/float(sampling_frq) * signal_frq) * signal_power;
    else
      point = (*funcp2)(2.0*M_PI * float(i)/float(sampling_frq)
                          * signal_frq - signal_phase/float(sampling_frq) * signal_frq) * signal_power;

    fft.AddPoint(point);
    signalFile<<float(i)/float(sampling_frq)<<","<<point<<endl;

#ifndef NDEBUG
    auto start = std::chrono::high_resolution_clock::now();
#endif
    if (fft.Process()) {
      fft.GetPower(powers);
      fft.GetPhase(phases, i);

#ifndef NDEBUG
      auto end = high_resolution_clock::now();
      nanoseconds ns = duration_cast<nanoseconds>(end - start);
      cout<<"Elapsed nanosecs: "<<ns.count()<<endl;
#endif

      for (size_t j=0; j<powers.size(); j++) {
        cout<<i<<"\t"<<j<<"\tpower: "<<powers[j]<<"\tphase: "<<phases[j]<<endl;
        spectroFile<<powers[j]<<",";
        phaseFile<<phases[j]<<",";
      }
      spectroFile<<endl;
      phaseFile<<endl;
    }
  }
  spectroFile.close();
  phaseFile.close();
}

int main()
{
  SignalSpectro(cos,cos);

  int ret = system("matlab -nodesktop -r \"plots;quit\"");
  ret = system("reset");

  return ret;
}
