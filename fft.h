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


#ifndef FFT_H
#define FFT_H

#include <fstream> // to read multitaper files
#include <boost/circular_buffer.hpp> // to hold recent samples
#include <fftw3.h> // fft library
#define _USE_MATH_DEFINES // pi
#include <cmath>

template <class T>
class Fft{
public:
  enum windowFunc {NONE=0, HAMMING, BLACKMAN_HARRIS, MULTITAPER};
  Fft(size_t winSize, windowFunc winf_, size_t frq_, size_t numTapers=5);
  ~Fft();
  void AddPoint(T p);
  bool Process();
  void GetPower(std::vector<T>& pow);
  void GetPhase(std::vector<T>& pha, unsigned long tidx,
                std::vector<T>* phaseShift=NULL);
private:
  boost::circular_buffer<T> buffer;
  fftw_complex *out;
  T* in;
  T* inTmp;
  T* winFunc;
  fftw_plan plan_forward;
  windowFunc winf;
  double winFuncSum;
  size_t frq;
  std::vector<T*> tapers;
  std::vector<T> tapersWeights;
  static const std::string dpssVectors;
  static const std::string dpssValues;
  size_t numTapers;
};

template<class T>
const std::string Fft<T>::dpssVectors = "dpss_E_102_5";
template<class T>
const std::string Fft<T>::dpssValues = "dpss_V_102_5";

template<class T>
Fft<T>::Fft(size_t winSize, windowFunc winf_, size_t frq_, size_t numTapers) :
  buffer(winSize), winf(winf_), frq(frq_) {
  out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ((winSize/2)+1) );
  in = (T*)fftw_malloc(sizeof(T) * winSize);
  inTmp = (T*)fftw_malloc(sizeof(T) * winSize);
  // FFTW_MEASURE: slow initialization, fast run time
  // FFTW_ESTIMATE: fast initialization, slow run time
  plan_forward = fftw_plan_dft_r2c_1d(winSize, in, out, FFTW_MEASURE);
  if (winf_ == windowFunc::MULTITAPER) {
    std::ifstream vecFile(dpssVectors);
    for (size_t i=0; i<numTapers; i++)
      tapers.push_back((T*)fftw_malloc(sizeof(T) * winSize));
    T val;
    size_t colNum = 0;
    bool eof = false;
    for (;!eof;) {
      for (size_t i=0; i<numTapers; i++) {
        vecFile>>val;
        if (!vecFile.good()) {
          eof = true;
          break;
        }
        tapers[i][colNum] = val;
      }
      colNum += 1;
    }
    std::ifstream valFile(dpssValues);
    for (size_t i=0; i<numTapers; i++) {
      valFile>>val;
      tapersWeights.push_back(val);
    }
  } else if (winf_ == windowFunc::HAMMING) {
    winFunc = (double *)fftw_malloc ( sizeof ( double ) * winSize );
    winFuncSum = 0.0;
    for (size_t i=0; i<winSize; i++) {
      winFunc[i] = 0.54 - 0.46 * cos(2.0*M_PI* float(i)/ (winSize - 1.0));
      winFuncSum += winFunc[i];
    }
  } else  if (winf_ == windowFunc::BLACKMAN_HARRIS) {
    winFunc = (double *)fftw_malloc ( sizeof ( double ) * winSize );
    winFuncSum = 0.0;
    for (size_t i=0; i<winSize; i++) {
      winFunc[i] = 0.35875 - 0.48829 * cos(2.0*M_PI* float(i)/(winSize-1.0))
                           + 0.14128 * cos(4.0*M_PI* float(i)/(winSize-1.0))
                           - 0.01168 * cos(6.0*M_PI* float(i)/(winSize-1.0));
      winFuncSum += winFunc[i];
    }
  } else {
    winFuncSum = winSize/2+1;
  }
}

template<class T>
Fft<T>::~Fft() {
  fftw_free(out);
  fftw_free(in);
  fftw_free(inTmp);
  if (winf != windowFunc::NONE)
    fftw_free(winFunc);
}

template<class T>
void Fft<T>::AddPoint(T p) {
  buffer.push_back(p);
}

template<class T>
bool Fft<T>::Process() {
  if (buffer.size() < buffer.capacity())
    return false;

  if (winf == windowFunc::MULTITAPER) {
    memcpy(inTmp, buffer.linearize(), sizeof(T)*buffer.capacity());

    const size_t nc = buffer.capacity()/2+1;
    std::vector<T> outTmpRe(nc);
    std::vector<T> outTmpIm(nc);
    for (size_t j = 0; j < nc; j++) {
      outTmpRe[j] = outTmpIm[j] = 0.0;
    }
    for (size_t i=0; i<tapers.size(); i++){
      for (size_t j=0; j<buffer.capacity(); j++)
        in[j] = inTmp[j] * tapers[i][j];
      fftw_execute(plan_forward);
      for (size_t j = 0; j < nc; j++) {
        outTmpRe[j] += out[j][0]*tapersWeights[i]/T(tapers.size());
        outTmpIm[j] += out[j][1]*tapersWeights[i]/T(tapers.size());
      }
    }
    for (size_t i = 0; i < nc; i++) {
      out[i][0] = outTmpRe[i];
      out[i][1] = outTmpIm[i];
    }
  } else if (winf != windowFunc::NONE) {
    memcpy(in, buffer.linearize(), sizeof(T)*buffer.capacity());
    for (size_t i=0; i<buffer.capacity(); i++)
      in[i] *= winFunc[i];
    fftw_execute(plan_forward);
  }

  return true;
}


template<class T>
void Fft<T>::GetPower(std::vector<T> &pow) {
  const size_t nc = buffer.capacity()/2+1;
  if (winf == windowFunc::MULTITAPER)
    for (size_t i = 0; i < nc; i++)
      pow[i] = sqrt(out[i][0]*out[i][0] + out[i][1]*out[i][1]) / 2.0;
  else if (winf != windowFunc::NONE)
    for (size_t i = 0; i < nc; i++)
      pow[i] = sqrt(out[i][0]*out[i][0] + out[i][1]*out[i][1]) / winFuncSum*2.0;
  else
    for (size_t i = 0; i < nc; i++)
      pow[i] = sqrt(out[i][0]*out[i][0] + out[i][1]*out[i][1]) / float(nc);
}

template<class T>
void Fft<T>::GetPhase(std::vector<T>& pha, unsigned long tidx,
                      std::vector<T>* phaseShift) {
  const size_t nc = buffer.capacity()/2+1;
  for (size_t i = 0; i < nc; i++ ) {
    pha[i] = atan2(out[i][1], out[i][0]) -
        (tidx % buffer.capacity()) / float(buffer.capacity()*buffer.capacity())
        * i * 2.0*M_PI;

    if (phaseShift) {
      float por = float(i)/float(nc-1);
      (*phaseShift)[i] = por*M_PI;
    }
  }
}

#endif // FFT_H
