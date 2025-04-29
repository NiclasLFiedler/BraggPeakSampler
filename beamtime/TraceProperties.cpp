#include "TraceProperties.h"

TraceProperties::TraceProperties()
    : amp(0), energy(0), energy_err(0), time_ns(0), time_ps(0),
      amp_idx(0), det_pileup(false), trace({}) {}

TraceProperties::TraceProperties(Double_t amp, Double_t energy, Double_t energy_err,
                                 Double_t time_ns, Double_t time_ps, int amp_idx, bool det_pileup, const std::vector<double>& trace)
    : amp(amp), energy(energy), energy_err(energy_err), time_ns(time_ns), time_ps(time_ps),
      amp_idx(amp_idx), det_pileup(det_pileup), trace(trace) {}

void TraceProperties::Clear() {
    amp = 0;
    energy = 0;
    energy_err = 0;
    amp_idx = 0;
    time_ns = 0;
    time_ps = 0;
    det_pileup = false;
    trace = {};
}

void TraceProperties::DetectPileup(double pulseDuration) {
    bool pileupDetected = false;
    double previousPeakTime = -1;

    for (int i = 1; i < trace.size() - 1; ++i) {
        if (trace[i] < trace[i-6] && trace[i] < trace[i+6] && (trace[i-6]-trace[i]) > 25 && (trace[i+6]-trace[i]) > 25){
            if (previousPeakTime != -1.0) {
                double timeSeparation = i - previousPeakTime;
                if (timeSeparation > pulseDuration) {
                    pileupDetected = true;
                    break;
                }
            }
            previousPeakTime = i;
        }
    }
    det_pileup = pileupDetected;
    return;
}

void TraceProperties::Process(){
    Double_t sum = 0;
    int pregate = 25;
    std::vector<Double_t>::iterator max_trace;
    for (int i = 0; i < pregate; ++i){
        sum += trace.at(i);
    }
    max_trace = std::min_element(trace.begin(), trace.end());
    
    amp = static_cast<Double_t>(-*max_trace+sum/pregate);
    amp_idx = static_cast<int>(std::distance(trace.begin(), max_trace));
    DetectPileup(20);
    if(amp_idx > 50) det_pileup = true;

    return;
}

void TraceProperties::SetTrace(const std::vector<double>& ftrace){
  trace=ftrace;
  return;
}
