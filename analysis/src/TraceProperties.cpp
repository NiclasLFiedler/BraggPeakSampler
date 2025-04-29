#include "../include/TraceProperties.h"

TraceProperties::TraceProperties(){}

void TraceProperties::Clear(){
    amp = 0;
    energy = 0;
    quenched_energy = 0;
    channel = 0;
    energy_err = 0;
    amp_idx = 0;
    time_ns = 0;
    time_ps = 0;
    discard_index = 0;
    pileupStatus = false;
    ampOffsetStatus = false;
    trace = {};
}

void TraceProperties::DetectPileup(double pulseDuration) {
    double previousPeakTime = -1;
    for (int i = 1; i < trace.size() - 1; ++i) {
        if (trace[i] < trace[i-6] && trace[i] < trace[i+6] && (trace[i-6]-trace[i]) > 50 && (trace[i+6]-trace[i]) > 50){
            if (previousPeakTime != -1.0) {
                double timeSeparation = i - previousPeakTime;
                if (timeSeparation > pulseDuration) {
                    pileupStatus = true;
                    break;
                }
            }
            previousPeakTime = i;
        }
    }
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
    if(amp_idx > discard_index) ampOffsetStatus = true;
    DetectPileup(30);
    return;
}

void TraceProperties::SetParameters(const std::vector<double>& ftrace, const double fchannel, const double ftime_ns, const double ftime_ps, int fdiscard_index, bool bsaveTrace){
    trace=ftrace;
    discard_index = fdiscard_index;
    time_ns = ftime_ns;
    time_ps = ftime_ps;
    channel = fchannel;
    Process();
    if(!bsaveTrace){
        trace = {};
    }
    return;
}
