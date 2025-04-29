#ifndef TRACE_PROPERTIES_H
#define TRACE_PROPERTIES_H

#include <Rtypes.h> // For Double_t
#include "Calibration.h"

class TraceProperties {
public:
    // Member variables
    Double_t amp = 0;
    Double_t energy = 0;
    int channel = 0;
    int discard_index = 0;
    Double_t quenched_energy = 0;
    Double_t energy_err = 0;
    Double_t time_ns = 0;
    Double_t time_ps = 0;
    std::vector<double> trace = {};

    int amp_idx = 0;
    bool pileupStatus = false;
    bool ampOffsetStatus = false;

    TraceProperties();

    ~TraceProperties() = default;

    void Clear();
    void DetectPileup(double pulseDuration);
    void Process();
    void SetParameters(const std::vector<double>& ftrace, const double fchannel, const double ftime_ns, const double ftime_ps, int fdiscard_index, bool bsaveTrace);
};

#endif // TRACE_PROPERTIES_H