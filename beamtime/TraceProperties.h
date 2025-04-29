#ifndef TRACE_PROPERTIES_H
#define TRACE_PROPERTIES_H

#include <Rtypes.h> // For Double_t

class TraceProperties {
public:
    // Member variables
    Double_t amp = 0;
    Double_t energy = 0;
    Double_t quenched_energy = 0;
    Double_t energy_err = 0;
    Double_t time_ns = 0;
    Double_t time_ps = 0;
    std::vector<double> trace = {};

    int amp_idx = 0;
    bool det_pileup = false;

    TraceProperties();
    TraceProperties(Double_t amp, Double_t energy, Double_t quenched_energy, Double_t energy_err, Double_t time_ns, Double_t time_ps, int amp_idx, bool det_pileup, const std::vector<double>& trace_data = {});

    ~TraceProperties() = default;

    // Member functions
    void Clear();
    void DetectPileup(double pulseDuration);
    void Process();
    void SetTrace(const std::vector<double>& trace);
};

#endif // TRACE_PROPERTIES_H