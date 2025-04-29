#ifndef ENERGY_CH_H
#define ENERGY_CH_H

#include <vector>
#include <iostream>

struct energy_ch {
    Double_t E;          // Energy value
    Double_t o_E;        // Energy uncertainty
    std::vector<Double_t> CH;  // Channel values
    std::vector<Double_t> o_CH;  // Channel uncertainties

    // Constructor to initialize values
    energy_ch(Double_t e = 0, Double_t o_e = 0)
        : E(e), o_E(o_e) {}

    void clear() {
        E = 0;
        o_E = 0;
        CH.clear();
        o_CH.clear();
    }
};

#endif
