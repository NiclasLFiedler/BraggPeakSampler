#ifndef ENERGY_CH_H
#define ENERGY_CH_H

#include <vector>
#include <iostream>

struct energy_ch {
    std::vector<Double_t> E;          // Energy value
    std::vector<Double_t> o_E;        // Energy uncertainty
    std::vector<Double_t> CH;  // Channel values
    std::vector<Double_t> o_CH;  // Channel uncertainties

    // Constructor to initialize values
    energy_ch() {}

    void clear() {
        E.clear();
        o_E.clear();
        CH.clear();
        o_CH.clear();
    }
};

#endif
