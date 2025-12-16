#ifndef ENERGY_CH_H
#define ENERGY_CH_H

#include <vector>
#include <iostream>

struct energy_ch {
    std::vector<Double_t> E = std::vector<Double_t>(32, 0);          // Energy value
    std::vector<Double_t> o_E = std::vector<Double_t>(32, 0);        // Energy uncertainty
    std::vector<Double_t> CH = std::vector<Double_t>(32, 0);  // Channel values
    std::vector<Double_t> o_CH = std::vector<Double_t>(32, 0);  // Channel uncertainties

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
