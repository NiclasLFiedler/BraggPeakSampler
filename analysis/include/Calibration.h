#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <vector>
#include <cmath>
#include <iostream>
#include "EnergyCh.h"
#include "DetectorProperties.h"

class Calibration {
public:
    std::vector<Double_t> slope;
    Double_t epsilon_1;
    Double_t epsilon_2;
    Double_t o_epsilon_1;
    Double_t o_epsilon_2;
    std::vector<Double_t> chi_1;
    std::vector<Double_t> chi_2;
    std::vector<Double_t> o_chi_1;
    std::vector<Double_t> o_chi_2;

    DetectorProperties* detProps = nullptr;

    Calibration(DetectorProperties* detProps)
        : detProps(detProps), slope(detProps->GetNLayers(), 0.0), chi_1(detProps->GetNLayers(), 0.0), chi_2(detProps->GetNLayers(), 0.0),
          o_chi_1(detProps->GetNLayers(), 0.0), o_chi_2(detProps->GetNLayers(), 0.0),
          epsilon_1(0), epsilon_2(0), o_epsilon_1(0), o_epsilon_2(0) {}

    const std::vector<Double_t>& getSlope() const { return slope; }
    void setSlope(const std::vector<Double_t>& new_slope) { slope = new_slope; }

    double GetVar(int channel, double e_i) {
        double ch = (e_i - epsilon_1) * (chi_1.at(channel) - chi_2.at(channel)) / (epsilon_1 - epsilon_2) + chi_1.at(channel);
        double sum1 = (epsilon_1 - epsilon_2) * ((chi_2.at(channel) - chi_1.at(channel)) + (chi_1.at(channel) - ch)) / (std::pow((chi_1.at(channel) - chi_2.at(channel)), 2)) * o_chi_1.at(channel);
        double sum2 = (epsilon_1 - epsilon_2) * (ch - chi_1.at(channel)) / (std::pow((chi_1.at(channel) - chi_2.at(channel)), 2)) * o_chi_2.at(channel);
        double sum3 = (1 + (ch - chi_1.at(channel)) / (chi_1.at(channel) - chi_2.at(channel))) * o_epsilon_1;
        double sum4 = (ch - chi_1.at(channel)) / (chi_1.at(channel) - chi_2.at(channel)) * o_epsilon_2;
        return std::pow(sum1, 2) + std::pow(sum2, 2) + std::pow(sum3, 2) + std::pow(sum4, 2);
    }

    double GetQuenchedVar(int channel, double energy) {
        return GetVar(channel, energy) / std::pow(1 - detProps->GetBirksConstant() * energy / detProps->GetLayerSizeZ(channel), 4);
    }

    double GetQuenchedEnergy(int channel, double trace_amp) {
        return epsilon_1 + slope.at(channel) * (trace_amp - chi_1.at(channel));
    }

    double GetEnergy(int channel, double QuenchedEnergy) {
        return QuenchedEnergy / (1 - detProps->GetBirksConstant() * QuenchedEnergy / detProps->GetLayerSizeZ(channel));
    }

    double FromEnergyToQuenched(int channel, double Energy) {
        return Energy / (1 + detProps->GetBirksConstant() * Energy / detProps->GetLayerSizeZ(channel));
    }

    double GetQuenchedStdDev(int channel, double Energy) {
        return sqrt(GetVar(channel, FromEnergyToQuenched(channel, Energy))) / (std::pow(1 - detProps->GetBirksConstant() * FromEnergyToQuenched(channel, Energy) / detProps->GetLayerSizeZ(channel), 2));
    }

    void energy_extrapolation(energy_ch source1, energy_ch source2){
        int j = 0;
        if(detProps->GetSimulationStatus()) return;
        if(detProps->GetReversedStatus()) j = detProps->GetNLayers()-1;
        for(int i = 0; i<source1.CH.size(); i++){
            epsilon_1 = source1.E;
            o_epsilon_1 = source1.o_E;
            epsilon_2 = source2.E;
            o_epsilon_2 = source2.o_E;
            chi_1.at(i) = source1.CH.at(std::abs(j-i));
            o_chi_1.at(i) = source1.o_CH.at(std::abs(j-i));
            chi_2.at(i) = source2.CH.at(std::abs(j-i));
            o_chi_2.at(i) = source2.o_CH.at(std::abs(j-i));
            slope.at(i) = (source2.E - source1.E) / (source2.CH.at(std::abs(j-i)) - source1.CH.at(std::abs(j-i)));
        }
        return;
    }
};
#endif