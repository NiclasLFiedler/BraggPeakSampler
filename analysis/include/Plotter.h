#ifndef PLOTTER_H
#define PLOTTER_H

#include <TLegend.h>
#include <TPaveStats.h>
#include <TH1D.h>
#include <TString.h>

class Plotter {
private:
    int coincidence_time;

public:
    Plotter(int coincidence_time);
    void Legend(TH1D* hist=nullptr);
    void Legend(TGraphErrors* graph=nullptr);
    void Legend(TH1D* hist, TLine* line, double bestfit, const char* param, const char* unit, double offset);
    void Histogram1D(TH1D* hist, TString x_title, TString y_title);
    void Histogram2D(TH2D* hist, TString x_title, TString y_title);
    void GraphError(TGraphErrors* graph, TString x_title, TString y_title, TString title);
    void AddVerticalLine(double x, double ymin, double ymax);
};

#endif
