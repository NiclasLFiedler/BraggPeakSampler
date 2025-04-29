#ifndef PLOTTER_H
#define PLOTTER_H

#include <TLegend.h>
#include <TPaveStats.h>
#include <TH1D.h>
#include <TString.h>

class Plotter {
private:
    bool bhetero;
    int coincidence_time;
    double tot_edep_cut_off;

public:
    Plotter(bool bhetero, int coincidence_time, double tot_edep_cut_off);
    void draw_legend(TH1D* hist1, TH1D* hist2, TH1D* hist3 = nullptr, TH1D* hist4 = nullptr);
    void draw_legend2(TH1D* hist1, TH1D* hist2);
    void plot_th1d(TH1D* hist, TString x_title, TString y_title);
};

#endif
