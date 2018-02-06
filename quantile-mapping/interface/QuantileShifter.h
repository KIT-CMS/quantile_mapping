#pragma once

#include <TSpline.h>


class QuantileShifter{
public:
    
    QuantileShifter(){};
    QuantileShifter(std::string filename, std::string source_name, std::string target_name, bool use_bisect_flag);
    void init(std::string filename, std::string source_name, std::string target_name, bool use_bisect_flag);

    double shift(double value) const;
    
private:
    
    TSpline3 source;
    TSpline3 target;
    bool initialized = false;
    bool use_bisect;

    double bisect(double quantile, double up, double down, int steps) const;
};
