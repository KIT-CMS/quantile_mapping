#pragma once

#include <TSpline.h>


class QuantileShifter{
public:
    
    QuantileShifter(std::string filename, std::string source_name, std::string target_name, bool use_bisect_flag=false);

    double shift(double value);
    
private:

    
    TSpline3* source;
    TSpline3* target;
    bool use_bisect;

    double bisect(double quantile, double up, double down, int steps);
};