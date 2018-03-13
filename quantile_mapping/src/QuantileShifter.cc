#include <TFile.h>
#include <math.h>
#include <iostream>
#include "../interface/QuantileShifter.h"

QuantileShifter::QuantileShifter(std::string filename, std::string source_name, std::string target_name, bool use_bisect_flag=false){
    init(filename, source_name, target_name, use_bisect_flag);
}

void QuantileShifter::init(std::string filename, std::string source_name, std::string target_name, bool use_bisect_flag=false){
    use_bisect = use_bisect_flag;
    TFile percentagefile(filename.c_str(), "READ");

    source = *((TSpline3*)(percentagefile.Get(source_name.c_str())));
    target = *((TSpline3*)(percentagefile.Get(target_name.c_str())));
    percentagefile.Close();
    initialized = true;
}

double QuantileShifter::shift(double value, double linear_interpolation_threshold) const{
    if (!initialized){
        std::cout << "QuantileShifter was not initialized!" << std::endl;
        throw std::exception();
    }
    if (linear_interpolation_threshold < 0.0){
        std::cout << "Threshold for linear interpolation must be positive!" << std::endl;
        throw std::exception();
    }
    
    double xup;
    double yup;
    double xdown;
    double ydown;
    
    int npoints = source.GetNp();
    source.GetKnot(0, xdown, ydown);
    source.GetKnot(npoints-1, xup, yup);
    if (value < xdown || value > xup){
        std::cout << "QuantileShifter - WARNING - Input value " << value << "out of range [" << xdown << ", " << xup << "]. No correction applied." << std::endl;
        return value;
    }
    
    double percentage = std::max(0.0, std::min(source.Eval(value), 1.0));
    
    if (percentage<linear_interpolation_threshold || percentage > 1.0 - linear_interpolation_threshold){
        int bin = source.FindX(value);
        if (bin >= npoints - 1){ //this can happen if an input exactly hits the upper boundary
            npoints = target.GetNp();
            target.GetKnot(npoints-1, xup, yup);
            return xup;
        }else{
            source.GetKnot(bin, xdown, ydown);
            source.GetKnot(bin + 1, xup, yup);
            percentage = ydown + (value - xdown) * (yup - ydown) / (xup - xdown);
            
            bin = FindY(percentage, true);
            target.GetKnot(bin, xdown, ydown);
            target.GetKnot(bin + 1, xup, yup);
            if (ydown == yup) return xdown;
            return xdown + (percentage - ydown) / (yup - ydown) * (xup - xdown);
        }
    }
    
    npoints = target.GetNp();
    target.GetKnot(0, xdown, ydown);
    for (int i=0; i < npoints-1; i++){
        target.GetKnot(i, xdown, ydown);
        target.GetKnot(i+1, xup, yup);
        if (percentage <= yup && yup != 0.0) break;
    }
    
    int steps = 5;
    if (use_bisect) return bisect(percentage, xup, xdown, steps);
    
    double result = xdown + (percentage - ydown) / (yup - ydown) * (xup - xdown);
    double derivative = target.Derivative(result);
    if (derivative == 0.0){
        std::cout << "QuantileShifter - WARNING - Default inversion method fails due to zero derivative. Bisecting is used instead." << std::endl;
        return bisect(percentage, xup, xdown, steps);
    }
    double correction = (percentage - target.Eval(result)) / target.Derivative(result);
    result += correction;
    if (std::abs(correction) > (xup - xdown) / 2.0 || result < xdown || result > xup){
        std::cout << "QuantileShifter - WARNING - Default inversion method yields too large corrections. Bisecting is used instead." << std::endl;
        return bisect(percentage, xup, xdown, steps);
    }
    return result;
}

double QuantileShifter::bisect(double percentage, double up, double down, int steps) const{
    if (steps <= 0) return down + (percentage - target.Eval(down)) / (target.Eval(up) - target.Eval(down)) * (up - down);
    double middle = (up + down) / 2.0;
    if (percentage > target.Eval(middle)) return bisect(percentage, up, middle, steps - 1);
    else return bisect(percentage, middle, down, steps - 1);
}

double QuantileShifter::FindY(double percentage, bool search_in_target) const{ //similar to FindX method of TSpline
    int bindown = 0;
    int binup;
    if (search_in_target) binup = target.GetNp() - 1;
    else binup = source.GetNp() - 1;
    
    while(binup - bindown > 1){
        double x;
        double y;
        int binhalf = (binup + bindown) / 2;
        if (search_in_target) target.GetKnot(binhalf, x, y);
        else source.GetKnot(binhalf, x, y);
        if (y < percentage) bindown = binhalf;
        else binup = binhalf;
    }
    return bindown;
}
