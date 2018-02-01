#include "quantile_shifter.h"

QuantileShifter:QuantileShifter(std::string filename, std::string source_name, std::string target_name, bool use_bisect_flag=false){
    use_bisect = use_bisect_flag;
    TFile quantilefile(inputfile.c_str(), "READ");
    source = *(quantilefile.Get(source_name.c_str()));
    target = *(quantilefile.Get(target_name.c_str()));
    quantilefile.Close();
}

double QuantileShifter:shift(double value){
    double *xup;
    double *yup;
    double *xdown;
    double *ydown;
    
    int npoints = target.GetNp()
    source.GetKnot(0, xdown, ydown)
    source.GetKnot(npoints-1, xup, yup)
    if (value < *xdown || value > *xup){
        std::cout << "QuantileShifter - WARNING - Input value " << value << "out of range [" << *xdown << ", " << *xup << "]. No correction applied." << std::endl;
        return value;
    }
    
    double quantile = std::max(0.0, std::min(source.Eval(value), 1.0))
    target.GetKnot(0, xdown, ydown)
    for (int i; i < npoints-1; i++){
        target.GetKnot(index, xdown, ydown)
        target.GetKnot(index+1, xup, yup)
        if (quantile <= *yup && *yup != 0.0) break;
    }
    
    int steps = 5
    if (use_bisect) return bisect(quantile, *xup, *xdown, steps);
    
    double result = *xdown + (quantile - *ydown) / (*yup - *ydown) * (*xup - *xdown);
    double derivative = target.Derivative(result);
    if (derivative == 0.0){
        std::cout << "QuantileShifter - WARNING - Default inversion method fails due to zero derivative. Bisecting is used instead." << std::endl;
        return bisect(quantile, *xup, *xdown, steps);
    }
    double correction = (quantile - target.Eval(result)) / target.Derivative(result);
    result += correction;
    if (std::abs(correction) > (*xup - *xdown) / 2.0 || result < *xdown || result > *xup){
        std::cout << "QuantileShifter - WARNING - Default inversion method yields too large corrections. Bisecting is used instead." << std::endl;
        return bisect(quantile, *xup, *xdown, steps);
    }
    return result;
}

double QuantileShifter:bisect(double quantile, double up, double down, int steps){
    if (steps <= 0) return down + (quantile - target.Eval(down)) / (target.Eval(up) - target.Eval(down)) * (up - down);
    double middle = (up + down) / 2.0;
    if (quantile > target.Eval(middle)) return bisect(quantile, up, middle, steps - 1);
    else return bisect(quantile, middle, down, steps - 1);
}