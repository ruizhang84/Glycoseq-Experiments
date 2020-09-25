#ifndef ENGINE_LEARN_NEURAL_NETWORK_H
#define ENGINE_LEARN_NEURAL_NETWORK_H

#include <numeric>
#include <cmath>
#include <vector>
#include <iostream>

namespace engine{
namespace learn{

class Classifier {
public:
    Classifier(std::vector<double> weight, double bias) :
        w_(weight), b_(bias) {}

    std::vector<double> Weight() { return w_; }
    double Bias() { return b_; }

    void set_weight(std::vector<double> weight) { w_ = weight; }
    void set_bias(double bias) { b_ = bias; }
    void set_debug(bool debug) { debug_ = debug; }
    
    const double Sigmoid (const double z) const 
    {
        return 1.0 / (1.0 + exp(-z));
    }
    
    const double Matmul (const std::vector<double>& x) const
    {
        return std::inner_product(x.begin(), x.end(),
            w_.begin(), b_);
    }

    const std::vector<double> Linear (const std::vector<std::vector<double>>& X) const
    {
        std::vector<double> y_pred;
        for (const auto& x : X)
        {
            y_pred.push_back(Sigmoid(Matmul(x)));
        }
        return y_pred;
    }

    const double BCELoss (const std::vector<std::vector<double>>& X, 
        const std::vector<double>& y) const
    {
        std::vector<double> y_pred = Linear(X);
        double loss = 0;
        int size = (int) y.size();
        for (int i = 0; i < size; i++)
        {
            loss += y[i] * log(y_pred[i]) + (1 - y[i]) * log(1- y_pred[i]);
        }

        return -loss * 1.0 / size; 
    }

    const std::vector<std::vector<double>> Transpose(const std::vector<std::vector<double>>& X) const
    {
        std::vector<std::vector<double>> X_T;
        for (int i = 0; i < (int) X[0].size(); i++)
        {
            std::vector<double> t;
            for (int j = 0; j < (int) X.size(); j++)
            {
                t.push_back(X[j][i]);
            }
            X_T.push_back(t);
        }
        return X_T;
    }

    const std::vector<double> Gradient(const std::vector<std::vector<double>>& X, 
        const std::vector<double>& y) const
    {
        std::vector<double> grad;
        std::vector<double> y_pred = Linear(X);
        int size = (int) y.size();
        for (int i=0; i < size; i++)
        {
            y_pred[i] -= y[i];
        }
        const std::vector<std::vector<double>> X_T = Transpose(X);
        for(const auto& it : X_T)
        {
            grad.push_back(std::inner_product(it.begin(), it.end(),
                y_pred.begin(), 0.0) * 1.0 / size);
        }
        // bias
        grad.push_back(std::accumulate(y_pred.begin(), y_pred.end(), 0.0) * 1.0 / size);
        return grad;
    }

    void Train(const std::vector<std::vector<double>>& X, 
        const std::vector<double>& y, int steps, double lr) 
    {
        for(int i=0; i < steps; i++)
        {
            if (debug_)
            {
                double loss = BCELoss(X, y);
                std::cout << "Iterations: " << i <<" loss: " << loss << std::endl;
            }
            std::vector<double> grad = Gradient(X, y);
            for(int i=0; i < (int) grad.size() - 1; i++)
            {
                w_[i] -= lr * grad[i];
            }
            b_ -= lr * grad[grad.size()-1];
        }
    }

protected:
    std::vector<double> w_;
    double b_;
    bool debug_ = false;
};


} // namespace engine
} // namespace learn

#endif