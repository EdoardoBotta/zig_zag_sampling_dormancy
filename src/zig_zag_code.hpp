#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <random>
#include "tree.hpp"
using namespace std;

int BinomialCoefficient(const int n, const int k) {
  std::vector<int> aSolutions(k);
  aSolutions[0] = n - k + 1;

  for (int i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}

/*
int sample_branch(Tree tree){
    sample branch from tree.v
    return index of the branch
} */

class SimState{
    public:
    Tree tree;
    std::vector<double> time_velocities;
    double theta1_velocity, theta2_velocity, c_velocity, K_velocity;
    double maxIncrement, timeLocalization;

    void UpdateTimeLocalization(){
        // TODO: Implement time localization computation
    }

    double PlusMinus(const double velocity){
        double timeLocalization = this->timeLocalization;

        if (velocity > 0){
            return std::max(velocity*timeLocalization, 0.0);
        } else return std::min(velocity*timeLocalization, 0.0);
    }
    
    public:
    double FlipVelocity(int ind){
        int nTimeIntervals = this->tree.CountTimeIntervals();
        
        switch(ind){
            case nTimeIntervals:
                this->c_velocity *= -1;
                break;
            case nTimeIntervals + 1:
                this->K_velocity *= -1;
                break;
            case nTimeIntervals + 2:
                this->theta1_velocity *= -1;
                break;
            case nTimeIntervals + 3:
                this->theta2_velocity *= -1;
                break;
            default:
                this->time_velocities[ind] *= -1;
                break;
        }
    }

    public:
    double GetVelocity(int ind){
        int nTimeIntervals = this->tree.CountTimeIntervals();
        switch(ind){
            case nTimeIntervals:
                return this->c_velocity;
            case nTimeIntervals + 1:
                return this->K_velocity;
            case nTimeIntervals + 2:
                return this->theta1_velocity;
            case nTimeIntervals + 3:
                return this->theta2_velocity;
            default:
                return this->time_velocities[ind];
        }
    }

    double ComputeLambda(float theta1, float theta2, float c, float K, int ind){
        int nTimeIntervals = this->tree.CountTimeIntervals();

        switch(ind){
            case nTimeIntervals:
                return this->ComputeLambdaC(theta1, theta2, c, K);
            case nTimeIntervals + 1:
                return this->ComputeLambdaK(c,K);
            case nTimeIntervals + 2:
                return this->ComputeLambdaTheta(theta1, 1);
            case nTimeIntervals + 3:
                return this->ComputeLambdaTheta(theta2, 2);
            default:
                return this->ComputeLambdaI(theta1, theta2, c, K, ind);
        }
    }

    double ComputeLambdaC(float theta1, float theta2, float c, float K){
        Tree curTree = this->tree;
        float brackets = 0;
        double prod;
        int nDormancyPeriods;

        for(int i=0; i<curTree.CountTimeIntervals(); i++){
            auto [nAwake, nDormient] = curTree.CountLineagesUsing(i);
            brackets += (nAwake + K*nDormient)*curTree.GetTime(i);
        }
        nDormancyPeriods = curTree.CountDormancyPeriods();
        prod = this->c_velocity*(brackets - 2*nDormancyPeriods/c);

        return std::max(prod, 0.0);
    }

    double ComputeLambdaK(float c, float K){
        Tree curTree = this->tree;
        float weightedSum = 0;
        double prod;
        int nDormancyPeriods;

        for(int i=0; i<curTree.CountTimeIntervals(); i++){
            auto [nAwake, nDormient] = curTree.CountLineagesUsing(i);
            cout << nDormient << "\n";
            weightedSum += nDormient*curTree.GetTime(i);
        }
        nDormancyPeriods = curTree.CountDormancyPeriods();
        cout << weightedSum;
        prod = this->K_velocity*(c*weightedSum - nDormancyPeriods/K);

        return std::max(prod, 0.0);
    }

    double ComputeLambdaTheta(float theta, int ind){
        // ind is expected to be 1 for theta_1 and 2 for theta_2
        Tree curTree = this->tree;
        Lineage lineage;
        double velocity, linLength, sigma;
        bool linFilter;

        if (ind == 1){
            velocity = this->theta1_velocity;
        }
        else velocity = this->theta1_velocity;

        for(int i=0; i<curTree.CountBranches(); i++){
            lineage = curTree.v[i];
            
            switch(ind){
                
                case 1: 
                linFilter = lineage.isActive();
                break;
                
                case 0:
                linFilter = !lineage.isActive();
                break;
            }
            
            if (linFilter) {
                linLength = curTree.GetBranchLength(i);
                sigma += (linLength/2 - lineage.n_mutations/theta);
            }   
        }
        return std::max(velocity*sigma, 0.0);
    }

    double ComputeLambdaI(float theta1, 
                          float theta2, 
                          float c, 
                          float K, 
                          int ind, 
                          bool upperBound = false){

        Tree curTree = this->tree;
        Lineage lineage;
        float sigma, linLength, theta;
        double prod;
        int linNumMutations;
        auto [nAwake, nDormient] = curTree.CountLineagesUsing(ind);
        float velocity = this->time_velocities[ind];

        // TODO: Should only consider brenaches using time interval i
        for(int i=0; i<curTree.CountBranches(); i++){
            
            lineage = curTree.v[i];
            linLength = curTree.GetBranchLength(i); 

            if (lineage.isActive()){
                theta = theta1;
            } else theta = theta2;
            
            sigma += (linNumMutations/linLength - theta/2);            
        }

        prod = velocity*(-sigma + BinomialCoefficient(nAwake, 2) + c*nAwake + c*K*nDormient);

        return std::max(prod, 0.0);
    }



};

class ZigZagSimulation{
    SimState state;
    int theta1, theta2, c, K;

    double SetParameter(int ind, double value){
        int nTimeIntervals = this->tree.CountTimeIntervals();
        
        switch(ind){
            case nTimeIntervals:
                this->c = value;
                break;
            case nTimeIntervals + 1:
                this->K = value;
                break;
            case nTimeIntervals + 2:
                this->theta1 = value;
                break;
            case nTimeIntervals + 3:
                this->theta2 = value;
                break;
            default:
                this->state.tree.times[ind] = value;
                break;
        }
    }

    public:
    double GetParameter(int ind){
        int nTimeIntervals = this->tree.CountTimeIntervals();
        
        switch(ind){
            case nTimeIntervals:
                return this->c;
            case nTimeIntervals + 1:
                return this->K;
            case nTimeIntervals + 2:
                return this->theta1;
            case nTimeIntervals + 3:
                return this->theta2;
            default:
                return this->state.tree.GetTime(ind);
        }
    }

    double SimulateNextFlip(const int ind, const double tMax){
        SimState curState = this->state;
        double rho = 0, alpha;
        double velocity, lambdaUpper;

        nTimeIntervals = curState.tree.CountTimeIntervals();

        // identify parameter
        // TODO: implement compute_x_upper() methods
        switch(ind){
            case nTimeIntervals:
                velocity = state.c_velocity;
                lambda_upper = state.compute_lambda_c_upper();
                break;
            case nTimeIntervals + 1:
                velocity = state.K_velocity;
                lambda_upper = state.compute_lambda_K_upper();
                break;
            case nTimeIntervals + 2:
                velocity = state.theta_1_velocity;
                lambda_upper = state.compute_lambda_theta_1_upper();
                break;
            case nTimeIntervals + 3:
                velocity = state.theta_2_velocity;
                lambda_upper = state.compute_lambda_theta_2_upper();
                break;
            default:
                velocity = state.time_velocities[ind];
                lambdaUpper = state.compute_lambda_i_upper(ind);  
        }

        std::default_random_engine generator;
        std::exponential_distribution<double> distribution(lambda_upper);

        std::mt19937 gen(1)
        std::uniform_real_distribution<double> unif_dist(0, 1);

        while(true){
            rho += distribution(generator);
            if(rho < tMax){
                alpha = ; // TODO: set new alpha
            } else alpha = 1;

            if (unif_dist(gen) < alpha){
                break;
            }
        }
        
        return rho;
    }
    
    
    void run(const double tEnd, const double maxIncrement){
        double t = 0, tau, nextToFlip, nextFlip, cur_value;
        
        this->InitializeState(); // TODO: Implement initialization
        
        while(t < tEnd){
            tau = maxIncrement;
            nextToFlip = 0;
            tree = this->state.tree;
            nTimeIntervals = tree.CountTimeIntervals();
            for(int i = 0; i < nTimeIntervals + 4; i++){
                if (this->state.GetVelocity(i) < 0){
                    // TODO: Implement flip rule
                }
            }
            for(i = 0; i < nTimeIntervals + 4; i++){
                nextFlip = this->SimulateNextFlip(i, maxIncrement);
                if (nextFlip < maxIncrement){
                    maxIncrement = nextFlip;
                    nextToFlip = i;
                }
            }
            t += maxIncrement;
            for(i = 0; i < nTimeIntervals + 4; i++){
                velocity = this->state.GetVelocity(i);
                cur_value = this->GetParameter(i);
                this->SetParameter(i, cur_value + velocity*maxIncrement);
            }
            if (nextToFlip != 0){
                this->state.FlipVelocity(nextToFlip);
            }

            /* TODO: Implement boundary conditions
               TODO: Implement MH step */
        }
    } 
         
};