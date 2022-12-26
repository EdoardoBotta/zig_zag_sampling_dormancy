#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <random>
#include "tree.hpp"
using namespace std;

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
    void FlipVelocity(const int ind){
        const int nTimeIntervals = this->tree.CountTimeIntervals();
        
        if(ind == nTimeIntervals){
            this->c_velocity *= -1;
        }
        else if(ind == nTimeIntervals + 1){
            this->K_velocity *= -1;
        }
        else if(ind == nTimeIntervals + 2){
            this->theta1_velocity *= -1;
        }
        else if(nTimeIntervals + 3){
            this->theta2_velocity *= -1;
        }
        else {
            this->time_velocities[ind] *= -1;
        }
    }

    public:
    double GetVelocity(const int ind){
        const int nTimeIntervals = this->tree.CountTimeIntervals();
        if(ind == nTimeIntervals){
            return this->c_velocity;
        }
        else if(ind == nTimeIntervals + 1){
            return this->K_velocity ;
        }
        else if(ind == nTimeIntervals + 2){
            return this->theta1_velocity;
        }
        else if(nTimeIntervals + 3){
            return this->theta2_velocity;
        }
        else {
            return this->time_velocities[ind];
        }
    }
    

    double ComputeLambda(const float theta1, const float theta2, const float c, 
                         const float K, const int ind, 
                         const bool upperBound = false){
        const int nTimeIntervals = this->tree.CountTimeIntervals();

        if(ind == nTimeIntervals){
            return this->ComputeLambdaC(theta1, theta2, c, K, upperBound);
        }
        else if(ind == nTimeIntervals + 1){
            return this->ComputeLambdaK(c,K, upperBound);
        }
        else if(ind == nTimeIntervals + 2){
            return this->ComputeLambdaTheta(theta1, 1, upperBound);
        }
        else if(nTimeIntervals + 3){
            return this->ComputeLambdaTheta(theta2, 2, upperBound);
        }
        else {
            return this->ComputeLambdaI(theta1, theta2, c, K, ind, upperBound);
        }
    }

    double ComputeLambdaC(const double theta1, 
                          const double theta2, 
                          const double c, 
                          const double K,
                          const bool upperBound = false){
        
        Tree curTree = this->tree;
        double brackets = 0, cOffset = 0, KOffset = 0;
        double prod;
        int nDormancyPeriods;

        if(upperBound){
            cOffset = this->PlusMinus(this->c_velocity);
            KOffset = this->PlusMinus(this->K_velocity);
        }

        for(int i=0; i<curTree.CountTimeIntervals(); i++){
            auto [nAwake, nDormient] = curTree.CountLineagesUsing(i);
            brackets += (nAwake + (K + KOffset)*nDormient)*curTree.GetTime(i);
        }
        nDormancyPeriods = curTree.CountDormancyPeriods();
        prod = this->c_velocity*(brackets - 2*nDormancyPeriods/(c + cOffset));

        return std::max(prod, 0.0);
    }

    double ComputeLambdaK(const double c, const double K,
                          const bool upperBound = false){
        
        Tree curTree = this->tree;
        double weightedSum = 0, branchOffset = 0.0, thetaOffset = 0.0;
        double prod;
        int nDormancyPeriods;

        if(upperBound){
            cOffset = this->PlusMinus(this->c_velocity);
            KOffset = this->PlusMinus(this->K_velocity);
        }

        for(int i=0; i<curTree.CountTimeIntervals(); i++){
            auto [nAwake, nDormient] = curTree.CountLineagesUsing(i);
            cout << nDormient << "\n";
            weightedSum += nDormient*curTree.GetTime(i);
        }
        nDormancyPeriods = curTree.CountDormancyPeriods();
        cout << weightedSum;
        prod = this->K_velocity*((c + cOffset)*weightedSum - nDormancyPeriods/(K + KOffset));

        return std::max(prod, 0.0);
    }

    double ComputeLambdaTheta(const float theta, const int ind,
                              const bool upperBound = false){
        // ind is expected to be 1 for theta_1 and 2 for theta_2
        
        Tree curTree = this->tree;
        double velocity, branchLength, branchVelocity, sigma;
        double branchOffset = 0, thetaOffset = 0;
        bool linFilter;
        std::vector<int> TimeIntervalsUsed;

        if (ind == 1){
            velocity = this->theta1_velocity;
        }
        else velocity = this->theta1_velocity;

        for(Lineage branch : curTree){
            
            switch(ind){
                
                case 1: 
                linFilter = branch.isActive();
                break;
                
                case 2:
                linFilter = !branch.isActive();
                break;
            }
            
            if (linFilter) {

                if(upperBound){
                    branchVelocity = 0.0;
                    timeIntervalsUsed = branch.GetTimeIntervalsUsed();
                
                    for (int i: timeIntervalsUsed){
                        branchVelocity += this->time_velocities[i];
                    }

                    branchOffset = this->PlusMinus(branchVelocity);
                    
                    switch(ind){
                    
                        case 1: 
                        thetaOffset = this->PlusMinus(this->theta1_velocity);
                        break;
                        
                        case 2:
                        thetaOffset = this->PlusMinus(this->theta2_velocity);
                        break;
                    }
                }

                branchLength = curTree.GetBranchLength(branch);
                sigma += (
                        (branchLength + branchOffset)/2 
                        - branch.GetNumMutations()/(theta + thetaOffset)
                    );
            }
        }
        return std::max(velocity*sigma, 0.0);
    }

    double ComputeLambdaI(const float theta1, 
                          const float theta2, 
                          const float c, 
                          const float K, 
                          const int ind, 
                          const bool upperBound = false){

        Tree curTree = this->tree;
        double sigma, branchLength, theta, prod;
        double branchVelocity = 0.0, branchOffset = 0.0, thetaOffset = 0.0;
        int branchNumMutations;
        auto [nAwake, nDormient] = curTree.CountLineagesUsing(ind);
        float velocity = this->time_velocities[ind];
        std::vector<int> timeIntervalsUsed;

        for(Lineage branch : curTree){
            
            if(branch.isUsingTime(ind)){
                branchLength = curTree.GetBranchLength(branch);
                branchNumMutations = branch.GetNumMutations(); 

                if (branch.isActive()){
                    theta = theta1;
                } else theta = theta2;

                if(upperBound){
                    branchVelocity = 0.0;
                    timeIntervalsUsed = branch.GetTimeIntervalsUsed();
                    
                    for (int i: timeIntervalsUsed){
                        branchVelocity += this->time_velocities[i];
                    }

                    branchOffset = this->PlusMinus(branchVelocity);
                    
                    if (branch.isActive()){
                        thetaOffset = this->PlusMinus(this->theta1_velocity);
                    } else thetaOffset = this->PlusMinus(this->theta2_velocity);
                    
                }
            
                sigma += (
                            (theta + thetaOffset)/2 
                            - branchNumMutations/(branchLength + branchOffset)
                        ); 
            }           
        }

        prod = velocity*(sigma + BinomialCoefficient(nAwake, 2) + c*nAwake + c*K*nDormient);

        return std::max(prod, 0.0);
    }
};

class ZigZagSimulation{
    SimState state;
    int theta1, theta2, c, K;

    void SetParameter(int ind, double value){
        const int nTimeIntervals = this->state.tree.CountTimeIntervals();
        
        if(ind == nTimeIntervals){
            this->c = value;
        }
        else if(ind == nTimeIntervals + 1){
            this->K = value;
        }
        else if(ind == nTimeIntervals + 2){
            this->theta1 = value;
        }
        else if(nTimeIntervals + 3){
            this->theta2 = value;
        }
        else {
            this->state.tree.times[ind] = value;
        }
    }

    public:
    double GetParameter(int ind){
        const int nTimeIntervals = this->state.tree.CountTimeIntervals();
        
        if(ind == nTimeIntervals){
            return this->c;
        }
        else if(ind == nTimeIntervals + 1){
            return this->K;
        }
        else if(ind == nTimeIntervals + 2){
            return this->theta1;
        }
        else if(nTimeIntervals + 3){
            return this->theta2;
        }
        else {
            return this->state.tree.GetTime(ind);
        }
    }

    double SimulateNextFlip(const int ind, const double tMax){
        
        SimState curState = this->state;
        double rho = 0, alpha;
        double velocity = this->state.GetVelocity(ind)
        double lambdaUpper = this->state.ComputeLambda(this->theta1, this->theta2, this->c, 
                                                       this->K, ind, true);

        int nTimeIntervals = curState.tree.CountTimeIntervals();

        std::default_random_engine generator;
        std::exponential_distribution<double> distribution(lambdaUpper);

        std::mt19937 gen(1);
        std::uniform_real_distribution<double> unifDist(0, 1);

        while(true){
            rho += distribution(generator);
            if(rho < tMax){
                alpha = 0; // TODO: set new alpha accordingly
            } else alpha = 1;

            if (unifDist(gen) < alpha){
                break;
            }
        }
        
        return rho;
    }
    
    
    void run(const double tEnd, double maxIncrement){
        double t = 0, tau, nextToFlip, nextFlip, cur_value, velocity;
        Tree tree;
        int nTimeIntervals;
        // this->InitializeState(); // TODO: Implement initialization
        
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
            for(int i = 0; i < nTimeIntervals + 4; i++){
                nextFlip = this->SimulateNextFlip(i, maxIncrement);
                if (nextFlip < maxIncrement){
                    maxIncrement = nextFlip;
                    nextToFlip = i;
                }
            }
            t += maxIncrement;
            for(int i = 0; i < nTimeIntervals + 4; i++){
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