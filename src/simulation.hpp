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

    SimState(Tree t, std::vector<double> tv, double theta1v, double theta2v,
             double cv, double Kv, double maxI, double timeLocalization) 
        : tree(t)
        , time_velocities(tv)
        , theta1_velocity(theta1v)
        , theta2_velocity(theta2v)
        , c_velocity(cv)
        , K_velocity(Kv)
        , maxIncrement(maxI)
        , timeLocalization(timeLocalization)
        {}
    
    void SetTimeLocalization(double timeLoc){
        this->timeLocalization = timeLoc;
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
        
        LineageCountResult countResult;
        Tree curTree = this->tree;
        double brackets = 0, cOffset = 0, KOffset = 0, tOffset = 0;
        double prod;
        int nDormancyPeriods;

        if(upperBound){
            cOffset = this->PlusMinus(this->c_velocity);
            KOffset = this->PlusMinus(this->K_velocity);
        }

        for(int i=0; i<curTree.CountTimeIntervals(); i++){
            countResult = curTree.CountLineagesUsing(i);
            nAwake = countResult.nAwake;
            nDormient = countResult.nDormient;

            if(upperBound){
                tOffset = this->PlusMinus(this->time_velocities[i]);
            }

            brackets += (nAwake + (K + KOffset)*nDormient)*(curTree.GetTime(i) + tOffset);
        }
        nDormancyPeriods = curTree.CountDormancyPeriods();
        prod = this->c_velocity*(brackets - 2*nDormancyPeriods/(c + cOffset));

        return std::max(prod, 0.0);
    }

    double ComputeLambdaK(const double c, const double K,
                          const bool upperBound = false){
        
        LineageCountResult countResult;
        Tree curTree = this->tree;
        double weightedSum = 0, cOffset = 0.0, KOffset = 0.0, tOffset = 0.0;
        double prod;
        int nDormancyPeriods;

        if(upperBound){
            cOffset = this->PlusMinus(this->c_velocity);
            KOffset = this->PlusMinus(this->K_velocity);
        }

        for(int i=0; i<curTree.CountTimeIntervals(); i++){
            countResult = curTree.CountLineagesUsing(i);
            nAwake = countResult.nAwake;
            nDormient = countResult.nDormient;

            if(upperBound){
                tOffset = this->PlusMinus(this->time_velocities[i]);
            }

            weightedSum += nDormient*(curTree.GetTime(i) + tOffset);
        }
        nDormancyPeriods = curTree.CountDormancyPeriods();
        prod = this->K_velocity*((c + cOffset)*weightedSum - nDormancyPeriods/(K + KOffset));

        return std::max(prod, 0.0);
    }

    double ComputeLambdaTheta(const float theta, const int ind,
                              const bool upperBound = false){
        // ind is expected to be 1 for theta_1 and 2 for theta_2
        
        Tree curTree = this->tree;
        double velocity, branchLength, branchVelocity;
        double branchOffset = 0, thetaOffset = 0, sigma = 0;
        bool linFilter;
        std::vector<int> timeIntervalsUsed;

        if (ind == 1){
            velocity = this->theta1_velocity;
        }
        else velocity = this->theta2_velocity;

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
        double sigma = 0, branchLength, theta, prod;
        double branchVelocity = 0.0, branchOffset = 0.0, thetaOffset = 0.0, cOffset, KOffset;
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
                    cOffset = this->PlusMinus(this->c_velocity);
                    KOffset = this->PlusMinus(this->K_velocity);
                    
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

        prod = velocity*(sigma + BinomialCoefficient(nAwake, 2) 
                            + (c + cOffset)*nAwake 
                            + (c + cOffset)*(K + KOffset)*nDormient
                        );

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

    private: Lineage SampleMHCandidateBranch(){
        int nCandidates = 0;
        int cumSum = 0;
        double randUnif;
        
        for (Lineage branch: this->state.tree){
            if (branch.isActive() || !branch.isOverlapping()){
                nCandidates += branch.GetTimeIntervalsUsed().size();
            }
        
        randUnif = GenerateRandomUniform(0, nCandidates);
        for (Lineage branch: this->state.tree){
            if (branch.isActive() || !branch.isOverlapping()){
                if (randUnif < cumSum + branch.GetTimeIntervalsUsed().size()){
                    return branch
                }
                cumSum += branch.GetTimeIntervalsUsed().size();
            }
        }    
    }

    private: void MHStep(){
        Lineage targetBranch;
        double accRate;
        double randUnif = GenerateRandomRealUniform();

        targetBranch = SampleMHCandidateBranch(); 

        move = BuildMHMove(targetBranch); // to be implemented

        accRate = ComputeAccRate(move); // to be implemented

        if (randUnif < accRate){
            return move.newTree()
        }
        return move.oldTree()
    }

    public:
    double ComputeTimeLocalization(){
        int nTimeIntervals = this->state.tree.CountTimeIntervals();
        double minLoc = -1;
        double curLoc;
        double velocity;
        int boundaryFactor;
        SimState curState = this->state;

        for(int i = 0; i < nTimeIntervals + 4; i++){
            boundaryFactor = 0;

            if(curState.tree.GetMutationsShorterChildAt(i) > 0){
                
                if(!curState.tree.IsSwapValid(i)){
                    boundaryFactor += 1;
                }
                
                if(i == 0 || i >= nTimeIntervals){
                    boundaryFactor += 1;
                }
            }

            velocity = curState.GetVelocity(i);
            if (velocity < 0){
                curLoc = this->GetParameter(i) / ((1 + 4*boundaryFactor)*velocity);
                if(curLoc < minLoc || minLoc == -1){
                    minLoc = curLoc;
                }
            }

            return std::min(minLoc, curState.maxIncrement);
        }
    }

    SimState ProgressTree(const double time){
        SimState curState = this->state;
        Tree newTree = curState.tree.ProgressTreeTimes(curState.time_velocities, time);

        return SimState(newTree, curState.time_velocities,
                        curState.theta1_velocity, curState.theta2_velocity,
                        curState.c_velocity, curState.K_velocity, curState.maxIncrement);
    }

    double SimulateNextFlip(const int ind, const double tMax){
        
        SimState curState = this->state;
        double rho = 0, alpha, progressedLambda;
        double velocity = this->state.GetVelocity(ind);
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
                SimState progressedState = this->ProgressTree(rho);
                progressedLambda = progressedState.ComputeLambda(theta1 + curState.theta1_velocity*rho,
                                                                 theta2 + curState.K_velocity*rho,
                                                                 c + curState.c_velocity*rho,
                                                                 K + curState.K_velocity*rho,
                                                                 ind, false); 
                alpha = progressedLambda/lambdaUpper; 
            } else alpha = 1;

            if (unifDist(gen) < alpha){
                break;
            }
        }
        
        return rho;
    }
    
    
    void run(const double tEnd, double maxIncrement){
        double t = 0, tau, nextToFlip, nextFlip, cur_value, velocity;
        int c, J, d, I;
        int nTimeIntervals, timeLoc;
        Tree tree({Lineage()}, {0.1});
        // this->InitializeState(); // TODO: Implement initialization
        
        while(t < tEnd){
            tau = maxIncrement;
            nextToFlip = -1;
            tree = this->state.tree;
            nTimeIntervals = tree.CountTimeIntervals();
            for(int i = 0; i < nTimeIntervals + 4; i++){
                velocity = this->state.GetVelocity(i);
                if (velocity < 0){
                    d = 1;
                    J = i;
                    if((i == 0 || i >= nTimeIntervals || this->state.tree.IsSwapValid(i)) && this->state.tree.GetMutationsShorterChildAt(i) > 0){
                        d = 1 + 4; // c = 4
                        J = 0;
                    }
                    if(-this->GetParameter(i)/(d*velocity) < tau){
                        tau = -this->GetParameter(i)/(d*velocity);
                        I = J;
                    }
                }
            }

            for(int i = 0; i < nTimeIntervals + 4; i++){
                timeLoc = this->ComputeTimeLocalization();
                this->state.SetTimeLocalization(timeLoc);
                nextFlip = this->SimulateNextFlip(i, tau);
                if (nextFlip < tau){
                    tau = nextFlip;
                    nextToFlip = i;
                }
            }
            t += tau;
            for(int i = 0; i < nTimeIntervals + 4; i++){
                velocity = this->state.GetVelocity(i);
                cur_value = this->GetParameter(i);
                this->SetParameter(i, cur_value + velocity*tau);
            }
            if (nextToFlip != -1){
                this->state.FlipVelocity(nextToFlip);
            }

            if(nextToFlip > 0 && nextToFlip < nTimeIntervals - 1){
                if(tree.IsSwapValid(nextToFlip)){
                    tree.SwapDown(nextToFlip);
                }
            }

            /* TODO: Implement MH step */
        }
    } 
         
};