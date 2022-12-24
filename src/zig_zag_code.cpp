#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <random>
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

class Lineage{
    public:
    bool active;
    int n_mutations;
    std::vector<int> time_accesses;
    Lineage* left_child;
    Lineage* right_child;
    Lineage* parent;

    public:
    bool is_active(){
        return this->active;
    }
};

class Tree{
    public:
    std::vector<Lineage> v;
    std::vector<int> times;


    public:
    int count_time_intervals(){
        return this->times.size();
    }

    public:
    int get_dormancy_periods(){
        int count_dormancy_periods = 0;
        
        for (int i=0; i < this->v.size(); i++){
            if (!v[i].active){
                count_dormancy_periods++;
            };    
        };
        return count_dormancy_periods;
    }
    
    public:
    float get_branch_length(int ind){
        int length = 0;
        std::vector<int> time_contribs = v[ind].time_accesses;
        for(int i=0; i < time_contribs.size(); i++){
            length += this->times[time_contribs[i]];
        }
        return length;
    }

    public:
    std::tuple<int,int> lineages_using_t_i(int ind){
        int n_awake = 0;
        int n_dormient = 0;
        int m = 0;

        for (int j=0; j < this->v.size(); j++){
            for (int k=0; k < this->v[j].time_accesses.size(); k++){
                if (this->v[j].time_accesses[k] == ind){
                    if (this->v[j].active){
                        n_awake++;
                    }
                    else n_dormient++;
                }
            }
        }
        return {n_awake, n_dormient};
    }


    public:
    float compute_likelihood(int c, int K){
        int n_awake, n_dormient;
        int sum_exp = 0;
        for (int i=0; i < this->times.size(); i++){
            n_awake = 0;
            n_dormient = 0;
            for (int j=0; j < this->v.size(); j++){
                int m = 0;
                for (int k=0; k < this->v[j].time_accesses.size(); k++){
                    if (this->v[j].time_accesses[k] == i){
                        if (this->v[j].active){
                            n_awake++;
                        }
                        else n_dormient++;
                    }
                    cout << n_awake;
                }
            }
            cout << "here";
            sum_exp = sum_exp + (BinomialCoefficient(n_awake, 2) + c*n_awake + c*K*n_dormient)*this->times[i];
        }
        cout << "Got here";
        return exp(-sum_exp)*pow(c, 2*this->get_dormancy_periods())*pow(K, this->get_dormancy_periods());
    }
};

class SimState{
    public:
    Tree tree;
    std::vector<double> time_velocities;
    double theta1_velocity, theta2_velocity, c_velocity, K_velocity;

    double compute_lambda_c(float theta1, float theta2, float c, float K){
        Tree cur_tree = this->tree;
        float brackets = 0;
        double prod;
        int n_dormancy_periods;

        for(int i=0; i<cur_tree.times.size(); i++){
            auto [n_awake, n_dormient] = cur_tree.lineages_using_t_i(i);
            brackets += (n_awake + K*n_dormient)*cur_tree.times[i];
        }
        n_dormancy_periods = cur_tree.get_dormancy_periods();
        prod = this->c_velocity*(brackets - 2*n_dormancy_periods/c);

        return std::max(prod, 0.0);
    }

    double compute_lambda_K(float c, float K){
        Tree cur_tree = this->tree;
        float weighted_sum = 0;
        double prod;
        int n_dormancy_periods;

        for(int i=0; i<cur_tree.times.size(); i++){
            auto [n_awake, n_dormient] = cur_tree.lineages_using_t_i(i);
            weighted_sum += n_dormient*cur_tree.times[i];
        }
        n_dormancy_periods = cur_tree.get_dormancy_periods();
        prod = this->K_velocity*(c*weighted_sum - n_dormancy_periods/K);

        return std::max(prod, 0.0);
    }

    double compute_lambda_theta(float theta, int ind){
        // ind is expected to be 1 for theta_1 and 2 for theta_2
        Tree cur_tree = this->tree;
        Lineage lineage;
        double velocity, lin_length, sigma;
        bool lin_filter;

        if (ind == 1){
            velocity = this->theta1_velocity;
        }
        else velocity = this->theta1_velocity;

        for(int i=0; i<cur_tree.v.size(); i++){
            lineage = cur_tree.v[i];
            
            switch(ind){
                
                case 1: 
                lin_filter = lineage.is_active();
                break;
                
                case 0:
                lin_filter = !lineage.is_active();
                break;
            }
            
            if (lin_filter) {
                lin_length = cur_tree.get_branch_length(i);
                sigma += (lin_length/2 - lineage.n_mutations/theta);
            }   
        }
        return std::max(velocity*sigma, 0.0);
    }

    double compute_lambda_i(float theta_1, float theta_2, float c, float K, int ind){
        Tree cur_tree = this->tree;
        Lineage lineage;
        float sigma, lin_length, theta;
        double prod;
        int lin_n_mutations;
        auto [n_awake, n_dormient] = cur_tree.lineages_using_t_i(ind);
        float velocity = this->time_velocities[ind];

        for(int i=0; i<cur_tree.v.size(); i++){
            
            lineage = cur_tree.v[i];
            lin_length = cur_tree.get_branch_length(i);

            if (lineage.is_active()){
                theta = theta_1;
            } else theta = theta_2;
            
            sigma += (lin_n_mutations/lin_length - theta/2);            
        }

        prod = velocity*(-sigma + BinomialCoefficient(n_awake, 2) + c*n_awake + c*K*n_dormient);

        return std::max(prod, 0.0);
    }



};

class ZigZagSimulation{
    SimState state;
    int theta1, theta2, c, K;

    double get_next_flip(const int ind, const double t_max){
        SimState cur_state = this->state;
        double rho = 0, alpha;
        double velocity, lambda_upper, U;

        n_time_intervals = cur_state.tree.count_time_intervals();

        // identify parameter
        // TODO: implement compute_x_upper() methods
        switch(ind){
            case n_time_intervals:
                velocity = state.c_velocity;
                lambda_upper = state.compute_lambda_c_upper();
                break;
            case n_time_intervals + 1:
                velocity = state.K_velocity;
                lambda_upper = state.compute_lambda_K_upper();
                break;
            case n_time_intervals + 2:
                velocity = state.theta_1_velocity;
                lambda_upper = state.compute_lambda_theta_1_upper();
                break;
            case n_time_intervals + 3:
                velocity = state.theta_2_velocity;
                lambda_upper = state.compute_lambda_theta_2_upper();
                break;
            default:
                velocity = state.time_velocities[ind];
                lambda_upper = state.compute_lambda_i_upper(ind);  
        }

        std::default_random_engine generator;
        std::exponential_distribution<double> distribution(lambda_upper);

        std::mt19937 gen(1)
        std::uniform_real_distribution<double> unif_dist(0, 1);

        while(true){
            rho += distribution(generator);
            if(rho < t_max){
                alpha = ; // set new alpha
            } else alpha = 1;

            if (unif_dist(gen) < alpha){
                break;
            }
        }
        
        return rho;
    }
    
    /* void run(){

    } */
         
};

int main(){
    Lineage lin;
    Tree tree;

    lin.active = true;
    lin.time_accesses = {0};
    tree.v = {lin};
    tree.times = {5};
    cout << tree.compute_likelihood(1, 1); 
    return 0;
}