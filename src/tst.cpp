#include "zig_zag_code.hpp"
#include <vector>
#include <cassert>
using namespace std;

Tree CreateSampleTree(){
    Lineage lin_1, lin_2, lin_3, lin_4, lin_5, lin_6, lin_7, lin_8;
    Tree tree;

    lin_1.active = true;
    lin_1.n_mutations = 0;
    lin_1.time_accesses = {0};

    lin_2.active = true;
    lin_2.n_mutations = 0;
    lin_2.time_accesses = {0};

    lin_3.active = true;
    lin_3.n_mutations = 0;
    lin_3.time_accesses = {1};
    lin_3.left_child = &lin_1;
    lin_3.right_child = &lin_2;
    
    lin_2.parent = &lin_3;
    lin_1.parent = &lin_3;

    lin_4.active = false;
    lin_4.n_mutations = 0;
    lin_4.time_accesses = {2};
    lin_4.left_child = &lin_3;

    lin_3.parent = &lin_4;
    
    lin_5.active = true;
    lin_5.n_mutations = 0;
    lin_5.time_accesses = {3,4};
    lin_5.left_child = &lin_4;

    lin_4.parent = &lin_5;

    lin_6.active = true;
    lin_6.n_mutations = 0;
    lin_6.time_accesses = {0,1,2,3};

    lin_7.active = true;
    lin_7.n_mutations = 0;
    lin_7.time_accesses = {0,1,2,3};

    lin_8.active = true;
    lin_8.n_mutations = 0;
    lin_8.time_accesses = {4};
    lin_8.left_child = &lin_7;
    lin_8.right_child = &lin_6;
    
    lin_7.parent = &lin_8;
    lin_6.parent = &lin_8;

    tree.v = {lin_1, lin_2, lin_3, lin_4, lin_5, lin_6, lin_7, lin_8};
    tree.times = {0.1,0.2,0.3,0.2,0.4};

    return tree;
}

void test_sim_state(){
    Tree t;
    SimState ss;
    
    t = CreateSampleTree();
    ss.tree = t;
    ss.time_velocities = {0.2,0.2,0.2,0.2,0.2};
    ss.theta1_velocity = 0.2;
    ss.theta2_velocity = 0.2;
    ss.c_velocity = 0.2;
    ss.K_velocity = 0.2;

    assert((ss.compute_lambda_c(1,1,1,1) - 0.26) < 1e-4);
    assert((ss.compute_lambda_c(1,1,0.0001,1) - 0) < 1e-4);
    
    ss.compute_lambda_K(2,2); // test lambda_K not passing (0.02)

}



int main(){
    Tree t;
    SimState ss;
    double ll;

    t = create_sample_tree();
    ll = t.compute_likelihood(1, 1); 
    
    assert(t.get_dormancy_periods() == 1);
    cout << t.get_branch_length(0);
    
    assert((t.get_branch_length(0) - 0.1) < 1e-4);
    assert((t.get_branch_length(3) - 0.3) < 1e-4);
    assert((t.get_branch_length(4) - 0.6) < 1e-4);
    assert((t.get_branch_length(6) - 0.8) < 1e-4);
    
    assert((ll - 0.003027555) < 1e-4);

    test_sim_state();




    return 0;
}