#include "zig_zag_code.hpp"
#include <vector>
#include <assert.h>
#include <math.h>
using namespace std;

Tree CreateSampleTree(){
    Lineage lin_1, lin_2, lin_3, lin_4, lin_5, lin_6, lin_7, lin_8;

    lin_1.active = true;
    lin_1.nMutations = 0;
    lin_1.time_accesses = {0};
    lin_1.leaf = true;

    lin_2.active = true;
    lin_2.nMutations = 0;
    lin_2.time_accesses = {0};
    lin_2.leaf = true;

    lin_3.active = true;
    lin_3.nMutations = 0;
    lin_3.time_accesses = {1};
    lin_3.leftChild = &lin_1;
    lin_3.rightChild = &lin_2;
    lin_3.leaf = false;

    lin_2.parent = &lin_3;
    lin_1.parent = &lin_3;

    lin_4.active = false;
    lin_4.nMutations = 0;
    lin_4.time_accesses = {2};
    lin_4.leftChild = &lin_3;
    lin_4.leaf = false;

    lin_3.parent = &lin_4;
    
    lin_5.active = true;
    lin_5.nMutations = 0;
    lin_5.time_accesses = {3,4};
    lin_5.leftChild = &lin_4;
    lin_5.leaf = false;

    lin_4.parent = &lin_5;

    lin_6.active = true;
    lin_6.nMutations = 0;
    lin_6.time_accesses = {0,1,2,3};
    lin_6.leaf = false;

    lin_7.active = true;
    lin_7.nMutations = 0;
    lin_7.time_accesses = {0,1,2,3};
    lin_7.leaf = true;

    lin_8.active = true;
    lin_8.nMutations = 0;
    lin_8.time_accesses = {4};
    lin_8.leftChild = &lin_7;
    lin_8.rightChild = &lin_6;
    lin_8.leaf = true;

    lin_7.parent = &lin_8;
    lin_6.parent = &lin_8;

    Tree tree({lin_1, lin_2, lin_3, lin_4, lin_5, lin_6, lin_7, lin_8}, {0.1,0.2,0.3,0.2,0.4});
    return tree;
}

void testTreeIterator(){
    Tree t = CreateSampleTree();
    for (Lineage lin : t){
        cout << lin.isActive();
    }
}

void testSimState(){
    Tree t = CreateSampleTree();
    SimState ss(t, {0.2,0.2,0.2,0.2,0.2}, 0.2, 0.2, 0.2, 0.2, 5);

    assert(fabs(ss.ComputeLambdaC(1,1,1,1) - 0.26) < 1e-4);
    assert(fabs(ss.ComputeLambdaC(1,1,0.0001,1) - 0) < 1e-4);
    
    assert(fabs(ss.ComputeLambdaK(2,2) - 0.02) < 1e-4); // test lambda_K not passing (0.02)
    assert(fabs(ss.ComputeLambdaTheta(2, 1) - 0.3) < 1e-4);
    assert(fabs(ss.ComputeLambdaTheta(2, 2) - 0.03) < 1e-4);

    assert(fabs(ss.ComputeLambdaI(2, 4, 1, 1, 2) - 1.6) < 1e-4);

}

void testLikelihood(){
    Tree t = CreateSampleTree();
    double ll = t.ComputeLikelihood(1, 1);
    
    
    cout << ll;
    cout << "\n";
    assert((ll - 0.003027555) < 1e-4);
}



int main(){

    // testTreeIterator();
    testLikelihood();
    testSimState();

    /* t = create_sample_tree(); 
    
    assert(t.get_dormancy_periods() == 1);
    cout << t.get_branch_length(0);
    
    assert((t.get_branch_length(0) - 0.1) < 1e-4);
    assert((t.get_branch_length(3) - 0.3) < 1e-4);
    assert((t.get_branch_length(4) - 0.6) < 1e-4);
    assert((t.get_branch_length(6) - 0.8) < 1e-4);

    test_sim_state(); */




    return 0;
}