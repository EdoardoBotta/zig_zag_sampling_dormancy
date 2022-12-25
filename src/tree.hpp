#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <random>
using namespace std;

class Lineage{
    public:
    bool active;
    int n_mutations;
    std::vector<int> time_accesses;
    Lineage* leftChild;
    Lineage* rightChild;
    Lineage* parent;

    public:
    bool isActive(){
        return this->active;
    }
};

class TreeIterator{
    public:
    using ValueType = Lineage;
    using PointerType = ValueType*;
    using ReferenceType = ValueType&;
    
    public:
    TreeIterator(PointerType ptr)
        : m_Ptr(ptr) {}

    TreeIterator&  operator++(){
        m_Ptr++;
        return *this;
    }

    TreeIterator operator++(int){
        TreeIterator iterator = *this;
        ++(*this);
        return iterator;
    }

    TreeIterator&  operator--(){
        m_Ptr--;
        return *this;
    }

    TreeIterator operator--(int){
        TreeIterator iterator = *this;
        --(*this);
        return iterator;
    }

    PointerType operator->(){
        return m_Ptr;
    }

    ReferenceType operator*(){
        return *m_Ptr;
    }

    ReferenceType operator[](int index){
        return *(m_Ptr + index);
    }

    bool operator==(const TreeIterator& other) const{
        return m_Ptr == other.m_Ptr;
    }

    bool operator==(const TreeIterator& other) const{
        return !(*this == other);
    }

    private:
    PointerType m_Ptr;
}

class Tree{
    public:
        using Iterator = TreeIterator;
    public:
    std::vector<Lineage> v;
    std::vector<double> times;

    public:
    Lineage GetBranch(int ind){
        return this->v[ind];
    }

    public:
    int CountTimeIntervals(){
        return this->times.size();
    }

    public:
    double GetTime(int ind){
        return this->times[ind];
    }

    public:
    int CountBranches(){
        return this->v.size();
    }

    public:
    int CountDormancyPeriods(){
        int nDormancyPeriods = 0;
        
        for (int i=0; i < this->CountBranches(); i++){
            if (!v[i].isActive()){
                nDormancyPeriods++;
            };    
        };
        return nDormancyPeriods;
    }
    
    public:
    float GetBranchLength(int ind){
        double length = 0;
        std::vector<int> timeContribs = v[ind].time_accesses;
        for(int i=0; i < timeContribs.size(); i++){
            length += this->times[timeContribs[i]];
        }
        return length;
    }

    public:
    std::tuple<int,int> CountLineagesUsing(int ind){
        int nAwake = 0;
        int nDormient = 0;
        int m = 0;

        for (int j=0; j < this->CountBranches(); j++){
            for (int k=0; k < this->v[j].time_accesses.size(); k++){
                if (this->v[j].time_accesses[k] == ind){
                    if (this->v[j].active){
                        nAwake++;
                    }
                    else nDormient++;
                }
            }
        }
        return {nAwake, nDormient};
    }


    public:
    float ComputeLikelihood(int c, int K){
        int nAwake, nDormient, m;
        double sumExp = 0;
        for (int i=0; i < this->CountTimeIntervals(); i++){
            nAwake = 0;
            nDormient = 0;
            for (int j=0; j < this->CountBranches(); j++){
                m = 0;
                for (int k=0; k < this->v[j].time_accesses.size(); k++){
                    if (this->v[j].time_accesses[k] == i){
                        if (this->v[j].active){
                            nAwake++;
                        }
                        else nDormient++;
                    }
                }
            }
            sumExp += (BinomialCoefficient(nAwake, 2) + c*nAwake + c*K*nDormient)*this->GetTime(i);
        }
        return exp(-sumExp)*pow(c, 2*this->CountDormancyPeriods())*pow(K, this->CountDormancyPeriods());
    }

    Iterator begin(){
        return Iterator(&(this->v[0]))
    }

    Iterator end(){
        return Iterator(&(this->v[0]) + this->v.size())
    }
};