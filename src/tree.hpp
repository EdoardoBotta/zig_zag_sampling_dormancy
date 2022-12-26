#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <random>
#include <unordered_map>
#include "utils.hpp"
using namespace std;

class Lineage{
    public:
    bool active;
    bool leaf;
    int nMutations;
    std::vector<int> time_accesses;
    Lineage* leftChild;
    Lineage* rightChild;
    Lineage* parent;

    public:
    bool isActive(){
        return this->active;
    }

    bool isLeaf(){
        return this->leaf;
    }

    bool isOverlapping(){
        return (this->time_accesses.size() != 1);
    }

    bool isUsingTime(const int ind){
        for(int i : this->time_accesses){
            if(i == ind){
                return true;
            }
        };
        return false;
    }

    int GetNumMutations(){
        return this->nMutations;
    }

    std::vector<int> GetTimeIntervalsUsed(){
        return this->time_accesses;
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

    TreeIterator& operator++(){
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

    bool operator!=(const TreeIterator& other) const{
        return !(*this == other);
    }

    private:
    PointerType m_Ptr;
};

class Tree{
    public:
        using Iterator = TreeIterator;
    public:
    std::vector<Lineage> v;
    std::vector<double> times;
    unordered_map<int, std::vector<Lineage>> endTimes;
    

    public:
    Tree(std::vector<Lineage> branches, std::vector<double> time_ints) 
        : v(branches)
        , times(time_ints)
        , endTimes(InititalizeTimeMap(branches))
        {}
    
    private:
    unordered_map<int, std::vector<Lineage>> InititalizeTimeMap(std::vector<Lineage> branches){
        int maxValue;
        unordered_map<int, std::vector<Lineage>> endMap;
        
        for(Lineage branch : branches){
            maxValue = -1;
            for(int i : branch.time_accesses){
                if(i > maxValue){
                    maxValue = i;
                }
            }
            if (endMap.find(maxValue) == endMap.end()){
                endMap[maxValue] = {branch};
            }
            else endMap[maxValue].push_back(branch);
        }

        return endMap;
    }
        
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
        
        for (Lineage branch : this->v){
            if (!branch.isActive()){
                nDormancyPeriods++;
            };    
        }
        return nDormancyPeriods;
    }
    
    public:
    float GetBranchLength(Lineage branch){
        double length = 0;
        std::vector<int> timeContribs = branch.time_accesses;
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

        for(Lineage branch : this->v){
            if(branch.isUsingTime(ind)){
                if(branch.isActive()){
                    nAwake++;
                } 
                else nDormient++;
            }
        }

        return {nAwake, nDormient};
    }


    public:
    float ComputeLikelihood(int c, int K){
        int nAwake, nDormient, m;
        double sumExp = 0;
    
        for (int i=0; i < this->CountTimeIntervals(); i++){
            auto [nAwake, nDormient] = this->CountLineagesUsing(i);
            sumExp += (BinomialCoefficient(nAwake, 2) + c*nAwake + c*K*nDormient)*this->GetTime(i);
        }
        
        return exp(-sumExp)*pow(c, 2*this->CountDormancyPeriods())*pow(K, this->CountDormancyPeriods());
    }

    public:
    Tree ProgressTreeTimes(std::vector <double> time_velocities, 
                            const double time){
        
        if (time_velocities.size() != this->times.size()){
            throw std::invalid_argument("sizes do not match");
        }                        

        std::vector<double> new_times = this->times;

        for(int i=0; i < time_velocities.size(); i++){
            new_times[i] += time_velocities[i]*time;
        }

        return Tree(this->v, new_times);
    }

    public:
    Tree PushDown(const int ind){
        if(ind == 0 || ind == this->CountTimeIntervals()){
            throw std::invalid_argument("Invalid index for pushdown");
        }
        std::vector<Lineage>& endingAtIndex = this->endTimes[ind];
        /* 
        find branches ending at ind-1.
        If any of these is child of endingAtIndex:
            for branch in endingAtIndex:
                drop 3 from time_accesses
                add 3 to time_accesses of its parent
            for branch in endingAtIndex(-1):
                add 3 to time_accesses
                drop 3 from time_accesses of its parent
            update hash map accordingly
        */

    }

    Iterator begin(){
        return Iterator(&(this->v[0]));
    }

    Iterator end(){
        return Iterator(&(this->v[0]) + this->v.size());
    }
};