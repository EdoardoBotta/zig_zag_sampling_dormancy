#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <random>
#include <unordered_map>
#include <algorithm>
#include "utils.hpp"
using namespace std;

struct lineageCounts {
    int nAwake, nDormient
};

typedef struct lineageCounts LineageCountResult;

class Lineage{
    public:
    bool active;
    bool leaf;
    int nMutations;
    int mergerTime; // TODO: Implement constructor with automatic mergerTime computation.
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

    int GetMergerTime(){
        return *std::max_element(this->time_accesses.begin(), this->time_accesses.end());
    }

    void DecreaseMergerTime(){
        std::vector<int>::iterator position = std::find(this->time_accesses.begin(), this->time_accesses.end(), this->mergerTime);
        if (position != this->time_accesses.end()){
            this->time_accesses.erase(position);
            this->mergerTime -= 1;
        }
    }

    void IncreaseMergerTime(){
        this->time_accesses.push_back(this->mergerTime+1);
        this->mergerTime += 1;
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
            maxValue = branch.GetMergerTime();
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
    double GetBranchLength(Lineage branch){
        double length;
        std::vector<int> time_accesses = branch.GetTimeIntervalsUsed();

        for(int i : time_accesses){
            length += this->times[i];
        }

        return length;
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
    int CountMutations(){
        int mutationCount = 0;
        for (Lineage branch : this->v){
            mutationCount += branch.GetNumMutations();
        }
        return mutationCount;
    }

    public:
    int GetMutationsShorterChildAt(int timeIdx){
        // should sum branch mutations if timeIdx == 0
        double minLength, branchLength;
        Lineage shortestBranch;
        std::vector<Lineage> mergingAtIdx;
        
        
        if(timeIdx > this->CountTimeIntervals()-1){
            return this->CountMutations();
        }
        
        mergingAtIdx = this->endTimes[timeIdx];
        for (Lineage branch : mergingAtIdx){
            branchLength = this->GetBranchLength(branch);
            if(branchLength < minLength){
                shortestBranch = branch;
                minLength = branchLength;
            }
        }
        return shortestBranch.GetNumMutations();
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
    LineageCountResult CountLineagesUsing(int ind){
        LineageCountResult result;
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
        
        result.nAwake = nAwake;
        result.nDormient = nDormient;

        return result;
    }


    public:
    float ComputeLikelihood(int c, int K){
        LineageCountResult countResult;
        int nAwake, nDormient, m;
        double sumExp = 0;
    
        for (int i=0; i < this->CountTimeIntervals(); i++){
            countResult = this->CountLineagesUsing(i);
            
            nAwake = countResult.nAwake;
            nDormient = countResult.nDormient;
            
            sumExp += (BinomialCoefficient(nAwake, 2) + c*nAwake + c*K*nDormient)*this->GetTime(i);
        }
        
        return exp(-sumExp)*pow(c, 2*this->CountDormancyPeriods())*pow(K, this->CountDormancyPeriods());
    }

    public:
    Tree ProgressTreeTimes(std::vector <double> time_velocities, 
                            const double time){
        
        int progressedTime;

        if (time_velocities.size() != this->times.size()){
            throw std::invalid_argument("sizes do not match");
        }                        

        std::vector<double> new_times = this->times;

        for(int i=0; i < time_velocities.size(); i++){
            progressedTime = this->times[i] + time_velocities[i]*time;
            new_times[i] = std::max(progressedTime, 0);
        }

        return Tree(this->v, new_times);
    }

    public:
    bool IsSwapValid(const int ind){
        if(ind <= 0 || ind >= this->CountTimeIntervals()-1){
            return false;
        }
        std::vector<Lineage>& mergingAtIndex = this->endTimes[ind];
        for (Lineage branch : mergingAtIndex){
            if (branch.leftChild->isActive() && branch.leftChild->GetMergerTime() == ind - 1){
                return false;
            }
            if (branch.rightChild->isActive() && branch.rightChild->GetMergerTime() == ind - 1){
                return false;
            }
        }
        return true;
    }

    public:
    void SwapDown(const int ind){
        if(ind == 0 || ind == this->CountTimeIntervals()){
            throw std::invalid_argument("Invalid index for pushdown");
        }
        std::vector<Lineage>& mergingAtIndex = this->endTimes[ind];
        std::vector<Lineage>& mergingAtPreviousIndex = this->endTimes[ind-1];
        
        for (Lineage branch: mergingAtIndex){
            branch.DecreaseMergerTime();
        }
        Lineage parent = mergingAtIndex[0].parent;
        parent.IncreaseMergerTime();

        for (Lineage branch: mergingAtPreviousIndex){
            branch.IncreaseMergerTime();
        }
        parent = mergingAtIndex[0].parent
        parent.DecreaseMergerTime();

        this->endTimes[ind] = mergingAtPreviousIndex;
        this->endTimes[ind-1] = mergingAtIndex;
        /* 
        find branches merging at ind-1.
        If any of these is child of mergingAtIndex:
            for branch in mergingAtIndex:
                drop 3 from time_accesses
                add 3 to time_accesses of its parent
            for branch in mergingAtIndex(-1):
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