#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
using namespace std;

int BinomialCoefficient(const int n, const int k) {
  std::vector<int> aSolutions(k);
  aSolutions[0] = n - k + 1;

  for (int i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}

int GenerateRandomUniform(int rangeStart, int rangeEnd){
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>    distr(rangeStart, rangeEnd);
    return distr(generator);
}

int GenerateRandomRealUniform(int rangeStart, int rangeEnd){
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_real_distribution<int>    distr(rangeStart, rangeEnd);
    return distr(generator);
}