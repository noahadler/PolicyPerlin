#include <iostream>
#include <boost/random.hpp>
#include "perlin.hpp"

int main() {
  boost::random::mt19937 mersenneTwister;
  adler::Perlin<boost::random::mt19937> perlin(mersenneTwister);

  // One dimensional
  std::cout << perlin.noise<1>(0.1) << std::endl;

  // Three dimensional
  std::cout << perlin.noise<3>(0.1,0.2,0.3) << std::endl;

  // FIve dimensional
  float v[5] = {0.1,0.2,0.3,0.4,0.5};
  std::cout << perlin.noise<5>(v) << std::endl;
}

