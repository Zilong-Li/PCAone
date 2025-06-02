// -*- compile-command: "g++ summarise_ld_r2bin.cpp -o summarise_ld_r2bin -O3 -std=c++17 -lz"; -*-
#include <zlib.h>
#include <string>
#include <functional>
#include <sstream>
#include <iostream>
#include <map>


static bool parseGzipLineByLine(const std::string& filename, 
                                std::function<void(const std::string&)> callback) {
  gzFile file = gzopen(filename.c_str(), "rb");
  if (!file) {
    return false;
  }
    
  char buffer[1024]; // should be enough
  std::string line;
    
  while (gzgets(file, buffer, sizeof(buffer))) {
    line = buffer;
        
    // Remove trailing newline if present
    if (!line.empty() && line.back() == '\n') {
      line.pop_back();
    }
        
    callback(line);
  }
    
  bool success = !gzeof(file) || gzerror(file, nullptr) == Z_OK;
  gzclose(file);
  return success;
}


int main (int argc, char** argv) {

  std::vector<std::string> args(argv + 1, argv + argc);
  if(argc <= 1 || args[0] == "-h" || args[0] == "-help" || args[0] == "--help"){
    std::cout << "-i input file pcaone.ld.gz\n";
    return 1;
  }

  std::string in;
  for(size_t i = 0; i < args.size(); i++){
    if(args[i] == "-i") in = args[++i];
  }


  std::map<int, std::pair<double, int>> iChr;
  // std::vector<std::map<int, std::pair<double, int>>> res;

  int dA{0}, dB{0}, d{0}, num{0};
  double r2{0};
  int n = 0; // the n-th line
  std::istringstream iss;
  std::string tok;
  
  parseGzipLineByLine(in, [&](const std::string& line) {
    n++;
    if(n==1) return;
    int j = 0;  // the j-th column
    iss.clear();
    iss.str(line);
    while (iss >> tok) {
      if(j == 6) r2 = std::stod(tok);
      if(j == 1) dA = std::stoi(tok);
      if(j == 4) dB = std::stoi(tok);
      j++;
    }
    d = dB - dA;
    iChr[d].first += r2;
    iChr[d].second += 1;
  });

  std::cout << "Distance\tR2Total\tnumPairs\n";
  for (const auto& [k, v] : iChr) {
    std::cout << k << "\t" << v.first << "\t" << v.second << "\n";
  }

  return 0;
}

