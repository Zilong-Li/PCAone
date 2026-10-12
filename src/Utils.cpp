/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Utils.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Utils.hpp"

#include <sys/utsname.h>

#include <cstring>
#include <fstream>
#include <stdexcept>

std::string get_machine() {
  struct utsname unameData;
  if (uname(&unameData) != 0) {
    perror("uname");
    exit(EXIT_FAILURE);
  }
  std::string machine{unameData.machine};
  std::string node{unameData.nodename};
  std::string release{unameData.release};
  std::string version{unameData.version};
  std::string sysname{unameData.sysname};
  return "Machine name: " + machine + "\nNode name: " + node + "\nOperating system release: " + release +
         "\nOperating system version: " + version + "\nOperating system name: " + sysname + "\n";
}

String1D split_string(const std::string& s, const std::string& separators) {
  String1D ret;
  bool is_seperator[256] = {false};
  for (auto& ch : separators) {
    is_seperator[(unsigned int)ch] = true;
  }
  int begin = 0;
  for (int i = 0; i <= (int)s.size(); i++) {
    if (is_seperator[(uint8_t)s[i]] || i == (int)s.size()) {
      ret.push_back(std::string(s.begin() + begin, s.begin() + i));
      begin = i + 1;
    }
  }
  return ret;
}

size_t count_lines(const std::string& fpath) {
  std::ifstream fin(fpath);
  if (!fin.is_open()) throw std::invalid_argument("can not open " + fpath);
  size_t count = 0;
  std::string line;
  while (getline(fin, line)) {
    count++;
  }
  return count;
}

int tgets(gzFile gz, char** buf, uint64* size) {
  int rlen = 0;
  char* tok = gzgets(gz, *buf + rlen, *size - rlen);  // return buf or NULL
  if (!tok) return rlen;
  int tmp = tok ? strlen(tok) : 0;
  if (tok[tmp - 1] != '\n') {
    // expand buf size if no end-of-line found
    rlen += tmp;
    *size *= 2;
    *buf = (char*)realloc(*buf, *size);
  }
  rlen += tmp;
  return rlen;
}

FILE* fopenOrDie(const char* filename, const char* instruction) {
  FILE* const inFile = fopen(filename, instruction);
  if (inFile) return inFile;
  /* error */
  perror(filename);
  exit(1);
}

size_t freadOrDie(void* buffer, size_t sizeToRead, FILE* file) {
  size_t const readSize = fread(buffer, 1, sizeToRead, file);
  if (readSize == sizeToRead) return readSize; /* good */
  if (feof(file)) return readSize;             /* good, reached end of file */
  /* error */
  perror("fread");
  exit(4);  // error fread
}

void fcloseOrDie(FILE* file) {
  if (!fclose(file)) {
    return;
  };
  /* error */
  perror("fclose error");
  exit(1);
}

bool isZstdCompressed(const char* filename) {
  FILE* file = fopen(filename, "rb");
  if (!file) return false;

  char magicNumber[4];
  if (fread(magicNumber, 1, 4, file) != 4) {
    fclose(file);
    return false;
  }

  bool isCompressed = (ZSTD_isFrame(magicNumber, 4) != 0);

  fclose(file);
  return isCompressed;
}
