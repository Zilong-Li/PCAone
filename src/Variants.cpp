/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Variants.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Variants.hpp"

#include <cstring>  // strtok_r
#include <fstream>

#include "Cmd.hpp"
#include "Utils.hpp"

namespace {

const std::string SEP{" \t"};

String1D pvar_tokens_to_bim_tokens(const String1D& tokens, const std::string& path, const std::string& line) {
  if ((int)tokens.size() < 5) cao.error("the input pvar file is not valid!\n => " + path + "\n" + line);
  std::string alt = tokens[4];
  size_t comma = alt.find(',');
  if (comma != std::string::npos) alt = alt.substr(0, comma);
  // PGEN reader uses ALT allele dosage (allele_idx=1), so expose ALT as BIM A1.
  return String1D{tokens[0], tokens[2], "0", tokens[1], alt, tokens[3]};
}

// CHR_BP_A1_A2 of a .bim line, and the same with the alleles swapped
std::string bim_match_key(const String1D& tokens, const std::string& path, const std::string& line) {
  if ((int)tokens.size() < 6) cao.error("the input bim file is not valid!\n => " + path + "\n" + line);
  return tokens[0] + "_" + tokens[3] + "_" + tokens[4] + "_" + tokens[5];
}
std::string bim_flip_key(const String1D& tokens, const std::string& path, const std::string& line) {
  if ((int)tokens.size() < 6) cao.error("the input bim file is not valid!\n => " + path + "\n" + line);
  return tokens[0] + "_" + tokens[3] + "_" + tokens[5] + "_" + tokens[4];
}

// BEAGLE may code the alleles 0..3
std::string decode_beagle_allele(const std::string& allele) {
  if (allele == "0") return "A";
  if (allele == "1") return "C";
  if (allele == "2") return "G";
  if (allele == "3") return "T";
  return allele;
}

// the keys of the variants of a .bim or .pvar, in its order
String1D bim_keys(const std::string& path, bool pvar) {
  std::ifstream fin(path);
  if (!fin.is_open()) cao.error("can not open " + path);
  String1D keys;
  std::string line;
  while (read_bim_line(fin, pvar, path, line)) keys.push_back(bim_match_key(split_string(line, SEP), path, line));
  return keys;
}

// The keys of the markers of a BEAGLE file, as marker_allele2_allele1. BEAGLE
// GLs count allele2 (AA, AB, BB), while the .mbim counts A1 (BED2GENO; PGEN
// exposes ALT as A1), so the counted alleles line up only with allele2 first.
// Keying allele1 first flipped exactly the sites that were already aligned.
String1D beagle_keys(const std::string& path) {
  gzFile fp = gzopen(path.c_str(), "r");
  if (!fp) cao.error("can not open " + path);
  uint64 bufsize = (uint64)128 * 1024 * 1024;
  char* original = (char*)calloc(bufsize, sizeof(char));
  char* buffer = original;
  const char* delims = "\t \n";
  tgets(fp, &buffer, &bufsize);  // skip header line
  if (buffer != original) original = buffer;
  buffer = original;
  String1D keys;
  while (tgets(fp, &buffer, &bufsize)) {
    if (buffer != original) original = buffer;
    char* marker = strtok_r(buffer, delims, &buffer);
    char* allele1 = strtok_r(NULL, delims, &buffer);
    char* allele2 = strtok_r(NULL, delims, &buffer);
    if (!marker || !allele1 || !allele2) cao.error("invalid BEAGLE record while matching markers:\n => " + path);
    keys.push_back(std::string(marker) + "_" + decode_beagle_allele(allele2) + "_" + decode_beagle_allele(allele1));
    buffer = original;
  }
  gzclose(fp);
  free(original);
  return keys;
}

// the input's keys in the .mbim (chr id cM bp a1 a2 maf), as they are or with
// the alleles swapped
BimMatch match_keys_to_mbim(const String1D& keys, const std::string& mbim_file) {
  std::ifstream fm(mbim_file);
  if (!fm.is_open()) cao.error("can not open " + mbim_file);
  String1D mbim_keys;
  std::unordered_map<std::string, int> mbim_lookup, mbim_flip_lookup;
  std::string line;
  while (getline(fm, line)) {
    auto tokens = split_string(line, SEP);
    if ((int)tokens.size() != 7) cao.error("the input file is not valid!\n => " + mbim_file);
    int idx = (int)mbim_keys.size();
    std::string key = bim_match_key(tokens, mbim_file, line);
    if (!mbim_lookup.insert({key, idx}).second)
      cao.error("duplicate SNP records found in " + mbim_file + " for key " + key);
    mbim_flip_lookup.insert({bim_flip_key(tokens, mbim_file, line), idx});
    mbim_keys.push_back(std::move(key));
  }

  BimMatch match;
  match.identical = keys == mbim_keys;
  const size_t nmax = std::min(keys.size(), mbim_keys.size());
  match.bim_indices.reserve(nmax);
  match.mbim_indices.reserve(nmax);
  match.flip.reserve(nmax);
  for (int i = 0; i < (int)keys.size(); ++i) {
    auto it = mbim_lookup.find(keys[i]);
    bool flip = false;
    if (it == mbim_lookup.end()) {
      it = mbim_flip_lookup.find(keys[i]);
      if (it == mbim_flip_lookup.end()) continue;
      flip = true;
    }
    match.bim_indices.push_back(i);
    match.mbim_indices.push_back(it->second);
    match.flip.push_back(flip);
  }
  return match;
}

}  // namespace

bool read_bim_line(std::istream& is, bool pvar, const std::string& path, std::string& line) {
  while (std::getline(is, line)) {
    if (!pvar) return true;
    if (line.empty() || line[0] == '#') continue;
    line = pvar_line_to_bim_line(line, path);
    return true;
  }
  return false;
}

std::string pvar_line_to_bim_line(const std::string& line, const std::string& path) {
  auto tokens = pvar_tokens_to_bim_tokens(split_string(line, SEP), path, line);
  return tokens[0] + "\t" + tokens[1] + "\t" + tokens[2] + "\t" + tokens[3] + "\t" + tokens[4] + "\t" + tokens[5];
}

BimMatch match_variants_to_mbim(const Param& params) {
  std::string path;
  String1D keys;
  if (params.file_t == FileType::BEAGLE) {
    path = params.filein;
    keys = beagle_keys(path);
  } else {
    const bool pvar = params.file_t == FileType::PGEN;
    path = params.filein + (pvar ? ".pvar" : ".bim");
    keys = bim_keys(path, pvar);
  }
  BimMatch match = match_keys_to_mbim(keys, params.filebim);
  if (match.bim_indices.empty()) cao.error("no overlapped SNPs found between " + path + " and " + params.filebim);
  return match;
}
