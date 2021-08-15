
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <bitset>
#include <cmath>
#include <cassert>
#include <cstring>

#if defined(__x86_64__)
  #include <immintrin.h>
#endif

#include "zstd.h"
#include <zlib.h>

#include "genotypes.h"
#include "utils.h"

namespace bgen {
  
// create lookup table of all probs for 8 bit integers. This has 511 entries, in
// order to allow for looking up values for the minor allele dosage, which can
// be up to 2.0 if a sample contains both copies are of the minor allele
// (2.0=255 + 255).
const float lut8[511] = {0.0000000, 0.0039216, 0.0078431, 0.0117647, 0.0156863, 0.0196078,
  0.0235294, 0.0274510, 0.0313725, 0.0352941, 0.0392157, 0.0431373, 0.0470588,
  0.0509804, 0.0549020, 0.0588235, 0.0627451, 0.0666667, 0.0705882, 0.0745098,
  0.0784314, 0.0823529, 0.0862745, 0.0901961, 0.0941176, 0.0980392, 0.1019608,
  0.1058824, 0.1098039, 0.1137255, 0.1176471, 0.1215686, 0.1254902, 0.1294118,
  0.1333333, 0.1372549, 0.1411765, 0.1450980, 0.1490196, 0.1529412, 0.1568627,
  0.1607843, 0.1647059, 0.1686275, 0.1725490, 0.1764706, 0.1803922, 0.1843137,
  0.1882353, 0.1921569, 0.1960784, 0.2000000, 0.2039216, 0.2078431, 0.2117647,
  0.2156863, 0.2196078, 0.2235294, 0.2274510, 0.2313725, 0.2352941, 0.2392157,
  0.2431373, 0.2470588, 0.2509804, 0.2549020, 0.2588235, 0.2627451, 0.2666667,
  0.2705882, 0.2745098, 0.2784314, 0.2823529, 0.2862745, 0.2901961, 0.2941176,
  0.2980392, 0.3019608, 0.3058824, 0.3098039, 0.3137255, 0.3176471, 0.3215686,
  0.3254902, 0.3294118, 0.3333333, 0.3372549, 0.3411765, 0.3450980, 0.3490196,
  0.3529412, 0.3568627, 0.3607843, 0.3647059, 0.3686275, 0.3725490, 0.3764706,
  0.3803922, 0.3843137, 0.3882353, 0.3921569, 0.3960784, 0.4000000, 0.4039216,
  0.4078431, 0.4117647, 0.4156863, 0.4196078, 0.4235294, 0.4274510, 0.4313725,
  0.4352941, 0.4392157, 0.4431373, 0.4470588, 0.4509804, 0.4549020, 0.4588235,
  0.4627451, 0.4666667, 0.4705882, 0.4745098, 0.4784314, 0.4823529, 0.4862745,
  0.4901961, 0.4941176, 0.4980392, 0.5019608, 0.5058824, 0.5098039, 0.5137255,
  0.5176471, 0.5215686, 0.5254902, 0.5294118, 0.5333333, 0.5372549, 0.5411765,
  0.5450980, 0.5490196, 0.5529412, 0.5568627, 0.5607843, 0.5647059, 0.5686275,
  0.5725490, 0.5764706, 0.5803922, 0.5843137, 0.5882353, 0.5921569, 0.5960784,
  0.6000000, 0.6039216, 0.6078431, 0.6117647, 0.6156863, 0.6196078, 0.6235294,
  0.6274510, 0.6313725, 0.6352941, 0.6392157, 0.6431373, 0.6470588, 0.6509804,
  0.6549020, 0.6588235, 0.6627451, 0.6666667, 0.6705882, 0.6745098, 0.6784314,
  0.6823529, 0.6862745, 0.6901961, 0.6941176, 0.6980392, 0.7019608, 0.7058824,
  0.7098039, 0.7137255, 0.7176471, 0.7215686, 0.7254902, 0.7294118, 0.7333333,
  0.7372549, 0.7411765, 0.7450980, 0.7490196, 0.7529412, 0.7568627, 0.7607843,
  0.7647059, 0.7686275, 0.7725490, 0.7764706, 0.7803922, 0.7843137, 0.7882353,
  0.7921569, 0.7960784, 0.8000000, 0.8039216, 0.8078431, 0.8117647, 0.8156863,
  0.8196078, 0.8235294, 0.8274510, 0.8313725, 0.8352941, 0.8392157, 0.8431373,
  0.8470588, 0.8509804, 0.8549020, 0.8588235, 0.8627451, 0.8666667, 0.8705882,
  0.8745098, 0.8784314, 0.8823529, 0.8862745, 0.8901961, 0.8941176, 0.8980392,
  0.9019608, 0.9058824, 0.9098039, 0.9137255, 0.9176471, 0.9215686, 0.9254902,
  0.9294118, 0.9333333, 0.9372549, 0.9411765, 0.9450980, 0.9490196, 0.9529412,
  0.9568627, 0.9607843, 0.9647059, 0.9686275, 0.9725490, 0.9764706, 0.9803922,
  0.9843137, 0.9882353, 0.9921569, 0.9960784, 1.0000000, 1.0039216, 1.0078431,
  1.0117647, 1.0156863, 1.0196078, 1.0235294, 1.0274510, 1.0313725, 1.0352941,
  1.0392157, 1.0431373, 1.0470588, 1.0509804, 1.0549020, 1.0588235, 1.0627451,
  1.0666667, 1.0705882, 1.0745098, 1.0784314, 1.0823529, 1.0862745, 1.0901961,
  1.0941176, 1.0980392, 1.1019608, 1.1058824, 1.1098039, 1.1137255, 1.1176471,
  1.1215686, 1.1254902, 1.1294118, 1.1333333, 1.1372549, 1.1411765, 1.1450980,
  1.1490196, 1.1529412, 1.1568627, 1.1607843, 1.1647059, 1.1686275, 1.1725490,
  1.1764706, 1.1803922, 1.1843137, 1.1882353, 1.1921569, 1.1960784, 1.2000000,
  1.2039216, 1.2078431, 1.2117647, 1.2156863, 1.2196078, 1.2235294, 1.2274510,
  1.2313725, 1.2352941, 1.2392157, 1.2431373, 1.2470588, 1.2509804, 1.2549020,
  1.2588235, 1.2627451, 1.2666667, 1.2705882, 1.2745098, 1.2784314, 1.2823529,
  1.2862745, 1.2901961, 1.2941176, 1.2980392, 1.3019608, 1.3058824, 1.3098039,
  1.3137255, 1.3176471, 1.3215686, 1.3254902, 1.3294118, 1.3333333, 1.3372549,
  1.3411765, 1.3450980, 1.3490196, 1.3529412, 1.3568627, 1.3607843, 1.3647059,
  1.3686275, 1.3725490, 1.3764706, 1.3803922, 1.3843137, 1.3882353, 1.3921569,
  1.3960784, 1.4000000, 1.4039216, 1.4078431, 1.4117647, 1.4156863, 1.4196078,
  1.4235294, 1.4274510, 1.4313725, 1.4352941, 1.4392157, 1.4431373, 1.4470588,
  1.4509804, 1.4549020, 1.4588235, 1.4627451, 1.4666667, 1.4705882, 1.4745098,
  1.4784314, 1.4823529, 1.4862745, 1.4901961, 1.4941176, 1.4980392, 1.5019608,
  1.5058824, 1.5098039, 1.5137255, 1.5176471, 1.5215686, 1.5254902, 1.5294118,
  1.5333333, 1.5372549, 1.5411765, 1.5450980, 1.5490196, 1.5529412, 1.5568627,
  1.5607843, 1.5647059, 1.5686275, 1.5725490, 1.5764706, 1.5803922, 1.5843137,
  1.5882353, 1.5921569, 1.5960784, 1.6000000, 1.6039216, 1.6078431, 1.6117647,
  1.6156863, 1.6196078, 1.6235294, 1.6274510, 1.6313725, 1.6352941, 1.6392157,
  1.6431373, 1.6470588, 1.6509804, 1.6549020, 1.6588235, 1.6627451, 1.6666667,
  1.6705882, 1.6745098, 1.6784314, 1.6823529, 1.6862745, 1.6901961, 1.6941176,
  1.6980392, 1.7019608, 1.7058824, 1.7098039, 1.7137255, 1.7176471, 1.7215686,
  1.7254902, 1.7294118, 1.7333333, 1.7372549, 1.7411765, 1.7450980, 1.7490196,
  1.7529412, 1.7568627, 1.7607843, 1.7647059, 1.7686275, 1.7725490, 1.7764706,
  1.7803922, 1.7843137, 1.7882353, 1.7921569, 1.7960784, 1.8000000, 1.8039216,
  1.8078431, 1.8117647, 1.8156863, 1.8196078, 1.8235294, 1.8274510, 1.8313725,
  1.8352941, 1.8392157, 1.8431373, 1.8470588, 1.8509804, 1.8549020, 1.8588235,
  1.8627451, 1.8666667, 1.8705882, 1.8745098, 1.8784314, 1.8823529, 1.8862745,
  1.8901961, 1.8941176, 1.8980392, 1.9019608, 1.9058824, 1.9098039, 1.9137255,
  1.9176471, 1.9215686, 1.9254902, 1.9294118, 1.9333333, 1.9372549, 1.9411765,
  1.9450980, 1.9490196, 1.9529412, 1.9568627, 1.9607843, 1.9647059, 1.9686275,
  1.9725490, 1.9764706, 1.9803922, 1.9843137, 1.9882353, 1.9921569, 1.9960784,
  2.0000000};

// uncompress a char array with zlib
void zlib_uncompress(char * input, int compressed_len, char * decompressed, int decompressed_len) {
  z_stream infstream;
  infstream.zalloc = Z_NULL;
  infstream.zfree = Z_NULL;
  infstream.opaque = Z_NULL;
  
  infstream.avail_in = compressed_len; // size of input
  infstream.next_in = (Bytef *) input; // input char array
  infstream.avail_out = decompressed_len; // size of output
  infstream.next_out = (Bytef *) decompressed; // output char array
  
  inflateInit(&infstream);
  inflate(&infstream, Z_NO_FLUSH);
  inflateEnd(&infstream);
  
  if (decompressed_len != (int) infstream.total_out) {
    throw std::invalid_argument("zlib decompression gave data of wrong length");
  }
}

// uncompress a char array with zstd
void zstd_uncompress(char * input, int compressed_len, char * decompressed,  int decompressed_len) {
  std::size_t total_out = ZSTD_decompress(decompressed, decompressed_len, input, compressed_len);
  if (decompressed_len != (int) total_out) {
    throw std::invalid_argument("zstd decompression gave data of wrong length");
  }
}

/// figure out the maximum number of probabilities across the individuals
///
/// @param max_ploidy maximum ploidy value across all samples
/// @param n_alleles number of alleles for the variant
/// @param phased whether the genotypes are phased
/// @return integer for maximum number of possible probabilites per sample
uint get_max_probs(int & max_ploidy, int & n_alleles, bool & phased) {
  uint max_probs;
  if (phased) {
    max_probs = n_alleles;
  } else {
    max_probs = n_choose_k(max_ploidy + n_alleles - 1, n_alleles - 1);
  }
  return max_probs;
}

/// get ploidy state for all samples (and missingness for layout2).
///
/// ploidy state is stored in ploidy member. Missingness is stored as indices of
/// samples This takes 100 microseconds for 500k samples. Layout 1 variants lack
/// ploidy info, and are assigned the same ploidy for all samples.
///
/// @param uncompressed char array possibly containing ploidy information
/// @param idx uint position where the ploidy data begins
void Genotypes::parse_ploidy(char * uncompressed, uint & idx) {
  if (has_ploidy) {
    if (layout == 2) {
      idx += n_samples;
    }
    return;
  }
  
  has_ploidy = true;
  ploidy = new std::uint8_t[n_samples];
  if (layout == 1) {
    std::memset(ploidy, max_ploidy, n_samples);
    return;
  }
  
  // we want to avoid parsing the ploidy states if  every sample has the same
  // ploidy. If we have a constant ploidy, set all entries to the same value
  std::uint8_t mask = 63;
  std::uint64_t mask_8 = std::uint64_t(0x8080808080808080);
  if (constant_ploidy) {
    std::memset(ploidy, max_ploidy, n_samples);
    for (uint x=0; x < (n_samples - (n_samples % 8)); x += 8) {
      // Simultaneously check if any of the next 8 samples are missing by casting
      // the data for the next 8 samples to an int64, and masking out all but
      // the bits which indicate missingness. Only check individual samples if
      // any are missing. This is ~3X quicker than looping across samples one by
      // one, provided the proportion of missing samples is low.
      if (*reinterpret_cast<const std::uint64_t*>(&uncompressed[idx + x]) & mask_8) {
        for (uint y=x; y < (x + 8); y++) {
          if (uncompressed[idx + x] & 0x80) {
            missing.push_back(x);
          }
        }
      }
    }
    // We looped through in batches of 8, so check the remainder not in an 8-batch
    for (uint x=(n_samples - (n_samples % 8)); x < n_samples; x++) {
      if (uncompressed[idx + x] & 0x80) {
        missing.push_back(x);
      }
    }
  } else {
    for (uint x=0; x < n_samples; x++) {
      ploidy[x] = mask & uncompressed[idx + x];
      if (uncompressed[idx + x] & 0x80) {
        missing.push_back(x);
      }
    }
  }
  idx += n_samples;
}

/// parse probabilities for layout1.
///
/// Layout 1 genotype probabilities are simple to parse, they are just 16-bit ints
/// one for each of the three possible genotypes. Missingness is encoded for a
/// sample when all probabilties are zero.
///
/// @param uncompressed char array containing genotype probabilities
/// @param idx uint position where the genotype probabilties begin
/// @return 1D float array of genotype probabilties (each from 0.0-1.0).
float * Genotypes::parse_layout1(char * uncompressed, uint & idx) {
  probs = new float[max_probs * n_samples];
  
  float factor = 1.0 / 32768;
  for (uint offset=0; offset<n_samples * max_probs; offset+=max_probs) {
    probs[offset] = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx]) * factor;
    probs[offset + 1] = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx + 2]) * factor;
    probs[offset + 2] = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx + 4]) * factor;
    idx += 6;
    
    if ((probs[offset] == 0.0) & (probs[offset + 1] == 0.0) & (probs[offset + 2] == 0.0)) {
      probs[offset] = std::nan("1");
      probs[offset + 1] = std::nan("1");
      probs[offset + 2] = std::nan("1");
    }
  }
  probs_parsed = true;
  return probs;
}

/// parse the initial data that defines the ploidy and phased status.
///
/// Layout 1 doesn't store any information before the genotype probabilities,
/// so layout 1 just receives default values.
///
/// @param uncompressed char array containing genotype probabilities
/// @param idx uint position where the genotype probabilties begin
void Genotypes::parse_preamble(char * uncompressed, uint & idx) {
  if (layout == 1) {
    phased = false;
    min_ploidy = 2;
    max_ploidy = 2;
    bit_depth = 16;
  } else if (layout == 2) {
    std::uint32_t nn_samples = *reinterpret_cast<const std::uint32_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint32_t);
    std::uint16_t allele_check = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint16_t);
    if (nn_samples != (std::uint32_t) n_samples) {
      throw std::invalid_argument("number of samples doesn't match!");
    }
    if (allele_check != n_alleles) {
      throw std::invalid_argument("number of alleles doesn't match!");
    }
    
    min_ploidy = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint8_t);
    max_ploidy = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint8_t);
  }
  
  constant_ploidy = (min_ploidy == max_ploidy);
  parse_ploidy(uncompressed, idx);
  
  if (layout == 2) {
    phased = (bool) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint8_t);
    bit_depth = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
    if ((bit_depth < 1) | (bit_depth > 32)) {
      throw std::invalid_argument("probabilities bit depth out of bounds");
    }
    idx += sizeof(std::uint8_t);
  }
  max_probs = get_max_probs(max_ploidy, n_alleles, phased);
}

/// parse probabilities for layout2
///
/// This is a fairly complex function, due to how the probabilties are encoded.
/// We need to extract probabilites per sample, but the final probability is
/// inferred as the remainder after all other genotype probabilties are parsed.
/// The number of bits used for encoding probabilties can differ between variants,
/// so there's a fast method for parsing 8-bit encoded probabilties, since that
/// is a common use case.
///
/// @param uncompressed char array containing genotype probabilities
/// @param idx uint position where the genotype probabilties begin
/// @return 1D float array of genotype probabilties (each from 0.0-1.0).
float * Genotypes::parse_layout2(char * uncompressed, uint & idx) {
  uint nrows = 0;
  if (!phased) {
    nrows = n_samples;
  } else {
    // phased probabilities require as many rows per sample as the ploidy
    if (constant_ploidy) {
      nrows = n_samples * max_ploidy;
    } else {
      for (uint n=0; n<n_samples; n++) { nrows += ploidy[n]; }
    }
  }
  probs = new float[max_probs * nrows];
  
  // get genotype/allele probabilities
  uint n_probs;
  uint max_less_1 = max_probs - 1;
  float prob = 0;
  float remainder;
  
  // define variables for parsing depths not aligned with 8 bit char array
  float factor = 1.0 / ((float) (std::pow(2, (int) bit_depth)) - 1);
  std::uint64_t probs_mask = std::uint64_t(0xFFFFFFFFFFFFFFFF) >> (64 - bit_depth);
  uint bit_idx = 0;  // index position in bits
  
  if (constant_ploidy & (max_probs == 3) & (bit_depth == 8)) {
    // A fast path for one scenario: all samples have ploidy=2, with 8 bits per
    // probability. This optimises memory accesses, and avoids looping over a
    // sample. This is ~2.5X faster than the standard route to compute
    // probabilities, and diploid samples with 8 bits/prob is likely the most
    // common use case, so the speed-up justifies this special case.
    std::uint64_t idx2 = 0;
    std::uint8_t first;
    std::uint8_t second;
    for (uint offset=0; offset < nrows * 3; offset += 3) {
      first = *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx + idx2]);
      second = *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx + idx2 + 1]);
      probs[offset] = lut8[first];
      probs[offset + 1] = lut8[second];
      probs[offset + 2] = lut8[255 - first - second];
      idx2 += 2;
    }
  } else {
    for (uint offset=0; offset < (nrows * max_probs); offset += max_probs) {
      // calculate the number of probabilities per sample (depends on whether the
      // data is phased, the sample ploidy and the number of alleles)
      if (constant_ploidy) {
        n_probs = max_less_1;
      } else if (phased) {
        n_probs = n_alleles - 1;
      } else if ((ploidy[offset / max_probs] == 2) && (n_alleles == 2)) {
        n_probs = 2;
      } else {
        n_probs = n_choose_k(ploidy[offset / max_probs] + n_alleles - 1, n_alleles - 1) - 1;
      }
      remainder = 1.0;
      for (uint x=0; x<n_probs; x++) {
        prob = ((*reinterpret_cast<const std::uint64_t* >(&uncompressed[idx + bit_idx / 8]) >> bit_idx % 8) & probs_mask) * factor ;
        bit_idx += bit_depth;
        remainder -= prob;
        probs[offset + x] = prob;
      }
      probs[offset + n_probs] = remainder;
      for (uint x=(n_probs + 1); x<max_probs; x++) {
        probs[offset + x] = std::nan("1");
      }
    }
  }
  
  uint offset;
  // for samples with missing data, just set values to NA
  for (auto n: missing) {
    offset = max_probs * n;
    for (uint x=0; x<max_probs; x++) {
      probs[offset + x] = std::nan("1");
    }
  }
  probs_parsed = true;
  return probs;
}

/// Read genotype data for a variant from disk and decompress.
///
/// The decompressed data is stored in the 'uncompressed' member. Decompression
/// is handled internally by either zlib_decompress, or zstd_decompress,
/// depending on compression scheme.
void Genotypes::decompress() {
  if (is_decompressed) {
    // don't decompress if already available
    return;
  }
  
  handle->seekg(offset);  // about 1 microsecond
  
  bool decompressed_field = false;
  std::uint32_t decompressed_len = 0;
  if (compression != 0) {
    if (layout == 1) {
      decompressed_len = n_samples * 6;
    } else if (layout == 2) {
      decompressed_field = true;
      handle->read(reinterpret_cast<char*>(&decompressed_len), sizeof(std::uint32_t));
    }
  }
  
  std::uint32_t compressed_len = next_var_offset - offset - decompressed_field * 4;
  char compressed[compressed_len];
  uncompressed = new char[decompressed_len];
  handle->read(&compressed[0], compressed_len); // about 20 microseconds
  
  if (compression == 0) { //no compression
    uncompressed = compressed;
  } else if (compression == 1) { // zlib
    zlib_uncompress(compressed, (int) compressed_len, uncompressed, (int) decompressed_len);  // about 2 milliseconds
  } else if (compression == 2) { // zstd
    zstd_uncompress(compressed, (int) compressed_len, uncompressed, (int) decompressed_len);
  }
  is_decompressed = true;
}

/// parse genotype data for a single variant
///
/// @return 1D float array of genotype probabilties (each from 0.0-1.0).
float * Genotypes::probabilities() {
  // avoid recomputation if called repeatedly for same variant
  if ((max_probs > 0) & probs_parsed) {
    return probs;
  }
  decompress();
  uint idx = 0;
  parse_preamble(uncompressed, idx);
  
  if (layout == 1) {
    probs = parse_layout1(uncompressed, idx);
  } else if (layout == 2) {
    probs = parse_layout2(uncompressed, idx);  // about 3 milliseconds
  }
  return probs;
}

/// find which allele corresponds to the minor allele
///
/// Rather than checking every individual to see which is the minor allele, we
/// check subsets, in batches of 100. We obtain alleles for individuals in the
/// batch, then check if a confidence interval for the frequency of the less
/// frequent allele could overlap 0.5. If not, we can be reasonably certain the
/// less frequent allele is the true minor allele, without having to check the
/// full cohort. This can be 600X faster than checking the full cohort in larger
/// populations.
///
/// @param dose float array of dosages for the reference (first) allele
/// @return index for minor allele (0 or 1)
int Genotypes::find_minor_allele(float * dose) {
  uint batchsize = 100;
  uint increment = std::max(n_samples / batchsize, (uint) 1);
  double total = 0;
  double freq = 0;
  
  // To make sure we don't hit weird groupings of alleles in individuals, this
  // picks samples uniformly thoughout the population, by using an appropriate
  // step size.
  for (uint idx2=0; idx2<increment; idx2++) {
    for (uint n=idx2; n<n_samples; n += increment) {
      total += dose[n];
    }
    freq = total / (batchsize * (idx2 + 1) * 2);
    if (minor_certain(freq, batchsize * (idx2 + 1), 5.0)) {
      break;
    }
  }
  
  if (freq <= 0.5) {
    return 0;
  } else {
    return 1;
  }
}

/// calculate dosage of the reference (first) allele for all samples.
///
/// Writes dosage float values to the dose member.
///
/// This is optimised for 8-bit genotypes and constant_ploidy. This function
/// pulls the homozygous minor allele genotype from the first genotype
/// probability. See fast_dosage_minor_second() for getting dosage when the minor
/// allele is the second allele.
///
/// This uses AVX2 vectorization to speed up calculations on x86_64 hardware.
///
/// @param uncompressed char array containing genotype probabilities
/// @param idx uint position where the genotype probabilties begin
void Genotypes::ref_dosage_fast(char * uncompressed, uint & idx) {
#if defined(__x86_64__)
  __m256i mask_odd = _mm256_set_epi8(0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0,
    -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0,  -1, 0, -1, 0, -1, 0, -1);
  __m256i mask_even = _mm256_set_epi8(-1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0,
    -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0);
  const float c = 1.0f / 255.0f;
  __m256 k = _mm256_set_ps(c, c, c, c, c, c, c, c);
  
  __m256i initial;
  __m256i first;
  __m256i second;
  __m128i lo16;
  __m256i lo;
  __m256 lo_float;
  __m128i hi16;
  __m256i hi;
  __m256 hi_float;
  for (uint n=0; n<(n_samples - (n_samples % 16)); n+=16) {
    initial = _mm256_loadu_si256((__m256i *) &uncompressed[idx]);
    
    // get heterozygous int dosage by masking out the even bytes, and right
    // shifting by one byte, to align homozygous and heterozygous counts. The
    // shift operation is from https://stackoverflow.com/a/25264853.
    second = _mm256_and_si256(initial, mask_even);
    second = _mm256_alignr_epi8(_mm256_permute2x128_si256(second, second, _MM_SHUFFLE(2, 0, 0, 1)), second, 1);
    
    // get homozygous int dosage by masking odd bytes, and left shift one bit to
    // multiply homozygous counts by 2 (since ploidy=2)
    first = _mm256_and_si256(initial, mask_odd);
    first = _mm256_slli_epi32(first, 1);
    
    // Now we have two 256 bit vectors with dosage counts adjusted for ploidy.
    // One for homozygous counts and one for heterozygous counts. Counts started
    // as 8-bit uints, interleaved with empty bytes (hom count spreads into
    // adjacent byte, as 9-bit uint). We can treat them as 16-bit ints for addition
    initial = _mm256_add_epi16(first, second);
    
    // convert the first half to floats, via 32 bit ints
    lo16 = _mm256_castsi256_si128(initial);
    lo = _mm256_cvtepi16_epi32(lo16);
    lo_float = _mm256_cvtepi32_ps(lo);
    
    // convert the second half to floats, via 32 bit ints
    hi16 = _mm256_extractf128_si256(initial, 1);
    hi = _mm256_cvtepi16_epi32(hi16);
    hi_float = _mm256_cvtepi32_ps(hi);
    
    _mm256_storeu_ps(&dose[n], _mm256_mul_ps(lo_float, k));
    _mm256_storeu_ps(&dose[n + 8], _mm256_mul_ps(hi_float, k));
    
    idx += 32;
  }
  // finish off the final unvectorized samples
  for (uint n=(n_samples - (n_samples % 16)); n<n_samples; n++) {
    dose[n] = lut8[*reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]) * 2 +
    *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx + 1])];
    idx += 2;
  }
#else
  // the unvectorized slow path. This is ~50% slower
  for (uint n=0; n<(n_samples - (n_samples % 2)); n+=2) {
    // speed up throughpot by calculating two samples at a time
    dose[n] = lut8[*reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]) * 2 +
      *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx + 1])];
    dose[n+1] = lut8[*reinterpret_cast<const std::uint8_t*>(&uncompressed[idx + 2]) * 2 +
      *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx + 3])];
    idx += 4;
  }
  // and finish off the final sample/s
  if (n_samples % 2) {
    dose[n_samples - 1] = lut8[*reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]) * 2 +
      *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx + 1])];
  }
#endif
}

/// calculate the dosage for the alternate (second) allele for all samples
///
/// This replaces the values in the dose member with 2.0 - value.
///
/// This uses AVX vectorization to speed up calculations on x86_64 hardware.
/// This assumes AVX is available on x86_64 machines, which isn't always true.
void Genotypes::alt_dosage() {
#if defined(__x86_64__)
  __m256 k = {2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f};
  __m256 batch;
  for (uint n=0; n<(n_samples - (n_samples % 8)); n+=8) {
    batch = _mm256_loadu_ps(dose + n);
    _mm256_storeu_ps(dose + n, _mm256_sub_ps(k, batch));
  }
  for (uint n=(n_samples - (n_samples % 8)); n<n_samples; n++) {
    dose[n] = 2.0f - dose[n];
  }
#else
  // TODO: add in vectorized version to speed up function on aarch64
  // alternative for when x86 vectorized oeprations are not available
  for (uint n=0; n<(n_samples - (n_samples % 8)); n+=8) {
    dose[n] = 2.0f - dose[n];
    dose[n+1] = 2.0f - dose[n+1];
    dose[n+2] = 2.0f - dose[n+2];
    dose[n+3] = 2.0f - dose[n+3];
    dose[n+4] = 2.0f - dose[n+4];
    dose[n+5] = 2.0f - dose[n+5];
    dose[n+6] = 2.0f - dose[n+6];
    dose[n+7] = 2.0f - dose[n+7];
  }
  for (uint n=(n_samples - (n_samples % 8)); n<n_samples; n++) {
    dose[n] = 2.0f - dose[n];
  }
#endif
}

/// calculate dosage of the reference (first) allele for all samples
///
/// The slow path, loops across samples. Figures out ploidy at each step,
/// and the homozygous genotype for the minor allele. Can use any bit depth.
/// For the most part we just have to use the first and second probabilities for
/// a sample, since they correspond to first hom, and het genotypes. This assumes
/// the first allele is for the ref genotype (but this isn't exposed externally,
/// since the calling function always converts to the minor allele dosage).
///
/// Variants encoded with layout 1 need a little bit more work, since we need to
/// skip over the second homozygous probability, and figure out missingness from
/// the probabilities.
///
/// @param uncompressed char array of genotype probabilities (encoding depends on layout)
/// @param idx uint index position in uncompressed where genotype probabilities start
void Genotypes::ref_dosage_slow(char * uncompressed, uint & idx) {
  uint ploidy = max_ploidy;
  uint half_ploidy = ploidy / 2;
  
  std::uint32_t maxval = std::pow(2, (std::uint32_t) (bit_depth)) - 1;
  float factor = (layout == 2) ? 1.0f / (float) maxval : 1.0f / 32768;
  std::uint32_t het;
  std::uint32_t hom;
  std::uint32_t hom_alt;
  std::uint64_t probs_mask = std::uint64_t(0xFFFFFFFFFFFFFFFF) >> (64 - bit_depth);
  uint bit_idx = 0;  // index position in bits
  for (uint n=0; n<n_samples; n++) {
    if (!constant_ploidy) {
      ploidy = this->ploidy[n];
      half_ploidy = ploidy / 2;
    }
    hom = ((*reinterpret_cast<const std::uint64_t* >(&uncompressed[idx + bit_idx / 8]) >> bit_idx % 8) & probs_mask);
    bit_idx += bit_depth;
    het = ((*reinterpret_cast<const std::uint64_t* >(&uncompressed[idx + bit_idx / 8]) >> bit_idx % 8) & probs_mask);
    bit_idx += bit_depth;
    dose[n] = ((hom * ploidy) + het * half_ploidy) * factor;
    if (layout == 1) {
      // layout1 stores hom alt probability, and all zeros indicates missingness
      hom_alt = *reinterpret_cast<const std::uint16_t* >(&uncompressed[idx + bit_idx / 8]);
      bit_idx += bit_depth;
      if ((hom == 0) & (het == 0) & (hom_alt == 0)) {
        missing.push_back(n);
      }
    }
  }
}

/// calculate minor allele dosage from the genotype probabilities
///
/// @return float array of dosage values (each from 0.0-2.0)
float * Genotypes::minor_allele_dosage() {
  if ((max_probs > 0) & dosage_parsed) {
    return dose;
  }
  decompress();
  uint idx = 0;
  parse_preamble(uncompressed, idx);
  
  if (n_alleles != 2) {
    throw std::invalid_argument("can't get allele dosages for non-biallelic var.");
  }
  
  // calculate the dosage for the first allele for all samples
  dose = new float[n_samples];
  if (constant_ploidy & (max_probs == 3) & (bit_depth == 8)) {
    // A fast path when we know the ploidy is constant and the bit depth is 8,
    // this avoids the bit shifts/masks used in the variable bit_depth path.
    ref_dosage_fast(uncompressed, idx);
  } else {
    ref_dosage_slow(uncompressed, idx);
  }
  
  minor_idx = find_minor_allele(dose);
  if (minor_idx != 0) {
    alt_dosage();
  }
  
  // for samples with missing data, just set values to NA
  for (auto n: missing) {
    dose[n] = std::nan("1");
  }
  
  dosage_parsed = true;
  return dose;
}

void Genotypes::clear_probs() {
  if (max_probs > 0) {
    delete[] ploidy;
    delete[] uncompressed;
  }
  if (probs_parsed) {
    delete[] probs;
  }
  if (dosage_parsed) {
    delete[] dose;
  }
  max_probs = 0;
}

} //namespace bgen
