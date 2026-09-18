/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/EvalAdmix.cpp
 * @author      Anders Albrechtsen
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef PCAONE_EVALADMIX_
#define PCAONE_EVALADMIX_

#include "Cmd.hpp"
#include "Data.hpp"

// correlation of residuals (van Waaij et al. 2023 projection estimator) from PCs
void run_evaladmix(Data* data, const Param& params);

#endif  // PCAONE_EVALADMIX_
