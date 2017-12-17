/*
 * This file is part of LGCA, an implementation of a Lattice Gas Cellular Automaton
 * (https://github.com/keva92/lgca).
 *
 * Copyright (c) 2015-2017 Kerstin Vater, Niklas Kühl, Christian F. Janßen.
 *
 * LGCA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * LGCA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with lgca. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LGCA_UTILS_H_
#define LGCA_UTILS_H_

#include "lgca_common.h"

#include "tclap/CmdLine.h"

namespace lgca {

// Gets values from the command line.
static inline void
get_vals_from_cmd(int argc, char **argv,
                  Real* Re, Real* Ma,
                  int* n_dir,
                  int* s_max,
                  int* coarse_graining_radius,
                  int* write_steps,
                  int* body_force_steps, int* body_force_intensity,
                  int* device,
                  int* max_block_size,
                  string* parallel_type,
                  string* output_format) {

    // Define the command line object.
    TCLAP::CmdLine cmd("Command description message", ' ', "0.9");

    // Define value arguments and add them to the command line.
    using TCLAP::ValueArg;

    ValueArg<Real>  ReArg("r", "Re", "Reynolds number.", false, 80.0, "real gt 0 (default: 80.0)");
    cmd.add(ReArg);

    ValueArg<Real> MaArg("m", "Ma", "Mach number.", false, 0.2, "real gt 0 (default: 0.2");
    cmd.add(MaArg);

    ValueArg<int> ndirArg("d", "n-dir", "Number of lattice directions.", false, 6, "int 4 (HPP) or int 6 (FHP) (default: FHP)");
    cmd.add(ndirArg);

    ValueArg<int> smaxArg("s", "steps", "Number of simulated time steps.", false, 50, "int gte 0 (default: 50)");
    cmd.add(smaxArg);

    ValueArg<int> cgArg("c", "cg-radius", "Coarse graining radius.", false, 15, "int gte 0 (default: 15)");
    cmd.add(cgArg);

    ValueArg<int> writeArg("w", "write-steps", "Number of steps after which the post-processed results are written to a file.", false, 10, "int gt 0 (default: 10)");
    cmd.add(writeArg);

    ValueArg<int> forceArg("", "bf-steps", "Number of steps after which a body force is applied to the particles.", false, 100, "int gt 0 (default: 100)");
    cmd.add(forceArg);

    ValueArg<int> bfIntArg("", "bf-int", "Intensity of the body force.", false, 250, "int gte 0 (default: 250)");
    cmd.add(bfIntArg);

    ValueArg<int> deviceArg("", "device", "Number of the device to use.", false, 0, "int gte 0 (default: 0)");
    cmd.add(deviceArg);

    ValueArg<int> blockSizeArg("", "blocksize", "Maximum block size in x direction.", false, 256, "int gt 0 (default: 256)");
    cmd.add(blockSizeArg);

    ValueArg<string> parallelArg("p", "parallel", "Parallelization type.", false, "OMP", "string (\"CUDA\" or \"OMP\") (default: \"OMP\")");
    cmd.add(parallelArg);

    ValueArg<string> outputArg("o", "output", "Output mode.", false, "live", "string (\"live\" or \"vti\") (default: \"live\")");
    cmd.add(outputArg);

    // Parse the args.
    cmd.parse(argc, argv);

    // Get the value parsed by each arg.
    *Re                     = ReArg.getValue();
    *Ma                     = MaArg.getValue();
    *n_dir                  = ndirArg.getValue();
    *s_max                  = smaxArg.getValue();
    *coarse_graining_radius = cgArg.getValue();
    *write_steps            = writeArg.getValue();
    *body_force_steps       = forceArg.getValue();
    *body_force_intensity   = bfIntArg.getValue();
    *device                 = deviceArg.getValue();
    *max_block_size         = blockSizeArg.getValue();
    *parallel_type          = parallelArg.getValue();
    *output_format          = outputArg.getValue();
}

// Prints a startup message.
static inline void print_startup_message() {

    printf("\n");
    printf("                  .::         .::::       .::         .:         \n");
    printf("                  .::       .:    .::  .::   .::     .: ::       \n");
    printf("                  .::      .::        .::           .:  .::      \n");
    printf("                  .::      .::        .::          .::   .::     \n");
    printf("                  .::      .::   .::::.::         .:::::: .::    \n");
    printf("                  .::       .::    .:  .::   .:: .::       .::   \n");
    printf("                  .::::::::  .:::::      .::::  .::         .::  \n");
    printf("\n");
    printf("This is a 2D Lattice Gas Cellular Automaton based on the HPP and FHP model.\nHave fun!\n\n");
}

// Returns a random real number within [0.0 and 1.0].
inline Real random_uniform() {

    return (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
}

// Generates random boolean values with 50:50 distribution.
inline bool random_bool() {

	return (random() > 0.5);
}

}

#endif /* LGCA_UTILS_H_ */
