/*
 * utils.h
 *
 *  Created on: Dec 10, 2015
 *      Author: Kerstin Vater
 * Description: This file defines some useful functions when working with
 *              cellular lattice gas automatons.
 *
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "lgca_common.h"

#include "tclap/CmdLine.h"

// Gets values from the command line.
static inline void
get_vals_from_cmd(int argc, char **argv,
                  string* test_case,
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
    ValueArg<string> caseArg("t", "testcase", "Test case.", false, "pipe",
            "string (\"pipe\", \"box\", \"karman\", \"periodic\", \"collision\", \"diffusion\") (default: \"pipe\")");
    cmd.add(caseArg);

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
    *test_case              = caseArg.getValue();
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

inline Real interpolate(Real val, Real y0, Real x0, Real y1, Real x1) {

    return ((val - x0) * (y1 - y0) / (x1 - x0) + y0);
}

inline std::vector<Real> rgb(Real min, Real max, Real val) {

    Real range = max - min;

    assert(range > 1.0e-06);

    std::vector<Real> rgb_code(3, 0.0);

    Real red   = 0.0;
    Real green = 0.0;
    Real blue  = 0.0;

    if (val <= min) {

        red   = 0.0;
        green = 0.0;
        blue  = 1.0;

    } else if (val < min + 0.25 * range) {

        red   = 0.0;
        green = interpolate(val, 0.0, min, 1.0, min + 0.25 * range);
        blue  = 1.0;

    } else if (val < min + 0.5 * range) {

        red   = 0.0;
        green = 1.0;
        blue  = interpolate(val, 1.0, min + 0.25 * range, 0.0, min + 0.5 * range);

    } else if (val < min + 0.75 * range) {

        red   = interpolate(val, 0.0, min + 0.5 * range, 1.0, min + 0.75 * range);
        green = 1.0;
        blue  = 0.0;

    } else if (val < max) {

        red   = 1.0;
        green = interpolate(val, 1.0, min + 0.75 * range, 0.0, max);
        blue  = 0.0;

    } else if (val >= max ) {

        red   = 1.0;
        green = 0.0;
        blue  = 0.0;
    }

    rgb_code[0] = red;
    rgb_code[1] = green;
    rgb_code[2] = blue;

    return rgb_code;
}

#endif /* UTILS_H_ */
