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

#ifndef LGCA_OMP_LATTICE_H_
#define LGCA_OMP_LATTICE_H_

#include "lattice.h"

namespace lgca {

template<Model model_>
class OMP_Lattice: public Lattice<model_> {

private:

    using ModelDesc = ModelDescriptor<model_>;

    // Auxiliary array on the CPU
    unsigned char* m_node_state_tmp_cpu;

    // Model-based values according to the number of lattice directions
    ModelDesc* m_model;

    // Computes cell quantities of interest as a post-processing procedure
	void cell_post_process();

    // Computes coarse grained quantities of interest as a post-processing procedure
	void mean_post_process();

    // Allocates the memory for the arrays on the host (CPU) and device (GPU)
    void allocate_memory();

    // Frees the memory for the arrays on the host (CPU) and device (GPU)
    void free_memory();

public:

	// Creates an openMP parallelized lattice gas cellular automaton object
	// of the specified properties.
    OMP_Lattice(const string m_test_case,
                const Real m_Re, const Real m_Ma_s,
                const int m_coarse_graining_radius);

	virtual ~OMP_Lattice();

    // Sets (proper) parallelization parameters.
    void setup_parallel();

    // Performs the collision and propagation step on the lattice gas automaton.
    void collide_and_propagate(const bool p);

    // Computes the mean velocity of the lattice.
    std::vector<Real> get_mean_velocity();

    // Applies a body force in the specified direction (x or y) and with the
    // specified intensity to the particles. E.g., if the intensity is equal 100,
    // every 100th particle changes it's direction, if feasible.
    void apply_body_force(const int forcing);

    // Computes quantities of interest as a post-processing procedure.
    void post_process();
};

} // namespace lgca

#endif /* LGCA_OMP_LATTICE_H_ */
