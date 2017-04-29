/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef WEDGE_CUT_GENERATOR_HPP_
#define WEDGE_CUT_GENERATOR_HPP_

#include "params.hpp"
#include "qxx/dmat.hpp"
#include "geometry.hpp"
#include "knapsack2.hpp"
#include "single_row_cut_generator.hpp"

/**
 * Models a two-dimensional intersection cut, and its associated convex set.
 * The convex set is given by
 * \f[
 *     C = \{ x \in R^n : (d^i)^T(x-f) \leq 1 \},
 * \f]
 * where f is a point in the interior of C.
 */
class IntersectionCut {
public:
	Point f;
	Point *d;
	const int n_faces;

	/**
	 * Constructs a new intersection cut and its associated convex set from
	 * the data provided.
	 *
	 * @param f A fractional point in the interior of the convex set.
	 * @param n_faces The number of faces of the convex set.
	 * @param lines (Optional) Array of lines that support each face of the
	 *   convex.
	 */
	IntersectionCut(Point f, int n_faces, const Line *lines = 0);
	~IntersectionCut();

	/**
	 * Sets the desired face of the convex set.
	 *
	 * @param index Index of the face.
	 * @param line Line that supports the face.
	 */
	void set_face(int index, Line line);

	rational get_continuous_coefficient(rational rx, rational ry);

	double get_trivial_lifting_coefficient_double(double rx, double ry);
	rational get_trivial_lifting_coefficient(rational rx, rational ry);
	void pre_lifting();

private:
	double d_r0x, d_p, d_d0x, d_d0y, d_d1x, d_d1y;
	rational p, r0x;
	bool pre_lifting_ready;
};

/**
 * Models an intersection cut associated with a wedge.
 */
class WedgeCut : public IntersectionCut {
public:
	/**
	 * Constructs a wedge cut from the provided data.
	 *
	 * @param f A fractional point in the interior of the wedge.
	 * @param left A point on the left face of the wedge.
	 * @param apex The apex of the wedge.
	 * @param right A point on the right face of the wedge.
	 */
	WedgeCut(Point f, Point left, Point apex, Point right);
};

/**
 * Models an intersection cut associated with a split.
 */
class SplitCut : public IntersectionCut {
public:
	/**
	 * Constructs a split cut from the provided data.
	 *
	 * @param f A fractional point in the interior of the split
	 * @param left A point on the left face of the split
	 * @param right A poit on the right face of the split
	 * @param direction The direction in which the split in unbounded.
	 */
	SplitCut(Point f, Point left, Point right, Point direction);
};

/**
 * This class can be used to generate wedge cuts.
 */
class WedgeCutGenerator: public SingleRowCutGenerator {
private:
	bool finished;
	q::dvec f, r1;
	int r1_offset;
	int cur_facet;
	Knapsack2 knapsack;
	static int max_depth;
	int n_knapsacks;

private:
	Constraint *cut;
	void eval_next();
	static q::dvec intersection(
		q::dvec a, q::dvec b, q::dvec c, q::dvec d);
	
public:
	WedgeCutGenerator(Row &row);
	~WedgeCutGenerator();

	bool has_next();
	Constraint* next();
};



#endif /* WEDGE_CUT_GENERATOR_HPP_ */
