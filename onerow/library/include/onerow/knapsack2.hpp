#ifndef KNAPSACK_2_HPP_
#define KNAPSACK_2_HPP_
#include <vector>
#include <utility>
#include "qxx/rational.hpp"
#include "geometry.hpp"
using std::vector;
using std::pair;

// **************************************************************************
// 
// **************************************************************************
#define KNAPSACK2_LEFT	0
#define KNAPSACK2_RIGHT	1
#define KNAPSACK2_BOTH	2
#define KNAPSACK2_RAY	3

class PairRational {
public:
	size_t operator()(const pair<rational, rational> &k) const;
};

struct Knapsack2Vertex {

	int side;
	q::dvec lower, upper, opposed;

};


// **************************************************************************
// 
// **************************************************************************
class Knapsack2 {

public:
	Knapsack2();
	Knapsack2(rational f, rational r1);
	~Knapsack2();

	void clear();
	void eval(rational f, rational r1);
	
private:
	void push(int side,
		const q::dvec &l, const q::dvec &u, const q::dvec &o);

public:

	vector<Knapsack2Vertex> list;

};

#endif
