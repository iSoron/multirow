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

#include <stdexcept>
#include <gtest/gtest.h>
#include <onerow/wedge_cut_generator.hpp>

using std::set;
using std::cout;
using std::endl;

TEST(IntersectionCutTest, set_face_test)
{
	IntersectionCut ic(Point(rational(1,2),rational(1,4)), 1);

	// should calculate d correctly
	ic.set_face(0, Line(rational(127,64), rational(0), rational(67,64), rational(1)));
	EXPECT_EQ(ic.d[0].x, rational(4,5));
	EXPECT_EQ(ic.d[0].y, rational(3,4));

	// shoud work with any sequence of points
	ic.set_face(0, Line(rational(67,64), rational(1), rational(127,64), rational(0)));
	EXPECT_EQ(ic.d[0].x, rational(4,5));
	EXPECT_EQ(ic.d[0].y, rational(3,4));


	// should throw exception on badly defined lines
	EXPECT_THROW(ic.set_face(0, Line(rational(1), rational(1), rational(1), rational(1))), std::invalid_argument);

	// should throw exception on lines crossing f
	EXPECT_THROW(ic.set_face(0, Line(rational(1,2), rational(1,4), rational(1), rational(1))), std::invalid_argument);

	// should not throw exception on good lines
	EXPECT_NO_THROW(ic.set_face(0, Line(rational(0), rational(0), rational(1), rational(1))));

	// should throw exception on index out of range
	EXPECT_THROW(ic.set_face(-1, Line()), std::out_of_range);
	EXPECT_THROW(ic.set_face( 1, Line()), std::out_of_range);
	EXPECT_THROW(ic.set_face(50, Line()), std::out_of_range);
}

TEST(IntersectionCutTest, set_face_test_2)
{
	IntersectionCut ic(Point(rational(1,2),0), 1);

	ic.set_face(0, Line(1, 0, 0, 1));
	EXPECT_EQ(ic.d[0].x, rational(2));
	EXPECT_EQ(ic.d[0].y, rational(2));

	ic.set_face(0, Line(0, 1, 1, 0));
	EXPECT_EQ(ic.d[0].x, rational(2));
	EXPECT_EQ(ic.d[0].y, rational(2));
}

TEST(IntersectionCutTest, get_continuous_coefficient_test)
{
	IntersectionCut ic(Point(rational(5,7),rational(0)), 2);
	ic.set_face(0, Line(3, 4, rational(29,7), rational(40,7)));
	ic.set_face(1, Line(rational(29,7), rational(40,7), rational(2), rational(2)));

	// should calculate d correctly
	EXPECT_EQ(rational(-21, 8), ic.d[0].x);
	EXPECT_EQ(rational(  7, 4), ic.d[0].y);
	EXPECT_EQ(rational( 91,12), ic.d[1].x);
	EXPECT_EQ(rational(-35, 8), ic.d[1].y);

	// should calculate continuous coefficients correctly
	EXPECT_EQ(rational(91,12), ic.get_continuous_coefficient(rational(1),rational(0)));
	EXPECT_EQ(rational( 7, 4), ic.get_continuous_coefficient(rational(0),rational(1)));
	EXPECT_EQ(rational( 7,40), ic.get_continuous_coefficient(rational(3,5),1));
	EXPECT_EQ(rational(77,72), ic.get_continuous_coefficient(rational(1,3),rational(1,3)));
}

TEST(IntersectionCutTest, get_trivial_lifting_coefficient_test)
{
	IntersectionCut ic(Point(rational(5,7),rational(0)), 2);
	ic.set_face(0, Line(3, 4, rational(29,7), rational(40,7)));
	ic.set_face(1, Line(rational(29,7), rational(40,7), rational(2), rational(2)));

	// should calculate d correctly
	EXPECT_EQ(rational(-21, 8), ic.d[0].x);
	EXPECT_EQ(rational(  7, 4), ic.d[0].y);
	EXPECT_EQ(rational( 91,12), ic.d[1].x);
	EXPECT_EQ(rational(-35, 8), ic.d[1].y);

	// should calculate lifted coefficients correctly
	EXPECT_EQ(rational(399,704), ic.get_trivial_lifting_coefficient(rational(119,264),rational(0)));
}

TEST(IntersectionCutTest, get_trivial_lifting_coefficient_test_2)
{
	IntersectionCut ic(Point(rational(1,2),rational(0)), 2);
	ic.set_face(0, Line(0,0,0,1));
	ic.set_face(1, Line(1,0,0,1));

	// should calculate d correctly
	EXPECT_EQ(rational(-2), ic.d[0].x);
	EXPECT_EQ( rational(0), ic.d[0].y);
	EXPECT_EQ( rational(2), ic.d[1].x);
	EXPECT_EQ( rational(2), ic.d[1].y);

	// should calculate lifted coefficients correctly
	EXPECT_EQ(rational(1), ic.get_trivial_lifting_coefficient(rational(-1,2),0));
}

TEST(WedgeCutGenerator, generate_test_1)
{
	Row r;
	r.basic_var_index = 999;
	r.is_integer = new bool[3];
	r.c.pi.resize(3);
	r.c.pi.push(0, rational(-3,5));
	r.c.pi.push(1, rational(-1));
	r.c.pi.push(2, rational(1));
	r.c.pi_zero = rational(5,7);

	std::set<Constraint> expected_cuts;
	rational cut_matrix[] =  {
			rational(-14,25), rational(-7,2), rational(-7,5),
			rational(-7,20), rational(-7,2), rational(-91,44),
			rational(-7,40), rational(-91,12), rational(-21,8),
			rational(0), rational(-35,3), rational(-35,4),
//			rational(-21,10), rational(0), rational(-7,5),
//			rational(-21,10), rational(-7,2), rational(0),
	};

	for(int i=0; i<3; i++)
		r.is_integer[i] = (i == 0);

	for(int i=0; i<4; i++)
	{
		Constraint c;
		c.pi.resize(3);
		c.pi_zero = -1;
		for(int j=0; j<3; j++)
			c.pi.push(j, cut_matrix[i*3+j]);
		expected_cuts.insert(c);
	}

	std::set<Constraint> actual_cuts;

	WedgeCutGenerator wcg(r);
	while(wcg.has_next())
		actual_cuts.insert(*wcg.next());

//	  cout << "EXPECTED:" << endl;
//	  for(Constraint c : expected_cuts)
//		cout << c << endl;
//
//	  cout << "ACTUAL:" << endl;
//	  for(Constraint c : actual_cuts)
//			cout << c << endl;

	EXPECT_TRUE(expected_cuts == actual_cuts);
}

TEST(WedgeCutGenerator, generate_test_2)
{
	Row r;
	r.basic_var_index = 999;
	r.is_integer = new bool[3];
	r.c.pi.resize(3);
	r.c.pi.push(0, rational(1,2));
	r.c.pi.push(1, rational(-1));
	r.c.pi.push(2, rational(1));
	r.c.pi_zero = rational(1,2);

	set<Constraint> expected_cuts;
	rational cut_matrix[] =  {
			rational(-1), rational(-2), rational(-2)
	};

	for(int i=0; i<5; i++)
		r.is_integer[i] = (i == 0);

	for(int i=0; i<1; i++)
	{
		Constraint c;
		c.pi.resize(3);
		c.pi_zero = -1;
		for(int j=0; j<3; j++)
		{
			c.pi.push(j, cut_matrix[i*3+j]);
		}
		expected_cuts.insert(c);
	}

	set<Constraint> actual_cuts;
	WedgeCutGenerator wcg(r);
	while(wcg.has_next())
		actual_cuts.insert(*wcg.next());

//	  cout << "EXPECTED:" << endl;
//	  for(Constraint c : expected_cuts)
//		cout << c << endl;
//
//	  cout << "ACTUAL:" << endl;
//	  for(Constraint c : actual_cuts)
//			cout << c << endl;

	EXPECT_TRUE(expected_cuts == actual_cuts);
}

TEST(WedgeCutGenerator, generate_test_3)
{
	bool is_integer[] = { true, false, false, true, true, true, true, true, false,
		false, false, false, false, false };
	Row r;
	r.basic_var_index = 0;
	r.is_integer = is_integer;
	r.c.pi.resize(14);
	r.c.pi.push(0, rational(1));
	r.c.pi.push(1, rational(7338,415411));
	r.c.pi.push(2, rational(1271,71708));
	r.c.pi.push(3, rational(19,577624));
	r.c.pi.push(4, rational(-19,577624));
	r.c.pi.push(5, rational(7,189263));
	r.c.pi.push(6, rational(-7,189263));
	r.c.pi.push(7, rational(-5,83513));
	r.c.pi.push(8, rational(157727,47951));
	r.c.pi.push(9, rational(-157727,47951));
	r.c.pi.push(10, rational(536435,145039));
	r.c.pi.push(11, rational(-536435,145039));
	r.c.pi.push(12, rational(237478,39665));
	r.c.pi.push(13, rational(237478,39665));
	r.c.pi_zero = rational(631975,707651);

	rational solution[14];
	solution[0] = rational(0);
	solution[1] = rational(0);
	solution[2] = rational(0);
	solution[3] = rational(0);
	solution[4] = rational(0);
	solution[5] = rational(0);
	solution[6] = rational(0);
	solution[7] = rational(1);
	solution[8] = rational(5531144031239411,576460752303423488);
	solution[9] = rational(8350617490688559,1152921504606846976);
	solution[10] = rational(6818955564856927,576460752303423488);
	solution[11] = rational(6377917030130699,4611686018427387904);
	solution[12] = rational(2547758494699757,18014398509481984);
	solution[13] = rational(0);

	cout << r.c.get_violation(solution).get_double() << endl;
	EXPECT_TRUE(r.c.get_violation(solution) < 0.001);

	WedgeCutGenerator wcg(r);
	while(wcg.has_next())
	{
		Constraint *c = wcg.next();
		cout << *c << endl;
		cout << c->get_violation(solution).get_double() << endl;
		EXPECT_TRUE(c->get_violation(solution) < 0.001);
	}
}
