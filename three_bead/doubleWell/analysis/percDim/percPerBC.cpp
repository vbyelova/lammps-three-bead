#include "percPerBC.h"


#pragma mark -
#pragma mark Debug spanning dimension
#pragma mark -

// 1D checks.
template<>
void percPerBC<1>::debugSpanningDimension( bool verbose )
{
	std::multimap< int, intVecSet > tests1D {
		{ 1, { {1},{2} } }
	};

	for( auto &test : tests1D ) __testSpanningDimension( test.second, test.first, verbose );
}

// 2D checks.
template<>
void percPerBC<2>::debugSpanningDimension( bool verbose )
{
	std::multimap< int, intVecSet > tests2D {
		{ 0, { {0,0},{0,0} } },
		{ 2, { {0,1},{1,0} } },
		{ 1, { {1,2},{2,4},{0,0} } },
		{ 2, { {1,2},{-1,0} } }
	};

	for( auto &test : tests2D ) __testSpanningDimension( test.second, test.first, verbose );
}

// 3D checks.
template<>
void percPerBC<3>::debugSpanningDimension( bool verbose )
{
	std::multimap< int, intVecSet > tests3D {
		{ 3, { {1,0,0},{0,1,0},{0,0,1} } },
		{ 2, { {1,2,3},{4,5,6},{1,0,-1} } },		// nb. 5 e_{1} - 2 e_{2} + 3 e_{3} = 0, so spanning dimension is 2.
		{ 1, { {1,2,0},{2,4,0},{-3,-6,0} } }
	};

	for( auto &test : tests3D ) __testSpanningDimension( test.second, test.first, verbose );
}


#pragma mark -
#pragma mark Debug independence
#pragma mark -

// 1D checks.
template<>
void percPerBC<1>::debugIsIndependent( bool verbose )
{
	__testIsIndependent( false, { {{1}} }, {{2}}, verbose );
	__testIsIndependent( true , { {{0}} }, {{-1}}, verbose );
}

// 2D checks.
template<>
void percPerBC<2>::debugIsIndependent( bool verbose )
{
	__testIsIndependent( true , { {0,1},{0,2} }, {1,0}, verbose );
	__testIsIndependent( false, { {0,1},{1,2} }, {1,0}, verbose );
}

// 3D checks.
template<>
void percPerBC<3>::debugIsIndependent( bool verbose )
{
	__testIsIndependent( true , { {1,2,3},{1,3,2} }, {1,0,0}, verbose );
	__testIsIndependent( false, { {1,2,3},{1,3,2} }, {1,2,3}, verbose );
}