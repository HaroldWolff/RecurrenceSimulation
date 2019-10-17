#include "rndm.hpp"

Rndm::Rndm(unsigned int seed)
{
    srand(seed);
}

Rndm::~Rndm()
{

}

double Rndm::Uniform()
{

    long double tempUniform = double(rand())/10000.;
    int i = 0;

	while ( tempUniform < 0. || 1. < tempUniform )
    {
        tempUniform = double(rand())/10000.;
        i++;
    }

    return ( tempUniform );
}

double Rndm::Uniform(double under, double upper)
{

    long double tempUniform = double(rand())/10000.;
    int i = 0;
    double difference = upper - under;

	while ( tempUniform < 0. || 1. < tempUniform )
    {
        tempUniform = double(rand())/10000.;
        i++;
    }

    tempUniform = under + (tempUniform * difference);

    return ( tempUniform );
}

double Rndm::Normal( double mean, double sd )
/*
 * This is a modification of the Kinderman + Monahan algorithm for
 * generating normal random numbers, due to Leva:
 *
 * J.L. Leva, Algorithm 712. A normal random number generator, ACM Trans.
 * Math. Softw.  18 (1992) 454--455.
 *
 * http://www.acm.org/pubs/citations/journals/toms/1992-18-4/p449-leva/
 *
 * Note: Some of the constants used below look like they have dubious
 * precision.  These constants are used for an approximate bounding
 * region test (see the paper).  If the approximate test fails,
 * then an exact region test is performed.
 *
 * Only 0.012 logarithm evaluations are required per random number
 * generated, making this method comparatively fast.
 *
 * Adapted to C++ by T. Veldhuizen.
 * http://pyprop.googlecode.com/svn/trunk/extern/blitz/source/random/normal.h
 */
{
	const double s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472;
	const double r1 = 0.27597, r2 = 0.27846;

	double u, v;

	for (;;)
	{
	// Generate P = (u,v) uniform in rectangle enclosing
	// acceptance region:
	//   0 < u < 1
	// - sqrt(2/e) < v < sqrt(2/e)
	// The constant below is 2*sqrt(2/e).

	u = Uniform();
	v = 1.715527769921413592960379282557544956242L
		* (Uniform() - .5 );

	// Evaluate the quadratic form
	double x = u - s;
	double y = fabs(v) - t;
	double q = x*x + y*(a*y - b*x);

	// Accept P if inside inner ellipse
	if (q < r1)
		break;

	// Reject P if outside outer ellipse
	if (q > r2)
		continue;

	// Between ellipses: perform exact test
	if (v*v <= -4.0 * log(u)*u*u)
		break;
	}


	return (mean + sd * (v/u) );
}

double Rndm::Quantile_Beta(double alpha, double beta)
{

    double x = INFINITY;


    if ( alpha<0 )
    {
        cout << "alpha in quantile_beta is smaller than 0\n";
        throw ( "alpha is under 0" );
    }
    if ( beta<0)
    {
        cout << "beta in quantile_beta is smaller than 0\n";
        throw ( "beta is under 0" );
    }

    //boost::math::beta_distribution<> b(alpha,beta);
    double runi = Uniform();
    x = boost::math::ibeta_inv(alpha,beta,runi);
    //x = Rndm::Quantile_Beta(alpha,beta);
    //x = b.beta_detail::quantile(alpha,beta,runi);

    return x;
}
