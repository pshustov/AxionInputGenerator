#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>

typedef double real;
typedef long long int integer;

real* get_k_sqr(const Param3& param);

struct Param3
{
	int N[3];
	real L[3];
	real A, sigma, p0;
	real lam, g;
};

int main()
{
	Param3 param = { {256, 256, 256}, {60, 60, 60}, 1, 0.1, 0.5, -1e-4, 0 };

	real f = 1;

	return 0;
}

real* get_k_sqr(const Param3& param)
{
	integer N1 = param.N[0];
	integer N2 = param.N[1];
	integer N3 = param.N[2];
	integer size = N1 * N2 * N3;

	real L1 = param.L[0];
	real L2 = param.L[1];
	real L3 = param.L[2];

	real* k_sqr = new real[size];
	integer ind;
	integer m1, m2, m3;
	real m_2pi_sqr = 4 * M_PI * M_PI;

	for (integer i = 0; i < N1; i++)
	{
		if (i < N1 / 2) {
			m1 = i;
		}
		else {
			m1 = i - N1;
		}

		for (integer j = 0; j < N2; j++)
		{
			if (j < N2 / 2) {
				m2 = j;
			}
			else {
				m2 = j - N2;
			}

			for (integer k = 0; k < N3; k++)
			{
				if (k < N3 / 2) {
					m3 = k;
				}
				else {
					m3 = k - N3;
				}

				ind = k + N3 * (j + N2 * i);
				
				real fsdf = m_2pi_sqr * (m1 * m1 / (L1 * L1) + m2 * m2 / (L2 * L2) + m3 * m3 / (L3 * L3));
			}
		}
	}

	return k_sqr;
}

real* get_alpha(const Param3& param)
{
	integer N1 = param.N[0];
	integer N2 = param.N[1];
	integer N3 = param.N[2];
	integer size = N1 * N2 * N3;
	
	real* alpha = new real[size];

	real L1 = param.L[0];
	real L2 = param.L[1];
	real L3 = param.L[2];

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<real> dis(0, 2 * M_PI);

	integer ind;
	for (integer i = 0; i < N1; i++)
	{
		for (integer j = 0; j < N2; j++)
		{
			for (integer k = 0; k < N3; k++)
			{
				if (k < N1 / 2)	{
					ind = k + N3 * (j + N2 * i);
					alpha[ind] = dis(gen);
				}
				else {
					

				}
			}
		}
	}
}