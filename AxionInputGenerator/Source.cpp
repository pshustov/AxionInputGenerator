#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <complex>
#include <fftw3.h>	
#include <fstream>	
#include <string>
#include <iostream>

typedef double Real;
typedef long double RealDouble;
typedef long int Integer;
typedef std::complex<Real> Complex;

struct Param3;
Param3 getInitialParam(int argc, char* argv[]);
void printData(const Param3& param, const Real* f, const Real* g);
Real* get_g(const Param3& param, const Complex* F);
Real* get_f(const Param3& param, const Complex* F);
Complex* get_F(const Param3& param);
Real* get_alpha(const Param3& param);
Real* get_p_sqr(const Param3& param, Real p0);

struct Param3
{
	Integer N[3];
	Real L[3];
	Real f0, sigma, p0;
};

int main(int argc, char* argv[])
{
	std::cout << "Generating the initial conditions file\n";

	//Param3 param = { {64, 64, 64}, {60, 60, 60}, 1, 16, 0 };
	Param3 param = getInitialParam(argc, argv);
	Integer N = param.N[0] * param.N[1] * param.N[2];

	Complex* F = get_F(param);
	Real* f = get_f(param, F);
	Real* g = get_g(param, F);

	printData(param, f, g);

	delete[] F;
	fftw_free(f);
	fftw_free(g);

	std::cout << "File created\n";
	return 0;
}

void printData(const Param3& param, const Real* f, const Real* g)
{
	Integer N = param.N[0] * param.N[1] * param.N[2];
	Real V = param.L[0] * param.L[1] * param.L[2];
	
	char strFilename[256];
	sprintf_s(strFilename, "initial_A%05.2f_s%05.2f_D%d_N%d_L%06.2f.txt", param.f0, param.sigma, 3, N, V);
	std::ofstream out(strFilename);
	out.precision(12);

	out << param.N[0] << "\n" << param.N[1] << "\n" << param.N[2] << "\n";
	out << param.L[0] << "\n" << param.L[1] << "\n" << param.L[2] << "\n";
	out << param.f0 << "\n";
	out << param.sigma << "\n";
	out << param.p0 << "\n";

	for (size_t i = 0; i < N; i++)
	{
		out << f[i] << "\n";
	}

	for (size_t i = 0; i < N; i++)
	{
		out << g[i] << "\n";
	}
}

Real* get_g(const Param3& param, const Complex* F)
{
	Integer N = param.N[0] * param.N[1] * param.N[2];
	Integer Nred = param.N[0] * param.N[1] * (param.N[2] / 2 + 1);

	if (typeid(Real) != typeid(double) || typeid(Complex) != typeid(std::complex<double>)) {
		throw;
	}

	double E;
	double* kSqr = get_p_sqr(param, 0);
	double* g = (double*)fftw_malloc(sizeof(double) * N);
	fftw_complex* F_copy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nred);
	for (Integer i = 0; i < Nred; i++)
	{
		E = sqrt(1 + kSqr[i]);
		F_copy[i][0] = E * F[i].imag();
		F_copy[i][1] = -E * F[i].real();
	}

	fftw_plan plan;
	plan = fftw_plan_dft_c2r_3d(param.N[0], param.N[1], param.N[2], F_copy, g, FFTW_ESTIMATE);

	fftw_execute(plan);

	fftw_destroy_plan(plan);

	Real V = param.L[0] * param.L[1] * param.L[2];
	for (Integer i = 0; i < N; i++)
	{
		g[i] /= V;
	}

	fftw_free(F_copy);

	return g;
}

Real* get_f(const Param3& param, const Complex* F)
{
	Integer N = param.N[0] * param.N[1] * param.N[2];
	Integer Nred = param.N[0] * param.N[1] * (param.N[2] / 2 + 1);

	if (typeid(Real) != typeid(double) || typeid(Complex) != typeid(std::complex<double>)) {
		throw;
	}

	double* f = (double*)fftw_malloc(sizeof(double) * N);
	fftw_complex* F_copy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nred);
	for (Integer i = 0; i < Nred; i++)
	{
		F_copy[i][0] = F[i].real();
		F_copy[i][1] = F[i].imag();
	}

	fftw_plan plan;
	plan = fftw_plan_dft_c2r_3d(param.N[0], param.N[1], param.N[2], F_copy, f, FFTW_ESTIMATE);
	
	fftw_execute(plan);

	fftw_destroy_plan(plan);

	Real V = param.L[0] * param.L[1] * param.L[2];
	for (Integer i = 0; i < N; i++)
	{
		f[i] /= V;
	}

	fftw_free(F_copy);

	//normalization
	Real meanSqr = 0;
	for (size_t i = 0; i < N; i++) {
		meanSqr += f[i] * f[i];
	}
	meanSqr /= N;
	Real mul = param.f0 / sqrt(meanSqr);
	for (size_t i = 0; i < N; i++)
	{
		f[i] *= mul;
	}

	return f;
}

Complex* get_F(const Param3& param)
{
	Real A = param.f0;
	Real sigma = param.sigma;
	Real V = param.L[0] * param.L[1] * param.L[2];
	Integer N = param.N[0] * param.N[1] * param.N[2];
	Integer Nred = param.N[0] * param.N[1] * (param.N[2] / 2 + 1);

	Real m = 2 * sqrt(M_PI) / sigma;
	Real F0 = A * sqrt(V * m * m * m);

	Real* alpha = get_alpha(param);
	Real* pSqr = get_p_sqr(param, param.p0);
	Complex* F = new Complex[Nred];

	Complex imagUnit(0., 1.);

	for (Integer i = 0; i < Nred; i++)
	{
		F[i] = F0 * exp(-pSqr[i] / (2 * sigma * sigma)) * exp(imagUnit * alpha[i]);
	}

	delete[] alpha;
	delete[] pSqr;

	return F;
}

Real* get_p_sqr(const Param3& param, Real p0)
{
	Integer N1 = param.N[0];
	Integer N2 = param.N[1];
	Integer N3 = param.N[2];
	Integer N3red = param.N[2] / 2 + 1;
	Integer size = N1 * N2 * N3;
	Integer sizeRed = N1 * N2 * N3red;

	Real L1 = param.L[0];
	Real L2 = param.L[1];
	Real L3 = param.L[2];

	Real* pSqr = new Real[sizeRed];
	Real p;


	Integer indRed;
	Integer m1, m2, m3;
	Real m1Sqr, m2Sqr, m3Sqr;

	for (Integer i = 0; i < N1; i++)
	{
		if (i < N1 / 2) {
			m1 = i;
		}
		else {
			m1 = i - N1;
		}
		m1Sqr = m1;
		m1Sqr *= m1Sqr;

		for (Integer j = 0; j < N2; j++)
		{
			if (j < N2 / 2) {
				m2 = j;
			}
			else {
				m2 = j - N2;
			}
			m2Sqr = m2;
			m2Sqr *= m2Sqr;

			for (Integer k = 0; k < N3red; k++)
			{
				if (k < N3 / 2) {
					m3 = k;
				}
				else {
					m3 = k - N3;
				}
				m3Sqr = m3;
				m3Sqr *= m3Sqr;

				indRed = k + N3red * (j + N2 * i);
				
				p = 2 * M_PI *sqrt((m1Sqr / (L1 * L1) + m2Sqr / (L2 * L2) + m3Sqr / (L3 * L3)));
				pSqr[indRed] = (p - p0) * (p - p0);
			}
		}
	}

	return pSqr;
}

Real* get_alpha(const Param3& param)
{
	Integer N1 = param.N[0];
	Integer N2 = param.N[1];
	Integer N3 = param.N[2];
	Integer N3red = param.N[2] / 2 + 1;
	Integer size = N1 * N2 * N3;
	Integer sizeRed = N1 * N2 * N3red;
	
	Real* alpha = new Real[sizeRed];

	Real L1 = param.L[0];
	Real L2 = param.L[1];
	Real L3 = param.L[2];

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<Real> dis(0, 2 * M_PI);

	Integer ind, indr, indRed;
	Integer ir, jr, kr;
	for (Integer i = 0; i < N1; i++)
	{
		ir = (N1 - i) % N1;
		for (Integer j = 0; j < N2; j++)
		{
			jr = (N2 - j) % N2;
			for (Integer k = 0; k < N3red; k++)
			{
				kr = (N3 - k) % N3;

				ind = k + N3 * (j + N2 * i);
				indr = kr + N3 * (jr + N2 * ir);
				
				indRed = k + N3red * (j + N2 * i);

				if (ind == indr)
				{
					alpha[ind] = 0;
				}
				else
				{
					alpha[indRed] = dis(gen);
				}
			}
		}
	}

	return alpha;
}

Param3 getInitialParam(int argc, char* argv[])
{
	Param3 param = { {0, 0, 0}, {0, 0, 0}, 0, 0, 0 };;
	std::cout << "This program works only with 3D grids\n";
	if (argc < 2)
	{
		Integer N;
		Real L;
		std::cout << "Input N: ";
		std::cin >> N;
		param.N[0] = N; param.N[1] = N; param.N[2] = N;
		std::cout << "Input L: ";
		std::cin >> L;
		param.L[0] = L; param.L[1] = L; param.L[2] = L;
		std::cout << "Input f0: ";
		std::cin >> param.f0;
		std::cout << "Input sigma: ";
		std::cin >> param.sigma;
		std::cout << "Input p0: ";
		std::cin >> param.p0;
	}
	else
	{
		int N;
		Real L;

		for (int i = 1; i < argc; ++i) {
			std::string arg = argv[i];

			if ((arg == "-N")) {
				if (i + 1 < argc) {
					N = atoi(argv[++i]);
					param.N[0] = N; param.N[1] = N; param.N[2] = N;
				}
				else {
					std::cerr << "-N option requires one argument." << std::endl;
					throw;
				}
			}

			if ((arg == "-L")) {
				if (i + 1 < argc) {
					L = atof(argv[++i]);
					param.L[0] = L; param.L[1] = L; param.L[2] = L;
				}
				else {
					std::cerr << "-N option requires one argument." << std::endl;
					throw;
				}
			}

			if ((arg == "-f0")) {
				if (i + 1 < argc) {
					param.f0 = atof(argv[++i]);
				}
				else {
					std::cerr << "-N option requires one argument." << std::endl;
					throw;
				}
			}

			if ((arg == "-s")) {
				if (i + 1 < argc) {
					param.sigma = atof(argv[++i]);
				}
				else {
					std::cerr << "-N option requires one argument." << std::endl;
					throw;
				}
			}

			if ((arg == "-p0")) {
				if (i + 1 < argc) {
					param.p0 = atof(argv[++i]);
				}
				else {
					std::cerr << "-N option requires one argument." << std::endl;
					throw;
				}
			}
		}

	}

	std::cout << "Initial parameters:\n";
	std::cout << "N: [" << param.N[0] << ", " << param.N[1] << ", " << param.N[2] << "]\n";
	std::cout << "L: [" << param.L[0] << ", " << param.L[1] << ", " << param.L[2] << "]\n";
	std::cout << "f0: " << param.f0 << std::endl;
	std::cout << "sigma: " << param.sigma << std::endl;
	std::cout << "p0: " << param.p0 << std::endl;

	return param;
}