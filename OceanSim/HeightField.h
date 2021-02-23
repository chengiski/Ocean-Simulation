#define gravity 9.796

#include <complex>
#include <random>
#include <ctime>
#include <glm/glm.hpp>

class HeightField
{
public:
	// mesh board resolution
	int xres, yres;
	int num;
	// real world size
	int Lx, Ly;
	// random number for gaussian distribution
	std::default_random_engine generator;
	std::normal_distribution<double> distribution;
	// Philips Spectrum Parameters
	float l;
	// global amplitude
	float A;
	// wind direction vector
	glm::vec2 w;
	// wind speed
	float v;
	// output map
	glm::vec3* displaceMap;
	glm::vec3* normalMap;
	glm::vec3* jacobianMap;
	// k vector input
	glm::vec2* k_vector;
	glm::vec2* k_vector_normalized;
	// lambda
	float lambda;
	/* map for fft */
	// store initial h0 and conj map data
	std::complex<float> *h0_map = NULL;
	std::complex<float> *h0_conj_map = NULL;
	// h_bar map
	std::complex<float> *h_bar_map = NULL;
	// D_x
	std::complex<float> *Dx = NULL;
	// D_y
	std::complex<float> *Dy = NULL;
	// slope_x
	std::complex<float> *slope_x = NULL;
	// slope_y
	std::complex<float> *slope_y = NULL;
	// Jacobian term
	std::complex<float> *Jxx = NULL;
	std::complex<float> *Jyy = NULL;
	std::complex<float> *Jxy = NULL;
	~HeightField();
	// initialize the object with predefined parameters
	void initialize();
	// overloading initializing
	void initialize(int resx, int resy, int Lx, int Ly, float A, float v, glm::vec2 w, float lambda);
	// Philips Spectrum
	float philipsSpectrum(glm::vec2 K);
	// h0
	std::complex<float> h0(glm::vec2 k);
	// dispersion relationship
	float dispersionRelation(glm::vec2 k);
	// apply euler formula
	std::complex<float> exponent_part(float omega, float t);
	// h_bar = h0*exponet_part + conj(h0*exponent_part)
	std::complex<float> h_bar(std::complex<float> h0, std::complex<float> h0_conj, float t, glm::vec2 k);
	void buildField(float t);
};

