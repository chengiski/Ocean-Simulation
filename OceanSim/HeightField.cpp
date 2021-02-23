
#include "HeightField.h"

#include <fftw3.h>

HeightField::~HeightField() {
	delete[] k_vector;
	delete[] k_vector_normalized;
	delete[] h0_map;
	delete[] h0_conj_map;
	delete[] h_bar_map;
	delete[] slope_x;
	delete[] slope_y;
	delete[] Dx;
	delete[] Dy;
	// Jacobian term 
	delete[] Jxx;
	delete[] Jxy;
	delete[] Jyy;
	// displace map and normal map
	delete[] displaceMap;
	delete[] normalMap;
	delete[] jacobianMap;
}

void HeightField::initialize() {
	// mesh board resolution
	xres = 256;
	yres = 256;
	num = xres * yres;
	// real world scale
	Lx = 1000;
	Ly = 1000;
	l = 0.05;
	// global amplitude
	A = 3e-7f;
	// wind speed
	v = 30;
	// wind vector
	w = glm::vec2(1.0, 1.0);
	// lambda
	lambda = 1.0;
	// pseudo random 
	generator.seed(time(NULL));
	// new storage for data map
	k_vector = new glm::vec2[num];
	k_vector_normalized = new glm::vec2[num];
	h0_map = new std::complex<float>[num];
	h0_conj_map = new std::complex<float>[num];
	h_bar_map = new std::complex<float>[num];
	Dx = new std::complex<float>[num];
	Dy = new std::complex<float>[num];
	slope_x = new std::complex<float>[num];
	slope_y = new std::complex<float>[num];
	// Jacobian term
	Jxx = new std::complex<float>[num];
	Jxy = new std::complex<float>[num];
	Jyy = new std::complex<float>[num];
	// output map
	displaceMap = new glm::vec3[num];
	normalMap = new glm::vec3[num];
	jacobianMap = new glm::vec3[num];
	for (int m = 0; m < yres; m++) {
		for (int n = 0; n < xres; n++) {
			int index = m * xres + n;
			k_vector[index] = glm::vec2(2 * float(std::_Pi) * (n - xres / 2) / Lx, 2 * float(std::_Pi) * (m - yres / 2) / Ly);
			if (glm::length(k_vector_normalized[index]) != 0) {
				k_vector_normalized[index] = glm::normalize(k_vector[index]);
			}
			else {
				k_vector_normalized[index] = k_vector[index];
			}

			h0_map[index] = h0(k_vector[index]);
			h0_conj_map[index] = std::conj(h0(k_vector[index]));
		}
	}
}


void HeightField::initialize(int resx, int resy, int Lx, int Ly, float A, float v, glm::vec2 w, float lambda) {
	// mesh board resolution
	xres = resx;
	yres = resy;
	num = xres * yres;
	// real world scale
	this->Lx = Lx;
	this->Ly = Ly;
	l = 0.1;
	// global amplitude
	this->A = A;
	// wind speed
	this->v = v;
	// wind vector
	this->w = w;
	// lambda
	this->lambda = lambda;
	// pseudo random 
	generator.seed(time(NULL));
	// new storage for data map
	k_vector = new glm::vec2[num];
	k_vector_normalized = new glm::vec2[num];
	h0_map = new std::complex<float>[num];
	h0_conj_map = new std::complex<float>[num];
	h_bar_map = new std::complex<float>[num];
	Dx = new std::complex<float>[num];
	Dy = new std::complex<float>[num];
	slope_x = new std::complex<float>[num];
	slope_y = new std::complex<float>[num];
	// Jacobian term
	Jxx = new std::complex<float>[num];
	Jxy = new std::complex<float>[num];
	Jyy = new std::complex<float>[num];
	// output map
	displaceMap = new glm::vec3[num];
	normalMap = new glm::vec3[num];

	for (int m = 0; m < yres; m++) {
		for (int n = 0; n < xres; n++) {
			int index = m * xres + n;
			k_vector[index] = glm::vec2(2 * float(std::_Pi) * (n - xres / 2) / Lx, 2 * float(std::_Pi) * (m - yres / 2) / Ly);
			if (glm::length(k_vector_normalized[index]) != 0) {
				k_vector_normalized[index] = glm::normalize(k_vector[index]);
			}
			else {
				k_vector_normalized[index] = k_vector[index];
			}

			h0_map[index] = h0(k_vector[index]);
			h0_conj_map[index] = std::conj(h0(k_vector[index]));
		}
	}
}

float HeightField::philipsSpectrum(glm::vec2 K) {
	if (K == glm::vec2(0.0f, 0.0f))
		return 0.0f;
	float L = v * v / gravity;
	// (k_head * w_head)^2
	float cos_part = glm::dot(glm::normalize(K), glm::normalize(w));
	cos_part = std::pow(cos_part, 2);
	// exp part
	// k^2
	float k = glm::length(K);
	// e^{-1/(kL)^2}
	float exp_part = std::exp(-1 / std::pow(k*L, 2));
	float phi_spec = A * exp_part*cos_part / std::pow(k, 4);
	phi_spec *= exp(-k * k*l*l);
	return phi_spec;
}
std::complex<float> HeightField::h0(glm::vec2 k) {
	float real = distribution(generator);
	float img = distribution(generator);
	float phi_part = std::sqrtf(0.5f)*std::sqrtf(philipsSpectrum(k));
	real *= phi_part;
	img *= phi_part;
	std::complex<float> res(real, img);
	return res;
}
float HeightField::dispersionRelation(glm::vec2 k) {
	//float wave_length = glm::length(float(std::_Pi) / k);
	// deep water
	float omega = std::sqrtf(glm::length(k)*gravity);
	return omega;
}
std::complex<float> HeightField::exponent_part(float omega, float t) {
	//apply euler formula
	float omega_t = omega * t;
	std::complex<float> res(std::cosf(omega_t), std::sinf(omega_t));
	return res;
}
std::complex<float> HeightField::h_bar(std::complex<float> h0, std::complex<float> h0_conj, float t, glm::vec2 k) {
	std::complex<float> exp_part, res;
	//h_twiddle = h0*exponent_part(dispersionRelation(k), t);
	exp_part = exponent_part(dispersionRelation(k), t);
	res = h0 * exp_part + h0_conj * std::conj(exp_part);
	return res;
}


void HeightField::buildField(float t) {
	// height, slope_x,slope_y,Dx,Dy
	fftwf_complex *in_z, *in_slope_x, *in_slope_y, *in_Dx, *in_Dy;
	fftwf_complex *out_z, *out_slope_x, *out_slope_y, *out_Dx, *out_Dy;
	// jacobian matrix
	fftwf_complex *in_Jxx, *in_Jxy, *in_Jyy;
	fftwf_complex *out_Jxx, *out_Jxy, *out_Jyy;
	// fftw plan 
	fftwf_plan p_z, p_slope_x, p_slope_y, p_Dx, p_Dy;
	fftwf_plan p_Jxx, p_Jxy, p_Jyy;

	// calculate 
	// h_bar_map 
	// slope_x slope_y
	// Dx Dy
	// Jxx Jxy Jyy
	for (int m = 0; m < yres; m++) {
		for (int n = 0; n < xres; n++) {
			int index = m * xres + n;
			// h_bar map
			h_bar_map[index] = h_bar(h0_map[index], h0_conj_map[index], t, k_vector[index]);
			slope_x[index] = std::complex<float>(0, k_vector[index].x)*h_bar_map[index];
			slope_y[index] = std::complex<float>(0, k_vector[index].y)*h_bar_map[index];
			// Dx, Dy 
			Dx[index] = std::complex<float>(0, -k_vector_normalized[index].x)*h_bar_map[index];
			Dy[index] = std::complex<float>(0, -k_vector_normalized[index].y)*h_bar_map[index];
			//jacobian term Jxx Jxy Jyy, (Jxy = Jyx)
			Jxx[index] = std::complex<float>(1, 0) + lambda * k_vector[index].x*Dx[index];
			Jxy[index] = std::complex<float>(1, 0) + lambda * k_vector[index].x*Dy[index];
			Jyy[index] = std::complex<float>(1, 0) + lambda * k_vector[index].y*Dy[index];
		}
	}

	// data convert into fftw input
	in_z = (fftwf_complex*)h_bar_map;
	in_slope_x = (fftwf_complex*)slope_x;
	in_slope_y = (fftwf_complex*)slope_y;
	in_Dx = (fftwf_complex*)Dx;
	in_Dy = (fftwf_complex*)Dy;
	in_Jxx = (fftwf_complex*)Jxx;
	in_Jxy = (fftwf_complex*)Jxy;
	in_Jyy = (fftwf_complex*)Jyy;
	// open new space for output 
	out_z = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num);
	out_slope_x = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num);
	out_slope_y = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num);
	out_Dx = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num);
	out_Dy = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num);
	out_Jxx = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num);
	out_Jxy = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num);
	out_Jyy = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * num);

	// fftw plan
	p_z = fftwf_plan_dft_2d(xres, yres, in_z, out_z, FFTW_BACKWARD, FFTW_ESTIMATE);
	p_slope_x = fftwf_plan_dft_2d(xres, yres, in_slope_x, out_slope_x, FFTW_BACKWARD, FFTW_ESTIMATE);
	p_slope_y = fftwf_plan_dft_2d(xres, yres, in_slope_y, out_slope_y, FFTW_BACKWARD, FFTW_ESTIMATE);
	p_Dx = fftwf_plan_dft_2d(xres, yres, in_Dx, out_Dx, FFTW_BACKWARD, FFTW_ESTIMATE);
	p_Dy = fftwf_plan_dft_2d(xres, yres, in_Dy, out_Dy, FFTW_BACKWARD, FFTW_ESTIMATE);
	p_Jxx = fftwf_plan_dft_2d(xres, yres, in_Jxx, out_Jxx, FFTW_BACKWARD, FFTW_ESTIMATE);
	p_Jxy = fftwf_plan_dft_2d(xres, yres, in_Jxy, out_Jxy, FFTW_BACKWARD, FFTW_ESTIMATE);
	p_Jyy = fftwf_plan_dft_2d(xres, yres, in_Jyy, out_Jyy, FFTW_BACKWARD, FFTW_ESTIMATE);

	fftwf_plan pp[10];
	pp[0] = p_z;
	pp[1] = p_slope_x;
	pp[2] = p_slope_y;
	pp[3] = p_Dx;
	pp[4] = p_Dy;
	pp[5] = p_Jxx;
	pp[6] = p_Jxy;
	pp[7] = p_Jyy;

#pragma omp parallel for
	for (int i = 0; i < 8; i++)
		fftwf_execute(pp[i]);

	// execute plan
	//fftwf_execute(p_z);
	//fftwf_execute(p_slope_x);
	//fftwf_execute(p_slope_y);
	//fftwf_execute(p_Dx);
	//fftwf_execute(p_Dy);
	//fftwf_execute(p_Jxx);
	//fftwf_execute(p_Jxy);
	//fftwf_execute(p_Jyy);
	// update fft result
	for (int n = 0; n < xres; n++) {
		for (int m = 0; m < yres; m++) {
			int index = m * xres + n;
			float sign = 1;
			// (-1)^(m+n)
			if ((m + n) % 2)
				sign = -1;
			normalMap[index] = glm::normalize(glm::vec3(
				sign * out_slope_x[index][0],
				-1,
				sign * out_slope_y[index][0]));

			displaceMap[index] = glm::vec3(
				(n - xres / 2) * Lx / xres - sign * lambda * out_Dx[index][0],
				sign * out_z[index][0],
				(m - yres / 2) * Ly / yres - sign * lambda * out_Dy[index][0]);
			jacobianMap[index] = glm::vec3(
				out_Jxx[index][0],
				out_Jyy[index][0],
				out_Jxy[index][0]
			);
		}
	}


	// destroy plan
	fftwf_destroy_plan(p_z);
	fftwf_destroy_plan(p_slope_x);
	fftwf_destroy_plan(p_slope_y);
	fftwf_destroy_plan(p_Dx);
	fftwf_destroy_plan(p_Dy);
	fftwf_destroy_plan(p_Jxx);
	fftwf_destroy_plan(p_Jxy);
	fftwf_destroy_plan(p_Jyy);
	// destroy dynamic memory
	fftwf_free(out_z);
	fftwf_free(out_slope_x);
	fftwf_free(out_slope_y);
	fftwf_free(out_Dx);
	fftwf_free(out_Dy);
	fftwf_free(out_Jxx);
	fftwf_free(out_Jxy);
	fftwf_free(out_Jyy);
}