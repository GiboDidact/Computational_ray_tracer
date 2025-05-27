#pragma once
#include "../pch.h"
#include "../ThirdParty/pbrv4/helpers.h"
#include "../ThirdParty/pbrv4/samplers.h"

/*
color space conversion:
we need code to convert all the color space stuff, spectral to XYZ, then XYZ to rgb, srgb, etc
even going from RGB-->spectral
most of it will be known functions and integrating

white balance


filters:
going to need filter for cameras if I want to filter the data from the pixel

sampling:
going to need to sample areas on a pixel sensor, directions on a point of the sensor or directions at a surface
also need to sample wavelengths of the spectral functions, and the conversion spaces

figure out what sampling I want to do, to lower variance, to get the most important rays with known information

monte carlos:
monte carlos integrators, estimators for all these integrals

*/

//do some monte carlos integrtation programming with variance reduction techniques

//program some functions to generate sampling for pdf/pmf from probability theory, some common techniques

//create own distributions for common scenarios and use cases for ray tracing and program their sample functions

//program some filtering algorithms

//figure out how to get denoising


/*
The problem we need to solve is just to have many computational solutions to sample from any random distribution
we need it to be accurate, and efficient
ultaimtely we have some domain of the integral for the monte carlos, and we need to just pick samples, like in statistics we need samples from the population
so that we can estimate the correct value. for monte carlos you need to generate a random variable, so create an algorithm that does that

theres many ways to compute distributions:
. One way is just given an analytic function pmf/pdf generate samples on a computer
. Quasi-random sequences, generating those and using those
. algorithmic methods, creating your own random process and sampling from the random process through an algorithm
. Theres many other ways to create a distribution and sample from it

to verify you can keep sampling and generate a histogram and the histogram should approach the correct pdf/pmf
or analytically prooving that the returned random variable has the distribution you wanted
*/

//multidimensional general I dont have anything: can do rejection method, can try the marginal integral method to take one variable out and the other is
//a conditional probability, can try to do inversion but I need to see how that works if you can do it at all, or just give a discrete grid and interpolate

//**marginal method is best, you just itnegrate out X = integral(p(x,y)dy), then sample that then do Y=p(x,y)/X
//you can solve for p(x) then use the 1d method to get it, then same with y, 


inline float VisibleWavelengthsPDF(float lambda) {
	if (lambda < 360 || lambda > 830)
		return 0;
	return 0.0039398042f / std::pow(std::cosh(0.0072f * (lambda - 538)),2);
}

inline float SampleVisibleWavelengths(float u) {
	return 538 - 138.888889f * std::atanh(0.85691062f - 1.82750197f * u);
}

//**for architecture I probably need a class that encapsalates a general pdf/pmf with its sampler, that way I can pass it to monte carlos and it knows
//how to sample the pdf and can access everything together

class Discrete_Inversion_Sampler
{
public:
	//this one just takes in the probability set, and it doesn't even need to be normalizd, we return index so we don't care about the exact x values
	//that'l be handled outside
	Discrete_Inversion_Sampler(std::vector<float> pmf_weights)
	{
		//create cdf
		weight_sum = 0;
		for (int i = 0; i < pmf_weights.size(); i++)
		{
			weight_sum += pmf_weights[i];
			cdf.push_back(weight_sum);
		}

		/*for (int i = 0; i < pmf_weights.size(); i++)
		{
			std::cout << cdf[i] << std::endl;
		}
		std::cout << "weight_sum: " << weight_sum << std::endl;
		*/
	}
	~Discrete_Inversion_Sampler() = default;

	int Sample() const
	{
		float U = Helper::GetRandomNumber(0, 1) * weight_sum;

		//find index, n<U<n+1, binary search since sorted list, orders of magnitude faster
		int index = -1;
		int low = 0;
		int high = cdf.size() - 1;
		while (low <= high)
		{
			int mid = low + (high - low) / 2.0f;

			if (cdf[mid] <= U && U < cdf[mid + 1])
			{
				index = mid + 1;
				break;
			}

			if (cdf[mid] < U)
				low = mid + 1;
			else
				high = mid - 1;

			//technically index 0 could be .3, so if U is <.3 it doesnt pick it up, this handles that case if U<cdf[0] then index = 0
			if (high == -1)
			{
				index = 0;
				break;
			}
		};
		if (index == -1)
		{
			std::cout << "couldn't find index\n";
			return 0;
		}

		return index;
	}

	std::vector<float> generateTestSampleHistogram(std::vector<float> x_values, int sample_count, float scaling_factor = 1.0f) const
	{
		//create histogram
		std::vector<float>index_samples;
		index_samples.resize(cdf.size(), 0);
		for (int i = 0; i < sample_count; i++)
		{
			int index_i = Sample();
			index_samples[index_i] += 1;
		}

		//normalize
		for (int i = 0; i < index_samples.size(); i++)
		{
			index_samples[i] = scaling_factor * (index_samples[i] / (float)sample_count);
		}

		std::vector<float> histogram_samples; //come in x,y pairs
		for (int i = 0; i < index_samples.size(); i++)
		{
			if (i == 0)
			{
				histogram_samples.push_back(x_values[i]);
				histogram_samples.push_back(0);
			}
			else
			{
				histogram_samples.push_back(x_values[i]);
				histogram_samples.push_back(index_samples[i - 1]);
			}
			histogram_samples.push_back(x_values[i]);
			histogram_samples.push_back(index_samples[i]);
		}
		histogram_samples.push_back(x_values[x_values.size() - 1] + 1);
		histogram_samples.push_back(index_samples[index_samples.size() - 1]);
		histogram_samples.push_back(x_values[x_values.size() - 1] + 1);
		histogram_samples.push_back(0);

		return histogram_samples;
	}
private:
	std::vector<float> cdf;
	float weight_sum;
};


//1D
//sample linear - the general function is: f(x) = (1-x)a + xb, x=[0,1], interpolates between a and b, a,b>=0

inline float LinearPdf(float x, float a, float b)
{
	if (x < 0 || x > 1)
		return 0;
	return 2 * pbrt::Lerp(x, a, b) / (a + b);
}

//u ~ uniform(0,1), a,b>=0
inline float SampleLinear(float a, float b)
{
	//debug check a,b>=0
	float u = Helper::GetRandomNumber(0, 1);
	if (u == 0 && a == 0) return 0;
	float x = (u * (a + b)) / (a + std::sqrtf(pbrt::Lerp(u, a * a, b * b)));
	return std::min(x, OneMinusEpsilon);
}

inline float SampleLinear(float u, float a, float b)
{
	//debug check a,b>=0
	if (u == 0 && a == 0) return 0;
	float x = (u * (a + b)) / (a + std::sqrtf(pbrt::Lerp(u, a * a, b * b)));
	return std::min(x, OneMinusEpsilon);
}

inline float InvertLinearSample(float x, float a, float b)
{
	return x * (a * (2 - x) + b * x) / (a + b);
}

//sample tent, radius is r
inline float SampleTent(float r)
{
	float u = Helper::GetRandomNumber(0, 1);
	Discrete_Inversion_Sampler coinflip({ .5,.5 });
	if (coinflip.Sample() == 0)
		return -r + r * SampleLinear(u, 0, 1);
	else
		return r * SampleLinear(u, 1, 0);
}
inline float SampleTent(float u, float r)
{
	Discrete_Inversion_Sampler coinflip({ .5,.5 });
	if (coinflip.Sample() == 0)
		return -r + r * SampleLinear(u, 0, 1);
	else
		return r * SampleLinear(u, 1, 0);
}

inline float TentPDF(float x, float r)
{
	if (std::abs(x) >= r)
		return 0;
	return 1 / r - std::abs(x) / (r * r);
}

inline float InvertTentSample(float x, float r)
{
	if (x <= 0)
		return (1 - InvertLinearSample(-x / r, 1, 0)) / 2.0f;
	else
		return 0.5f + InvertLinearSample(x / r, 1, 0) / 2.0f;
}

inline float SampleTentSumMethod(float r)
{
	//Z=X+Y creates a tent, X,Y~uniform(a,b)
	float X = Helper::GetRandomNumber(-r/2.0f, r/2.0f);
	float Y = Helper::GetRandomNumber(-r/2.0f, r/ 2.0f);
	return X + Y;
}

//exponential

//e^-ax, x[0,infinity]
inline float SampleExponential(float a)
{
	//debug make sure a!= 0
	float u = Helper::GetRandomNumber(0, 1);
	return -std::log(1 - u) / a;
}

inline float ExponentialPDF(float x, float a)
{
	return a * std::exp(-a * x);
}

inline float InvertExponentialSample(float x, float a)
{
	return 1 - std::exp(-a * x);
}

//normal
inline float GaussianIntegral(float x0, float x1, float mu = 0, float sigma = 1)
{
	float sigmaRoot2 = sigma * float(1.414213562373095);
	return 0.5f * (std::erf((mu - x0) / sigmaRoot2) -
		std::erf((mu - x1) / sigmaRoot2));
}
/*inline float Gaussian(float x, float mu = 0, float sigma = 1) {
	return 1 / std::sqrt(2 * std::numbers::pi * sigma * sigma) *
		FastExp(-((x - mu)*(x-mu)) / (2 * sigma * sigma));
}

inline float NormalPDF(float x, float mu = 0, float sigma = 1)
{
	return Gaussian(x, mu, sigma);
}*/
//inverse

inline float SampleNormal(float mu = 0, float sigma = 1)
{
	float u = Helper::GetRandomNumber(0, 1);
	return mu + Sqrt2 * sigma * pbrt::ErfInv(2 * u - 1);
}

//logistic
inline float LogisticPDF(float x, float s)
{
	x = std::abs(x);
	return std::exp(-x / s) / (s * std::sqrt(1 + std::exp(-x / s)));
}

inline float SampleLogistic(float s)
{
	float u = Helper::GetRandomNumber(0, 1);
	return -s * std::log(1 / u - 1);
}

inline float InvertLogisticSample(float x, float s)
{
	return 1 / (1 + std::exp(-x / s));
}


//2D
//sample bilinear, f(x,y) = (1-x)(1-y)w0 + x(1-y)w1 + y(1-x)w2 + xyw3
inline float BilinearPDF(glm::vec2 p, std::vector<float> w)
{
	if (p.x < 0 || p.x > 1 || p.y < 0 || p.y > 1)
		return 0;
	if (w[0] + w[1] + w[2] + w[3] == 0)
		return 1;
	return 4 * ((1 - p[0]) * (1 - p[1]) * w[0] + p[0] * (1 - p[1]) * w[1] +
		(1 - p[0]) * p[1] * w[2] + p[0] * p[1] * w[3]) /
		(w[0] + w[1] + w[2] + w[3]);
}

inline glm::vec2 InvertBilinearSample(glm::vec2 p, std::vector<float> w)
{
	return { InvertLinearSample(p.x, pbrt::Lerp(p.y,w[0],w[2]), pbrt::Lerp(p.y,w[1],w[3])),
		InvertLinearSample(p.y,w[0] + w[1],w[2] + w[3]) };
}

//x[0,1], y[0,1]
inline glm::vec2 SampleBilinear(std::vector<float> w)
{
	//debug w size == 4
	glm::vec2 u(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1));
	glm::vec2 p;
	p.y = SampleLinear(u.y, w[0] + w[1], w[2] + w[3]);
	p.x = SampleLinear(u.x, pbrt::Lerp(p.y, w[0], w[2]), pbrt::Lerp(p.y, w[1], w[3]));
	
	return p;
}

//uniformdiskpolar
inline glm::vec2 SampleUniformDiskPolar()
{
	glm::vec2 u(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1));
	float r = std::sqrt(u.x);
	float theta = 2 * std::numbers::pi * u.y;
	return { r * std::cos(theta), r * std::sin(theta) };
}

inline glm::vec2 SampleDiskNaive()
{
	glm::vec2 u(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1));
	float r = u.x;
	float theta = 2 * std::numbers::pi * u.y;
	return { r * std::cos(theta), r * std::sin(theta) };
}

inline glm::vec2 RejectionSampleDisk()
{
	glm::vec2 p;
	do
	{
		p.x = 1 - 2*Helper::GetRandomNumber(0,1);
		p.y = 1 - 2*Helper::GetRandomNumber(0,1);
	} while (p.x*p.x + p.y*p.y > 1);
	return p;
}

//uniformdiskconcentric
inline glm::vec2 SampleUniformDiskConcentric(glm::vec2 u)
{
	//glm::vec2 u(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1));
	
	// Map _u_ to $[-1,1]^2$ and handle degeneracy at the origin
	glm::vec2 uOffset = glm::vec2(2,2) * u - glm::vec2(1, 1);
	if (uOffset.x == 0 && uOffset.y == 0)
		return { 0, 0 };

	// Apply concentric mapping to point
	float theta, r;
	if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
		r = uOffset.x;
		theta = PiOver4 * (uOffset.y / uOffset.x);
	}
	else {
		r = uOffset.y;
		theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
	}
	return r * glm::vec2(std::cos(theta), std::sin(theta));
}

//uniform hemisphere
inline glm::vec3 SampleUniformHemisphere()
{
	glm::vec2 u(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1));
	float z = u.x;
	float r = Helper::SafeSqrt(1 - z * z);
	float phi = 2 * std::numbers::pi * u.y;

	return { r * std::cos(phi), r * std::sin(phi), z };
}

inline float UniformHemispherePDF() { return Inv2Pi; };

inline glm::vec2 InvertUniformHemisphereSample(glm::vec3 w)
{
	float phi = std::atan2(w.y, w.x);
	if (phi < 0)
		phi += 2 * std::numbers::pi;
	return glm::vec2(w.z, phi / (2 * std::numbers::pi));
}

//uniform sphere
static glm::vec3 SampleUniformSphere()
{
	glm::vec2 u(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1));
	float z = 1 - 2*u.x;
	float r = Helper::SafeSqrt(1 - z * z);
	float phi = 2 * std::numbers::pi * u.y;

	return { r * std::cos(phi), r * std::sin(phi), z };
}

inline float UniformSpherePDF() { return Inv4Pi; }

inline glm::vec2 InvertUniformSphereSample(glm::vec3 w)
{
	float phi = std::atan2(w.y, w.x);
	if (phi < 0)
		phi += 2 * std::numbers::pi;
	return glm::vec2((1-w.z)/2.0f, phi / (2 * std::numbers::pi));
}


//cosine weighted hemisphere
inline glm::vec3 SampleCosineHemisphere(glm::vec2 u)
{
	glm::vec2 d = SampleUniformDiskConcentric(u);
	float z = Helper::SafeSqrt(1 - d.x * d.x - d.y * d.y);
	return glm::vec3(d.x, d.y, z);
}

inline float CosineHemispherePDF(float cosTheta)
{
	return cosTheta * InvPi;
}

inline glm::vec2 InvertCosineHemisphereSample(glm::vec3 w)
{
	//return InvertUniformDiskConcentricSample({ w.x,w.y });
}

//cone
inline float UniformConePDF(float cosThetaMax)
{
	return 1 / (2 * std::numbers::pi * (1 - cosThetaMax));
}

//a cone at (0,0,0) pointing (0,0,1), with angle cosThetaMax
inline glm::vec3 SampleUniformCone(float cosThetaMax)
{
	glm::vec2 u(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1));
	float cosTheta = (1 - u.x) + u.x * cosThetaMax;
	float sinTheta = Helper::SafeSqrt(1 - cosTheta * cosTheta);
	float phi = u.y * 2 * std::numbers::pi;
	return Helper::SphericalDirection(sinTheta, cosTheta, phi);
}

//SHAPES - do later
//light: sphere, cylinder, disk, triangle

inline std::function<float(float)> normalizepdf(std::function<float(float)> pdf, float a, float b)
{
	//all you need is the area under the curve and then you divide the pdf by area
	int N = 10000;
	float area = 0;
	float delta_x = (b - a) / (float)N;
	for (int n = 0; n < N; n++)
	{
		//reiman method, just add the rectangle area: width*height = delta_x*pdf(x_n)
		float current_x = std::clamp(a + delta_x * n, a, b);
		//reiman method, just add the rectangle area: width*height = delta_x*pdf(x_n)
		area += delta_x * pdf(current_x);
	}
	std::cout << area << std::endl;
	auto normalized_pdf = [=](float x) -> float { return pdf(x) / area; };
	
	return normalized_pdf;
}

inline int Discrete_Rejection_Sample(std::vector<float> pmf_weights, float max_height = 0)
{
	//find max height
	if (max_height == 0)
	{
		max_height = std::numeric_limits<float>::min();
		for (int i = 0; i < pmf_weights.size(); i++)
		{
			max_height = std::max(max_height, pmf_weights[i]);
		}
	}

	//box dimension
	float x_0 = 0;
	float x_1 = pmf_weights.size(); //add extra one for last histogram right side
	float y_0 = 0;
	float y_1 = max_height;

	//check efficiency area ratios
	if (false)
	{
		float Area_box = (x_1 - x_0) * (y_1 - y_0);
		float Area_histogram = 0;
		for (int i = 0; i < pmf_weights.size(); i++)
		{
			Area_histogram += 1 * pmf_weights[i];
		}
		float efficiency = Area_histogram / Area_box;
		std::cout <<"max height: "<<max_height << " box area: " << Area_box << " histogram area: " << Area_histogram << " efficiency: " << efficiency << std::endl;
	}

	//now pick random uniform variables X~uniform(x_0,x_1), Y~uniform(y_0,y_1) and see if its inside the histogram or not
	int tries = 0;
	while (true)
	{
		tries++;
		float X = Helper::GetRandomNumber(x_0, x_1);
		float Y = Helper::GetRandomNumber(y_0, y_1);

		float x_index = std::floor(X);
		float height = pmf_weights[x_index];

		if (Y <= height) {
			//std::cout << "tries: " << tries << std::endl;
			return x_index;
		}
	}
}

inline std::vector<float> generateRejectionTestSampleHistogram(std::vector<float> x_values, std::vector<float> pmf_weights, int sample_count, float scaling_factor = 1.0f)
{
	//get max height for optimization
	float max_height = std::numeric_limits<float>::min();
	for (int i = 0; i < pmf_weights.size(); i++)
	{
		max_height = std::max(max_height, pmf_weights[i]);
	}
	//create histogram
	std::vector<float>index_samples;
	index_samples.resize(pmf_weights.size(), 0);
	for (int i = 0; i < sample_count; i++)
	{
		int index_i = Discrete_Rejection_Sample(pmf_weights, max_height);
		index_samples[index_i] += 1;
	}

	//normalize
	for (int i = 0; i < index_samples.size(); i++)
	{
		index_samples[i] = scaling_factor * (index_samples[i] / (float)sample_count);
	}

	std::vector<float> histogram_samples; //come in x,y pairs
	for (int i = 0; i < index_samples.size(); i++)
	{
		if (i == 0)
		{
			histogram_samples.push_back(x_values[i]);
			histogram_samples.push_back(0);
		}
		else
		{
			histogram_samples.push_back(x_values[i]);
			histogram_samples.push_back(index_samples[i - 1]);
		}
		histogram_samples.push_back(x_values[i]);
		histogram_samples.push_back(index_samples[i]);
	}
	histogram_samples.push_back(x_values[x_values.size() - 1] + 1);
	histogram_samples.push_back(index_samples[index_samples.size() - 1]);
	histogram_samples.push_back(x_values[x_values.size() - 1] + 1);
	histogram_samples.push_back(0);

	return histogram_samples;
}

//continuous reject sampler
inline glm::vec2 Continuous_Rejection_Sample2D(std::function<float(glm::vec2)> _pdf, float x_1, float x_2, float y_1, float y_2, float max_height = 0)
{
	int N = 100;
	float buffer = 1.0f;
	float delta_x = (x_2 - x_1) / (float)N;
	float delta_y = (y_2 - y_1) / (float)N;
	if (max_height == 0)
	{
		max_height = std::numeric_limits<float>::min();
		for (float x = x_1; x <= x_2; x += delta_x)
		{
			for (float y = y_1; y <= y_2; y += delta_y)
			{
				max_height = std::max(max_height, _pdf(glm::vec2(x, y)));
			}
		}
		max_height += buffer;
	}

	float z_1 = 0;
	float z_2 = max_height;

	while (true)
	{
		float X = Helper::GetRandomNumber(x_1, x_2);
		float Y = Helper::GetRandomNumber(y_1, y_2);
		float Z = Helper::GetRandomNumber(z_1, z_2);

		if (Z <= _pdf(glm::vec2(X, Y)))
			return { X,Y };
	}
}

//continuous rejection histogram maker
inline float Continuous_Rejection_Sample(std::function<float(float)> _pdf, float a, float b, float max_height = 0)
{
	int N = 100;
	float buffer = 1.0f;
	float delta_x = (b - a) / (float)N;
	if (max_height == 0)
	{
		//we have to approximate the max height, discretize the function and pick the max one plus some buffer
		max_height = std::numeric_limits<float>::min();
		for (float x = a; x <= b; x += delta_x)
		{
			max_height = std::max(max_height, _pdf(x));
		}
		max_height += buffer;
	}

	//box dimensions
	float x_0 = a;
	float x_1 = b;
	float y_0 = 0;
	float y_1 = max_height;

	if (false)
	{
		float Area_box = (x_1 - x_0) * (y_1 - y_0);
		float Area_graph = 0;
		for (float x = a; x <= b; x += delta_x)
		{
			//reiman area approximation
			Area_graph += delta_x * _pdf(x);
		}
		float efficiency = Area_graph / Area_box;
		std::cout << "area box: " << Area_box << " area graph: " << Area_graph << " efficiency: " << efficiency << std::endl;
	}

	while (true)
	{
		float X = Helper::GetRandomNumber(x_0, x_1);
		float Y = Helper::GetRandomNumber(y_0, y_1);

		if (Y <= _pdf(X))
			return X;
	}
}

//given a bunch of the sampled data, it will generate a histogram to form a continuous function curve
inline std::vector<float> generateContinuousSampleHistogram(std::vector<float> sampledata, float a, float b, int bin_count,
	float scaling_factor = 1.0f, int rounding_method = 0)
{
	//generate samples and create histogram
	std::vector<float>index_samples;
	index_samples.resize(bin_count + 1);//add 1 because we need last index at N, not N-1 to get to b
	float delta_x_1 = (b - a) / (float)bin_count;
	for (int i = 0; i < sampledata.size(); i++)
	{
		float X_i = sampledata[i];
		int index = 0;
		switch (rounding_method)
		{
		case -1: index = std::floor((X_i - a) / (delta_x_1));
			break;
		case 0: index = std::round((X_i - a) / (delta_x_1));
			break;
		case 1: index = std::ceil((X_i - a) / (delta_x_1));
			break;
		default:index = std::round((X_i - a) / (delta_x_1));
		}
		if (index < 0 || index >= index_samples.size())
			continue;
		index_samples[index] += 1;
	}

	//normalize
	for (int i = 0; i < index_samples.size(); i++)
	{
		index_samples[i] = scaling_factor * (index_samples[i] / (float)sampledata.size());
	}

	std::vector<float> histogram_samples; //come in x,y pairs
	//ignore last one, it messes up with histogram
	for (int i = 0; i < index_samples.size() - 1; i++)
	{
		float x = a + (delta_x_1)*i;
		histogram_samples.push_back(x);
		histogram_samples.push_back(index_samples[i]);
	}

	return histogram_samples;
}

inline std::vector<float> generateContinuousRejectionTestSampleHistogram(std::function<float(float)> _pdf, float a, float b, int bin_count,
	                                                              int sample_count, float scaling_factor = 1.0f, int rounding_method = 0)
{
	//approximate max_height
	int N = 100;
	float buffer = 1.0f;
	float delta_x = (b - a) / (float)N;
	float max_height = std::numeric_limits<float>::min();
	for (float x = a; x <= b; x += delta_x)
	{
		max_height = std::max(max_height, _pdf(x));
	}
	max_height += buffer;

	//generate samples and create histogram
	std::vector<float>index_samples;
	index_samples.resize(bin_count + 1);//add 1 because we need last index at N, not N-1 to get to b
	float delta_x_1 = (b - a) / (float)bin_count;
	for (int i = 0; i < sample_count; i++)
	{
		float X_i = Continuous_Rejection_Sample(_pdf, a, b, max_height);
		int index = 0;
		switch (rounding_method)
		{
		case -1: index = std::floor((X_i - a) / (delta_x_1));
			break;
		case 0: index = std::round((X_i - a) / (delta_x_1));
			break;
		case 1: index = std::ceil((X_i - a) / (delta_x_1));
			break;
		default:index = std::round((X_i - a) / (delta_x_1));
		}

		index_samples[index] += 1;
	}

	//normalize
	for (int i = 0; i < index_samples.size(); i++)
	{
		index_samples[i] = scaling_factor * (index_samples[i] / (float)sample_count);
	}

	std::vector<float> histogram_samples; //come in x,y pairs
	//ignore last one, it messes up with histogram
	for (int i = 0; i < index_samples.size() - 1; i++)
	{
		float x = a + (delta_x_1)*i;
		histogram_samples.push_back(x);
		histogram_samples.push_back(index_samples[i]);
	}

	return histogram_samples;
}

//create an inversion_continuous sampler per function, itll calculate the CDF once with great precision and when you call sample it be dummy quick
//can even add a sample tester where it gives you an array of N sampled so you can easily plot them
class Continuous_Inversion_Sampler
{
public:
	Continuous_Inversion_Sampler(std::function<float(float)> _pdf, float _a, float _b, float precision_N) : pdf(_pdf), a(_a), b(_b), N(precision_N)
	{
		//compute CDF discretely with fundamental theorem of calculus and trapezoid integration, make sure first entry is 0 and last is 1
		cdf.resize(N + 1);
		float delta_x = (b - a) / (float)N;
		float sum = 0;
		cdf[0] = 0.0f;
		for (int n = 1; n < N + 1; n++)
		{
			//reiman method, just add the rectangle area: width*height = delta_x*pdf(x_n)
			float current_x = std::clamp(a + delta_x * n, a, b);
			//float current_x = a + delta_x * n;
			sum += delta_x * pdf(current_x);
			cdf[n] = sum;
		}
		//we have values from 0 to some X, since its an approximation it most likely overshot 1 maybe undershot 1. renormalize it, everything scales the same
		float scaling_term = 1.0f / cdf[N];
		for (int n = 1; n < N; n++)
		{
			cdf[n] *= scaling_term;
		}
		cdf[N] = 1.0f;
	}
	~Continuous_Inversion_Sampler() = default;

	float Sample(float U = 0, bool use_U = false) const
	{
		if(!use_U)
			U = Helper::GetRandomNumber(0, 1);

		//find index, n<U<n+1, binary search since sorted list, orders of magnitude faster
		int index = -1;
		int low = 0;
		int high = N;
		while (low <= high)
		{
			int mid = low + (high - low) / 2.0f;

			if (cdf[mid] < U && U <= cdf[mid + 1])
			{
				index = mid;
				break;
			}

			if (cdf[mid] < U)
				low = mid + 1;
			else
				high = mid - 1;
		};
		
		if (index == -1)
		{
			std::cout << "couldn't find index\n";
			return 0;
		}

		//now find the x, interpolating between the heights of index i and i + 1
		//Y= x_n + u*(x_n+1 - x_n)
		//u = e-f(index_n) / f(index_n+1) - f(index_n), if its a horizontal line then its not 1 to 1 but whatever itll give you x_n 
		float t = std::clamp((U - cdf[index]) / (cdf[index + 1] - cdf[index]), 0.0f, 1.f);
		float delta_x = (b - a) / (float)N;
		float X = (a + delta_x * index) + t * (delta_x * (index + 1) - delta_x * index);
		
		return X;
	}
	
	std::vector<float> generateTestSampleHistogram (int bin_count, int sample_count, float scaling_factor = 1.0f, int rounding_method = 0) const
	{
		std::vector<float>exponential_samples;
		exponential_samples.resize(bin_count + 1);//add 1 because we need last index at N, not N-1 to get to b
		float delta_x_1 = (b - a) / (float)bin_count;
		for (int i = 0; i < sample_count; i++)
		{
			float X_i = Sample();
			int index = 0;
			switch (rounding_method)
			{
			case -1: index = std::floor((X_i - a) / (delta_x_1));
				break;
			case 0: index = std::round((X_i - a) / (delta_x_1));
				break;
			case 1: index = std::ceil((X_i - a) / (delta_x_1));
				break;
			default:index = std::round((X_i - a) / (delta_x_1));
			}

			exponential_samples[index] += 1;
		}

		//normalize
		for (int i = 0; i < exponential_samples.size(); i++)
		{
			exponential_samples[i] = scaling_factor*(exponential_samples[i] / (float)sample_count);
		}

		std::vector<float> histogram_samples; //come in x,y pairs
		//ignore last one, it messes up with histogram
		for (int i = 0; i < exponential_samples.size() - 1; i++)
		{
			float x = a + (delta_x_1)*i;
			histogram_samples.push_back(x);
			histogram_samples.push_back(exponential_samples[i]);
		}

		return histogram_samples;
	}

	float getA() const { return a; }
	float getb() const { return b; }
	float getN() const { return N; }

	float PDF(float x) const { return pdf(x); }

private:
	int N;
	float a, b;
	std::function<float(float)> pdf;
	std::vector<float> cdf;
};

//this is for constant time inversion discrete, can implement later
class AliasTable
{
public:
	AliasTable() = default;
	AliasTable(std::vector<float> weights) : bins(weights.size())
	{
		//normalize the weights
		float sum = std::accumulate(weights.begin(), weights.end(), 0.);
		for (size_t i = 0; i < weights.size(); i++)
		{
			bins[i].p = weights[i] / sum;
		}
		
		//see if each bin is above or below the avg, 1/N
		struct Outcome {
			float scaled_prob;
			size_t index;
		};
		std::vector<Outcome> over, under;
		for (size_t i = 0; i < bins.size(); i++)
		{
			float scaled_prob = bins[i].p * bins.size();
			if (scaled_prob < 1)
				under.push_back(Outcome{ scaled_prob, i });
			else
				over.push_back(Outcome{ scaled_prob, i });
		}

		//probabilities are scaled by N, so probability of avg is 1


	}
	int Sample(float u, float* pmf, float* uRemapped = nullptr) const
	{
		//pick a random variable i
		//go to that bin and pick another random variable, and choose i or its alias depending on q_i, 1-q_i
	}
	std::string ToString() const;
	size_t size() const { return bins.size(); }
	float PMF(int index) const { return bins[index].p; }
private:
	struct Bin
	{
		float q, p; //p is probability from pmf given, q is the conditional probability, P(uniform=i and pmf=i)
		int alias;
	};
	std::vector<Bin> bins;
};











/*
//for inversion you can analytically get CDf, or compute it with trapezoid discrete, you can analytically inverse or use the interpolation method or root find
float inversion_continuous(std::function<float(float)> pdf, float a, float b, int precision_N)
{
	int N = precision_N;
	//compute CDF discretely with fundamental theorem of calculus and trapezoid integration, make sure first entry is 0 and last is 1
	std::vector<float> cdf(N+1);
	float delta_x = (b - a) / (float)N;
	if (a + delta_x * (float)N > b)
		std::cout << "inversion_continuous: the x goes over the b value with delta_x\n";
	float sum = 0;
	cdf[0] = 0.0f;
	for (int n = 1; n < N + 1; n++)
	{
		//reiman method, just add the rectangle area: width*height = delta_x*pdf(x_n)
		float current_x = a + delta_x * n;
		sum += delta_x * pdf(current_x);
		cdf[n] = sum;
		//std::cout << delta_x * pdf(current_x) << std::endl;
	}
	//we have values from 0 to some X, since its an approximation it most likely overshot 1 maybe undershot 1. renormalize it, everything scales the same
	float scaling_term = 1.0f / cdf[N];
	//std::cout << scaling_term << std::endl;
	for (int n = 1; n < N; n++)
	{
		cdf[n] *= scaling_term;
	}
	cdf[N] = 1.0f;
	//float current_x = a + delta_x * N;
	//std::cout<<"last sum: "<<delta_x * pdf(current_x)<<std::endl;
	//std::cout << 1.0f - sum << std::endl;
	//cdf[N] = 1.0f;

	//for (int i = 1; i < cdf.size(); i++)
	//{
		//std::cout << cdf[i]-cdf[i-1] << ", ";
	//}
	//std::cout << std::endl;

	//x=a + ((b-a)/N)*index
	//y=cdf(index)

	//generate uniform value between 0 and 1
	float U = Helper::GetRandomNumber(0, 1);

	//find index, n<U<n+1
	int index = -1;
	for (int i = 1; i < N + 1; i++)
	{
		if (cdf[i - 1] < U && U <= cdf[i])
		{
			index = i - 1;
			break;
		}
	}
	if (index == -1)
	{
		std::cout << "couldn't find index\n";
		return 0;
	}
	//std::cout <<"U: "<<U<<" index: "<< index << std::endl;
	//now find the x, interpolating between the heights of index i and i + 1
	//Y= x_n + u*(x_n+1 - x_n)
	//u = e-f(index_n) / f(index_n+1) - f(index_n), if its a horizontal line then its not 1 to 1 but whatever itll give you x_n
	float t = (U - cdf[index]) / (cdf[index + 1] - cdf[index]);
	if (t < 0 || t > 1)
	{
		std::cout << "inverse sampling t not 0 to 1\n";
	}
	float X = (a + delta_x*index) + t*(delta_x*(index+1) - delta_x*index);
	//std::cout << "t: " << t << "X: " << X << std::endl;
	float t_2 = t;
	float X_2 = a + (b - a) * (index + t_2) * (1.0f/(float)N);
	return X;
}
*/