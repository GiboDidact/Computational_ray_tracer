#pragma once
#include "../pch.h"
#include "Sampling.h"
/*
if you have some integral(f(x)), then you can choose an estimator to approximate f(x)
luckily for monte carlos method we are trying to estimate a mean which is the simplest one, so we use mean estimators
remember E[f(x)]=integral(f(x)p(x)dx), the E[f(x)]=u, u is the mean of the random variable f(x). But we know by 
law of large numbers that, F being an avg of i.i.d random variable f(x), E[F]=u in the limit of samples-->infinity
So the limit of the avg estimator F = u = E[f(x)] = integral(f(x)p(x)dx)

we have the standard, unbiased, monte carlos estimator that you can choose any pdf(x) on to sample from.
with uniform: (b-a)/n*sumf(x_i)
with nonuniform: 1/n*sumf(x_i)/pdf(x_i)
or you could potentially create your own biased estimator that tries to approximate the expected value of f(x)

for the standard unbiased monte carlo estimator: the expected value is always equal to the integral, the variance changes for different pdf(x) thats the 
whole point. changing pdf(x) still always keeps it consistent and unbiased, however variance changes. the whole problem is finding ones with lowest variance

then based on the estimator you have you can study some of its theoritcal analysis:
. the bias, or expected value - desired value
. if its consistent or not
. the variance
. the efficiency: 1 / variance * running_time
. the mean squared error: variance if unbiased
. even the error bounds with chebyshev inequality, central limit theorem

*though most of these are hard to figure out, or you have to do it analytically, you might be able to numerically figure some out though

then you can extent to multiple dimensions, and its pretty much easy to extent its the same thing its analogous in n dimensions


variance:
the main thing we want to do is try different pdfs(x) and see which has the least variance. But how do we calculate variance?
Mathematically its 1/N*variance(g/pdf) or E[X^2] - E[X]^2. However from here theres no general reduction and if you had a specific case
you could try to analytically solve it though E[X^2] term can be hard as well as variance(g/pdf). It probably is hard to ever solve these
analytically, so how to approximate it? Well from statistics the answer can be from estimators. It seems reasonable that we can approximate it with a bunch
of samples from the distribution by some estimator, and it will get better with more samples. 

*/

//really the only thing I need is a way to sample generally given distribution functions
 
//**create different estimators and compare variance, variacne should change with different distribution samples and how many samples N
class distributionEstimator
{
public:
	distributionEstimator() = default;
	~distributionEstimator() = default;

	//integral_function and pdf need same range of a to b. pdf should be normalized
	//**in order for this to work pdf has to: pdf(x)>=0, and integral(pdf(x)) = 1, or else this isn't mathematically correct
	float EstimateFunctionN(float N, std::function<float(float)> integral_function, std::function<float(float)> pdf, float a, float b, 
		                    float sampling_precision = 500)
	{
		float sum = 0;
		Continuous_Inversion_Sampler inversion_sampler(pdf, a, b, sampling_precision);
		for (int i = 0; i < N; i++)
		{
			float X_i = inversion_sampler.Sample();
			float g_x = integral_function(X_i);
			float p_x = pdf(X_i);
			sum += g_x / p_x; //p_x can't be zero
		}

		return  sum / N;
	}

	float bias() { return 0; }
	bool consistent() { return true; }

	float approximate_expected_value(std::function<float(float)> integral_function, std::function<float(float)> pdf, float a, float b)
	{
		return EstimateFunctionN(approximate_samples, integral_function, pdf, a, b);
	}

	//*pdf has to be valid, always positive and integrate to 1
	float approximate_variance(std::function<float(float)> integral_function, std::function<float(float)> pdf, float a, float b)
	{
		//generate all samples and calculate the g(x)/p(x) estimator mean
		Continuous_Inversion_Sampler inversion_sampler(pdf, a, b, 2000);
		float sample_variance_size = approximate_samples;
		std::vector<float> sample_points(sample_variance_size);
		float sum = 0;
		for (int i = 0; i < sample_variance_size; i++)
		{
			float X_i = inversion_sampler.Sample();
			float X_transform =  integral_function(X_i) / pdf(X_i);

			sum += X_transform;
			sample_points[i] = X_transform;
		}
		float EV = sum / sample_variance_size;
		
		//calculate the g(x)/p(x) estimator variance
		sum = 0;
		for (int i = 0; i < sample_variance_size; i++)
		{
			float X_i = sample_points[i];
			sum += std::powf(X_i - EV, 2.0f);
		}
		float VAR_GF = sum / sample_variance_size;
		
		//don't forget for variance of monte carlo estimator to divide by N. Var(monte estimator) = (1/N)*variance(g(x)/f(x))
		return VAR_GF / approximate_samples;
	}

	float approximate_efficiency(std::function<float(float)> integral_function, std::function<float(float)> pdf, float a, float b)
	{
		float Var = approximate_variance(integral_function, pdf, a, b);
		Timer time1;
		time1.Begin();
		EstimateFunctionN(approximate_samples, integral_function, pdf, a, b);
		return 1.0f / (Var * time1.getTimeNano());
	}

private:
	int approximate_samples = 100000;
};

class uniformEstimator
{
public:
	float EstimateFunctionN(float N, std::function<float(float)> integral_function, float a, float b)
	{
		float sum = 0;
		for (int i = 0; i < N; i++)
		{
			float X_i = Helper::GetRandomNumber(a, b);
			sum += integral_function(X_i);
		}

		return ((b - a) / N) * sum;
	}

	float bias() { return 0; }
	bool consistent() { return true; }
	float theoretical_expected_value()
	{
		//E[estimator] = integral(f(x))
	}
	float theoretical_variance()
	{
		//V[estimator] = (b-a)^2/N*V[f(x)], x ~ uniform(a,b)
	}
	float approximate_expected_value(std::function<float(float)> integral_function, float a, float b)
	{
		return EstimateFunctionN(approximate_samples, integral_function, a, b);
	}
	//variance... I think for variance its relative to the same function and domain, if you have the same thing your trying to calculate with the
	//estimator then you can compare the variances, given the same function and domain. Or else there just different and you can't compare.
	//also variance changes with sample count too 1/N
	//Variance[estimator] = 1/N*Variance(f(x)/p(x))
	float approximate_variance(std::function<float(float)> integral_function, float a, float b)
	{
		float sample_sum = 0;
		std::vector<float> sample_points;
		float sample_variance_size = approximate_samples;
		sample_points.resize(sample_variance_size);
		for (int i = 0; i < sample_variance_size; i++)
		{
			float X_i = Helper::GetRandomNumber(a, b);

			sample_sum += integral_function(X_i);
			sample_points[i] = integral_function(X_i);
		}

		float E_V = (1.0f / sample_variance_size) * sample_sum;
		float sum = 0;
		for (int i = 0; i < sample_variance_size; i++)
		{
			float X_i = sample_points[i];
			sum += std::powf(X_i - E_V, 2);
		}

		//this is variance[f(x)] (the p(x) of uniform is 1/b-a so you take it out and the variance squares it)
		float variance_f = (1.0f / (sample_variance_size - 1.0f)) * sum;
		//for variance of actual estimator: divide by N again and add the (b-a)^2 on top
		return (std::powf(b-a,2)/ approximate_samples) * variance_f;
	}

	//efficiency = 1/variance(estimator)*time_to_run(estimator)
	float approximate_efficiency(std::function<float(float)> integral_function, float a, float b)
	{
		float Var = approximate_variance(integral_function, a, b);
		Timer time1;
		time1.Begin();
		EstimateFunctionN(approximate_samples, integral_function, a, b);
		return 1.0f / (Var * time1.getTimeNano());
	}
	void PrintChebychevInequality(float delta, std::function<float(float)> integral_function, float a, float b)
	{
		//float G = EstimateFunctionN(1, 0);
		//float E_G = approximate_expected_value();
		float V_G = approximate_variance(integral_function, a, b);

		//float error = std::abs(G - E_G);
		float floor = std::sqrtf(V_G / delta);
		std::cout << "P[" << "|G - E[G]|" << ">=" << floor << "] <= " << delta << std::endl;
	}
	void PrintChebychevInequalityRanges(std::function<float(float)> integral_function, float a, float b)
	{
		float V_G = approximate_variance(integral_function, a, b);

		for (float delta = 0.05; delta <= 1.0f; delta += .05)
		{
			float floor = std::sqrtf(V_G / delta);
			std::cout << "P[" << "|G - E[G]|" << ">=" << floor << "] <= " << delta << std::endl;
		}
	}
	//technically for more error confidence intervals you could use central limit theorem of the 68/95/99.8 rule

private:
	int approximate_samples = 100000;
};
