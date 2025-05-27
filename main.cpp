#include "pch.h"

//#define POLYGON_APP
//#define MUSCLE_APP
//#define PLOTTER_TEST
//#define ODESOLVER_TEST

//RAY Tracing
//#define MONTECARLOS_TEST
#define RAYTRACER_TEST
//#define SAMPLER_TEST
//#define SPECTRUM_TEST
//#define SHAPE_TEST
//#define FILTER_TEST

#ifdef POLYGON_APP
#include "Applications/PolygonApp1.h"
#endif

#ifdef MUSCLE_APP
#include "Applications/MuscleCrossbridgeApp.h"
#endif

#ifdef PLOTTER_TEST
#include "Applications/PlotterTestApp.h"
#endif

#ifdef ODESOLVER_TEST
#include "Applications/ODESolverApp.h"
#endif

#ifdef MONTECARLOS_TEST
#include "Applications/RayTracingTests/MonteCarlosTestApp.h"
#endif

#ifdef RAYTRACER_TEST
#include "Applications/RayTracerTestApp.h"
#endif

#ifdef SAMPLER_TEST
#include "Applications/RayTracingTests/SamplerTestApp.h"
#endif

#ifdef SPECTRUM_TEST
#include "Applications/RayTracingTests/SpectrumColorTestApp.h"
#endif

#ifdef SHAPE_TEST
#include "Applications/RayTracingTests/ShapeTestApp.h"
#endif

#ifdef FILTER_TEST
#include "Applications/RayTracingTests/FilterFilmTestApp.h"
#endif

int main()
{
#ifdef POLYGON_APP
	PolygonApp1 app1;
	app1.Start();
#endif

#ifdef MUSCLE_APP
	MuscleCrossbridgeApp app2;
	app2.Start();
#endif

#ifdef PLOTTER_TEST
	PlotterTestApp app3;
	app3.Start();
#endif

#ifdef ODESOLVER_TEST
	ODESolverApp app4;
	app4.Start();
#endif

#ifdef MONTECARLOS_TEST
	MonteCarlosTestApp app5;
	app5.Start();
#endif

#ifdef RAYTRACER_TEST
	RayTracerTestApp app6;
	app6.Start();
#endif

#ifdef SAMPLER_TEST
	SamplerTestApp app7;
	app7.Start();
#endif

#ifdef SPECTRUM_TEST
	SpectrumColorTestApp app8;
	app8.Start();
#endif

#ifdef SHAPE_TEST
	ShapeTestApp app9;
	app9.Start();
#endif

#ifdef FILTER_TEST
	FilterFilmTestApp app10;
	app10.Start();
#endif

	return EXIT_SUCCESS;
}
