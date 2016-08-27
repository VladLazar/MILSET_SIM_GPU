#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/count.h>
#include <thrust/transform.h>
#include <thrust/random.h>
#include <iostream>
#include <time.h>
#include <utility>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <thrust/tuple.h>

int N = 1000;


struct Point{
	long double x, y, z;
};

struct randomPoint{ //generez puncte random intre 0 si 1.0
	__host__ __device__ thrust::tuple<double, double, double, double, double, double> operator()(const unsigned int n){
		thrust::default_random_engine rng; // fac generatoru
		rng.discard(n); //???
		return thrust::make_tuple(
			(double)rng() / thrust::default_random_engine::max,//x
			(double)rng() / thrust::default_random_engine::max,//y
			(double)rng() / thrust::default_random_engine::max,//aux pt unghi
			(double)rng() / thrust::default_random_engine::max,//alt aux
			(double)rng() / (thrust::default_random_engine::max / 0.87) + 0.174,//unghi mic de bounce
			(double)rng() / (thrust::default_random_engine::max / 2.094) + 1.047);//unghi mare de bounce
	}
};


struct functor{

	const float3 nucleus = make_float3(0.5, 0.5, 100);
	const int distAtomCerc = 50;
	const double razaAtomuluiAur = 0.7;
	const long long distAtomPlan = 100000;
	const int Zp = 2;
	const int ZAu = 79;//for the nucleus
	const long double e = 1.602;
	const int v = 15; // m/s v standard = 15000000 m/s; val * 10 ^ 6
	const long double epsilon0 = 8.85;
	const long double pi = 3.14;//add more
	const long double alphaMass = 6.644; //in kg
	const long double eps = 0.0000000000001;

	__host__ __device__ thrust::tuple<double, double, double, double> operator()(thrust::tuple<double, double, double, double, double, double> punct){
		long double x = thrust::get<0>(punct);
		long double y = thrust::get<1>(punct);
		long double z = 0;
		long double w = thrust::get<2>(punct);

		long double b = sqrt((x - nucleus.x) * (x - nucleus.x) + (y - nucleus.y) * (y - nucleus.y));
		//x >= 0.45 && x <= 0.55 && y >= 0.45 && y <= 0.55
		if (b < 0.1)
		{
			if (w < 0.7)
			{
				long double kk = thrust::get<4>(punct);

				Point A, B;
				double radiusOfCircle = distAtomCerc * tan(kk);
				long double angle = thrust::get<3>(punct) * 2 * pi;

				A.z = nucleus.z;
				A.x = x;
				A.y = y;
				B.z = nucleus.z + distAtomCerc;
				//B.x = x + cos(angle * pi / 180.0);
				//B.y = y + sin(angle * pi / 180.0);
				B.x = x + cos(angle) * radiusOfCircle;
				B.x = y + sin(angle) * radiusOfCircle;

				long double t = ((nucleus.z + distAtomPlan) - A.z) / (B.z - A.z); //parameter;

				return thrust::make_tuple(
					(double)A.x + (B.x - A.x) * t,
					(double)A.y + (B.y - A.y) * t,
					(double)nucleus.z + distAtomPlan,
					(double)(kk * 180) / pi);
			}
			else
			{
				long double kk = thrust::get<5>(punct);
				return thrust::make_tuple(0, 0, 0, (kk * 180) / pi);
			}
		}

		//long double tg2 = (Zp * ZAu * e * e / (b * 2 * pi * epsilon0 * alphaMass * v * v)) * pow(10.0, 4);//tg (unghi de dev / 2)
		double tg2 = 0.0004993 / b;
		double tg = 2 * tg2 / (1 - tg2 * tg2);//tg deviation angle
		double angle = w * 2 * pi;

		//if (w > 0.5) tg = fabs(tg);

		double radiusOfCircle = distAtomCerc * tg;

		if (distAtomCerc <= distAtomPlan)
		{
			Point A, B;

			A.z = nucleus.z;
			A.x = x;
			A.y = y;
			B.z = nucleus.z + distAtomCerc;
			B.x = x + cos(angle) * radiusOfCircle;
			B.y = y + sin(angle) * radiusOfCircle;

			double t = ((nucleus.z + distAtomPlan) - A.z) / (B.z - A.z); //parameter;

			return thrust::make_tuple(
				(double)A.x + (B.x - A.x) * t,
				(double)A.y + (B.y - A.y) * t,
				(double)nucleus.z + distAtomPlan,
				(double)atan(tg) * 180 / pi);
		}
		else
		{
			return thrust::make_tuple(69, 69, 69, 69);//the reflected atom does not hit the screen
		}
	}
};

int main()
{
	//std::cin >> N;

	thrust::device_vector< thrust::tuple<double, double, double, double, double, double> > particles(N);
	thrust::device_vector< thrust::tuple<double, double, double, double> > intersectPlan(N);
	thrust::counting_iterator<unsigned int> counter(0);
	thrust::transform(counter, counter + N, particles.begin(), randomPoint());

	thrust::transform(particles.begin(), particles.end(), intersectPlan.begin(), functor());

	//std::cout << std::setprecision(3) << std::fixed;

	for (int i = 0; i < intersectPlan.size(); ++i)
	{
		thrust::tuple<double, double, double, double> aux = intersectPlan[i];
		std::cout << "NUM: " << i + 1 << " " << thrust::get<0>(aux) << " " << thrust::get<1>(aux) << " " << thrust::get<2>(aux) << " " << thrust::get<3>(aux) << "\n";
	}

	std::cout << "\n\n-----> DATA FOR GRAPH <-----\n\n";

	for (int i = 0; i < intersectPlan.size(); ++i)
	{
		thrust::tuple<double, double, double, double> aux = intersectPlan[i];
		std::cout << std::setiosflags(std::ios::fixed) << thrust::get<0>(aux) << ',' << thrust::get<1>(aux) << "\n";
	}

	std::cout << "\n\n-----> DATA FOR GRAPH2 <-----\n\n";

	for (int i = 0; i < particles.size(); ++i)
	{
		thrust::tuple<double, double, double, double, double, double> aux = particles[i];
		std::cout << std::setiosflags(std::ios::fixed) << thrust::get<0>(aux) << ',' << thrust::get<1>(aux) << "\n";
	}

	return 0;
}
