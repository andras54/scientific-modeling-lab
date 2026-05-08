#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>

#include "specialfunctions.h"

class DensityG
{
public:
	DensityG() = delete;
	DensityG(DensityG const&) = default;
	DensityG(DensityG&&) = default;
	DensityG& operator=(DensityG const&) = default;
	DensityG& operator=(DensityG&&) = default;
	DensityG(int _dim, int _sidelength, double _side, double _epsilon) : dim(_dim), sidelengthPowers(_dim, 1), alpha(std::sqrt(alglib::invincompletegammac((double)_dim * 0.5, _epsilon))),
		deltaRinv(_sidelength / _side), deltaR(_side / _sidelength), buffer(_dim), actualCube(_dim), sideLength(_sidelength), centralCube(_dim)
	{
		int length = 1;
		for (int i = 0; i < _dim; i++)
		{
			length *= _sidelength;
		}
		for (int i = _dim - 1; i > 0; i--)
		{
			sidelengthPowers[i - 1] = _sidelength * sidelengthPowers[i];
		}
		density = std::vector<double>(length);
		relativeCubeHolder = PermutationsAround0(_dim);
	}
	~DensityG() = default;

	void EvaluateDensity(std::vector<std::vector<double>>::iterator coordinates1, std::vector<std::vector<double>>::iterator const& coordinates2) // make sure that coordinates are inside sidelength!!!
	{
		std::fill(density.begin(), density.end(), 0.);
		for (; coordinates1 != coordinates2; ++coordinates1)
		{
			std::transform((*coordinates1).begin(), (*coordinates1).end(), centralCube.begin(), [&](double coord) {return (int)(coord * deltaRinv); });
			for (int i = 0; i < (int)relativeCubeHolder.size(); i = i + dim)
			{
				std::transform(centralCube.begin(), centralCube.end(), relativeCubeHolder.begin() + i, actualCube.begin(), [](int cent, int rel) {return cent + rel; });
				if (IsValid(actualCube))
				{
					density[EvaluateCoordinate(actualCube)] += Fraction(actualCube, *coordinates1);
				}
			}
		}
		//double averageDensityInv = density.size() / std::accumulate(density.begin(), density.end(), 0.);
		//std::transform(density.begin(), density.end(), density.begin(), [&](double d) {return d * averageDensityInv - 1.; });
	}

	std::vector<double> const& GetDensity()
	{
		return density;
	}

private:
	bool IsNeighbour(std::vector<int> const& cube, std::vector<double> const& coordinates)
	{
		std::transform(coordinates.begin(), coordinates.end(), cube.begin(), buffer.begin(), [&](double co, int c) {return co - std::max(c * deltaR, std::min((c + 1) * deltaR, co)); });
		return std::sqrt(std::inner_product(buffer.begin(), buffer.end(), buffer.begin(), 0.)) < deltaR;
	}

	bool IsValid(std::vector<int> const& cube)
	{
		for (int i : cube)
		{
			if (i < 0 || i >= sideLength)
			{
				return false;
			}
		}
		return true;
	}

	double Fraction(std::vector<int> const& cube, std::vector<double> const& coordinates)
	{
		if (IsNeighbour(cube, coordinates))
		{
			std::transform(coordinates.begin(), coordinates.end(), cube.begin(), buffer.begin(), [&](double co, int c) {return alpha * (c - co * deltaRinv); });
			return std::accumulate(buffer.begin(), buffer.end(), 1., [&](double b1, double b2) {return b1 * (alglib::errorfunction(b2 + alpha) - alglib::errorfunction(b2)); });
		}
		else
		{
			return 0.;
		}
	}

	int EvaluateCoordinate(std::vector<int> const& cube)
	{
		return std::inner_product(cube.begin(), cube.end(), sidelengthPowers.begin(), 0);
	}

	std::vector<int> PermutationsAround0(int dim)
	{
		auto Step = [](int step) {if (step == 1) { return 0; } else if (step == 0) { return -1; } else { return 1; }};

		int length = 1;
		for (int i = 0; i < dim; i++)
		{
			length *= 3;
		}
		std::vector<int> result(length * dim);
		int step = 1;
		int lengthr = length;
		for (int i = 0; i < dim; i++)
		{
			lengthr /= 3;
			for (int j = 0; j < length; j++)
			{
				if ((j - 1) / lengthr != j / lengthr)
				{
					step = Step(step);
				}
				result[dim * j + i] = step;
			}
			step = Step(step);
		}
		return result;
	}

	std::vector<double> density;
	std::vector<int> sidelengthPowers;
	const double deltaRinv;
	const double deltaR;
	const int dim;
	const double alpha;
	const int sideLength;

	std::vector<double> buffer;
	std::vector<int> centralCube;
	std::vector<int> actualCube;
	std::vector<int> relativeCubeHolder;
};