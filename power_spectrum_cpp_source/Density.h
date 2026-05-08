#include <vector>
#include <cmath>
#include <numeric>

#include "Divide.h"

class Density
{
public:
	Density() = delete;
	Density(Density const&) = default;
	Density(Density&&) = default;
	Density& operator=(Density const&) = default;
	Density& operator=(Density&&) = default;
	Density(int _dim, int _sidelength, double _side) : div(_dim, _sidelength), t(_dim), r(_dim), dim(_dim)
	{
		deltaRinv = _sidelength / _side;
		//deltaRinv3 = deltaRinv * deltaRinv * deltaRinv;
		int length = 1;
		for (int i = 0; i < _dim; i++)
		{
			length *= _sidelength;
		}
		density = std::vector<double>(length);
		sidelengthPowers = div.GetSidelengthPowers();
	}
	~Density() = default;

	void EvaluateDensity(std::vector<std::vector<double>>::iterator coordinates1, std::vector<std::vector<double>>::iterator const& coordinates2)//(std::vector<std::vector<double>> const& coordinates) // make sure that coordinates are inside sidelength!!!
	{
		std::fill(density.begin(), density.end(), 0.);
		for (; coordinates1 != coordinates2; ++coordinates1)//(std::vector<double> const& c : coordinates)
		{
			SetTR(*coordinates1);
			div.EvaluateDivisions(t, r);
			for (int i = 0; i < (int)div.GetCoordinates().size(); i++)
			{
				if (div.GetValid()[i])
				{
					density[div.GetCoordinates()[i]] += div.GetFractions()[i];
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
	std::vector<double> density;
	std::vector<int> t;
	std::vector<double> r;
	std::vector<int> sidelengthPowers;
	double deltaRinv;
	//double deltaRinv3;
	Divide div;
	const int dim;
	double remainder = 0;

	void SetTR(std::vector<double> const& coordinates)
	{
		for (int i = 0; i < dim; i++)
		{
			remainder = coordinates[i] * deltaRinv;
			t[i] = (int)std::round(remainder) - 1;
			remainder -= (int)remainder;
			if (remainder < 0.5)
			{
				r[i] = 0.5 - remainder;
			}
			else
			{
				r[i] = 1.5 - remainder;
			}
		}
	}
};