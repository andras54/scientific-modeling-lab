#include <vector>
#include <algorithm>

class Divide
{
public:
	Divide() = delete;
	Divide(Divide const&) = default;
	Divide(Divide&&) = default;
	Divide& operator=(Divide const&) = default;
	Divide& operator=(Divide &&) = default;
	Divide(int _dim, int _sidelength) : sidelengthPowers(_dim, 1), dim(_dim), sidelength(_sidelength)
	{
		int length = 2;
		for (int i = _dim - 1; i > 0; i--)
		{
			sidelengthPowers[i - 1] = _sidelength * sidelengthPowers[i];
			length *= 2;
		}
		coordinates = std::vector<int>(length);
		fractions = std::vector<double>(length);
		valid = std::vector<bool>(length);
	}
	~Divide() = default;
	
	void EvaluateDivisions(std::vector<int> const& t, std::vector<double> const& r)
	{
		std::fill(coordinates.begin(), coordinates.end(), 0);
		std::fill(fractions.begin(), fractions.end(), 1.);
		std::fill(valid.begin(), valid.end(), true);
		step = 1;
		for (int i = 0; i < dim; i++)
		{
			std::transform(coordinates.begin(), coordinates.begin() + step, coordinates.begin(), [&](int c0) {return c0 + t[i] * sidelengthPowers[i]; });
			std::transform(coordinates.begin(), coordinates.begin() + step, coordinates.begin() + step, [&](int c1) {return c1 + sidelengthPowers[i]; });

			std::transform(fractions.begin(), fractions.begin() + step, fractions.begin() + step, [&](double f1) {return f1 * (1. - r[i]); });
			std::transform(fractions.begin(), fractions.begin() + step, fractions.begin(), [&](double f0) {return f0 * r[i]; });

			std::transform(valid.begin(), valid.begin() + step, valid.begin() + step, [&](bool v1) {return v1 && (t[i] + 1 <= sidelength - 1); });
			std::transform(valid.begin(), valid.begin() + step, valid.begin(), [&](bool v0) {return v0 && (t[i] >= 0); });

			step *= 2;
		}
	}

	std::vector<int> const& GetCoordinates() const
	{
		return coordinates;
	}

	std::vector<int> const& GetSidelengthPowers() const
	{
		return sidelengthPowers;
	}

	std::vector<double> const& GetFractions() const
	{
		return fractions;
	}

	std::vector<bool> const& GetValid() const
	{
		return valid;
	}

private:
	std::vector<int> coordinates;
	std::vector<double> fractions;
	std::vector<bool> valid;
	std::vector<int> sidelengthPowers;
	int step = 1;
	const int dim;
	const int sidelength;
};