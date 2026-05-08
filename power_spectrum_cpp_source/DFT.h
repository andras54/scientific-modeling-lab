#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>

typedef std::complex<double> comp;

const comp euler(2.71828182845904523536, 0.);
const comp pi(3.14159265358979323846264, 0.);
const comp im(0., 1.);

class DFT
{
public:
	DFT() = delete;
	DFT(DFT const&) = default;
	DFT(DFT&&) = default;
	DFT& operator=(DFT const&) = default;
	DFT& operator=(DFT&&) = default;

	DFT(int _dim, int _sidelength) : dim(_dim), sidelength(_sidelength), length(1), sidelengthPowers(_dim)
	{
		if (dim < 1)
		{
			std::cout << "DFT object with dimension less than 1 can not be constructed!" << std::endl;
			exit(-1);
		}
		for (int i = 0; i < _dim; i++)
		{
			sidelengthPowers[dim - 1 - i] = length;
			length *= sidelength;
		}
		E = std::vector<comp>(sidelength * sidelength, comp(1., 0.));
		comp exp = std::pow(euler, -im * 2. * pi / (double)sidelength);
		for (int i = 1; i < sidelength; i++)
		{
			E[sidelength + i] = E[sidelength + i - 1] * exp;
		}
		for (int i = 2; i < sidelength; i++)
		{
			std::transform(E.begin() + (i - 1) * sidelength, E.begin() + i * sidelength, E.begin() + sidelength, E.begin() + i * sidelength, [&](comp x1, comp x2) {return x1 * x2;});
		}
	}

	void CheckE()
	{
		for (int i = 0; i < sidelength;i++)
		{
			for (int j = 0; j < sidelength; j++)
			{
				std::cout << E[i * sidelength + j] << " ";
			}
			std::cout << std::endl;
		}
	}

	std::vector<comp> FourierTransform(std::vector<comp> const& delta)
	{
		which = false;

		if ((int)delta.size() != length)
		{
			std::cout << "Not compatible size!" << std::endl;
		}

		std::vector<comp> buffer0(delta);
		std::vector<comp> buffer1(length);

		int eRow;
		int indexbaseForA;
		int indexstepForA;
		int actualJ;

		for (int i = 0; i < dim; i++)
		{
			if (!which)
			{
				for (int j = 0; j < length; j++)
				{
					actualJ = j;
					for (int k = 0; k < dim - i - 1; k++)
					{
						actualJ /= sidelength;
					}
					eRow = actualJ % sidelength;
					indexbaseForA = j - eRow * sidelengthPowers[i];

					eRow *= sidelength;
					indexstepForA = 0;

					for (int k = 0; k < sidelength; k++)
					{
						buffer1[j] += E[eRow + k] * buffer0[indexbaseForA + indexstepForA];
						indexstepForA += sidelengthPowers[i];
					}
				}
				std::fill(buffer0.begin(), buffer0.end(), comp(0., 0.));
			}
			else
			{
				for (int j = 0; j < length; j++)
				{
					actualJ = j;
					for (int k = 0; k < dim - i - 1; k++)
					{
						actualJ /= sidelength;
					}
					eRow = actualJ % sidelength;
					indexbaseForA = j - eRow * sidelengthPowers[i];

					eRow *= sidelength;
					indexstepForA = 0;

					for (int k = 0; k < sidelength; k++)
					{
						buffer0[j] += E[eRow + k] * buffer1[indexbaseForA + indexstepForA];
						indexstepForA += sidelengthPowers[i];
					}
				}
				std::fill(buffer1.begin(), buffer1.end(), comp(0., 0.));
			}
			which = !which;
		}

		if (!which)
		{
			return buffer0;
		}
		else
		{
			return buffer1;
		}
	}
private:
	int length;
	int sidelength;
	int dim;

	bool which = false;

	std::vector<comp> E;
	std::vector<int> sidelengthPowers;
};