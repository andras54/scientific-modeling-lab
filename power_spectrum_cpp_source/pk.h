#include <tuple>
#include <stdexcept>

void SaveVector(std::string name, std::vector<double> v)
{
	std::ofstream save(name + ".dat");
	for (double x : v)
	{
		save << std::setprecision(16) << x << std::endl;
	}
}

void SaveVectors(std::string name, std::vector<double> v1, std::vector<double> v2)
{
	std::ofstream save(name);
	for (int i = 0; i<v1.size(); i++)
	{
		save << std::setprecision(16) << v1[i] << " " << v2[i] << std::endl;
	}
}

int addbool(int a, bool b) {
	return a + (int)b;
}

template<typename t>
std::vector<t> vflip(std::vector<t> const& v)
{
	std::vector<t> result(v.size());
	int k = v.size() - 1;
	for (int i = 0; i < v.size();i++)
	{
		result[i] = v[k];
		--k;
	}
	return result;
}

template<typename t>
std::vector<t> boolidx(std::vector<t> const& v, std::vector<bool> const& b)
{
	std::vector<t> result(std::accumulate(b.begin(), b.end(), 0, addbool));
	int k = 0;
	for (int i = 0; i < v.size();i++)
	{
		if (b[i])
		{
			result[k] = v[i];
			k += 1;
		}
	}
	return result;
}

template<typename t>
std::vector<bool> operator==(std::vector<t> const& v1, std::vector<t> const& v2)
{
	if (v1.size() != v2.size())
	{
		throw std::runtime_error("Vector sizes don't match!");
	}
	std::vector<bool> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] == v2[i];
	}
	return result;
}

template<typename t>
std::vector<bool> operator<=(std::vector<t> const& v1, std::vector<t> const& v2)
{
	if (v1.size() != v2.size())
	{
		throw std::runtime_error("Vector sizes don't match!");
	}
	std::vector<bool> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] <= v2[i];
	}
	return result;
}

template<typename t>
std::vector<bool> operator>=(std::vector<t> const& v1, std::vector<t> const& v2)
{
	if (v1.size() != v2.size())
	{
		throw std::runtime_error("Vector sizes don't match!");
	}
	std::vector<bool> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] >= v2[i];
	}
	return result;
}

template<typename t>
std::vector<bool> operator<(std::vector<t> const& v1, std::vector<t> const& v2)
{
	if (v1.size() != v2.size())
	{
		throw std::runtime_error("Vector sizes don't match!");
	}
	std::vector<bool> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] < v2[i];
	}
	return result;
}

template<typename t>
std::vector<bool> operator>(std::vector<t> const& v1, std::vector<t> const& v2)
{
	if (v1.size() != v2.size())
	{
		throw std::runtime_error("Vector sizes don't match!");
	}
	std::vector<bool> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] > v2[i];
	}
	return result;
}

std::vector<bool> operator&&(std::vector<bool> const& v1, std::vector<bool> const& v2)
{
	if (v1.size() != v2.size())
	{
		throw std::runtime_error("Vector sizes don't match!");
	}
	std::vector<bool> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] && v2[i];
	}
	return result;
}

template<typename t>
std::vector<t> operator-(std::vector<t> const& v, t num)
{
	std::vector<t> result(v.size());
	for (int i = 0; i < v.size(); i++)
	{
		result[i] = v[i] - num;
	}
	return result;
}

template<typename t>
std::vector<t> operator*(std::vector<t> const& v1, std::vector<t> const& v2)
{
	std::vector<t> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] * v2[i];
	}
	return result;
}

std::vector<double> operator/(std::vector<double> const& v1, std::vector<double> const& v2)
{
	std::vector<double> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] / v2[i];
	}
	return result;
}

template<typename t>
std::vector<bool> operator==(std::vector<t> const& v, t num)
{
	std::vector<bool> result(v.size());
	for (int i = 0; i < v.size(); i++)
	{
		result[i] = v[i] == num;
	}
	return result;
}

template<typename t>
std::vector<bool> operator<(std::vector<t> const& v, t num)
{
	std::vector<bool> result(v.size());
	for (int i = 0; i < v.size(); i++)
	{
		result[i] = v[i] < num;
	}
	return result;
}

template<typename t>
std::vector<bool> operator>(std::vector<t> const& v, t num)
{
	std::vector<bool> result(v.size());
	for (int i = 0; i < v.size(); i++)
	{
		result[i] = v[i] > num;
	}
	return result;
}

template<typename t>
std::vector<bool> operator>=(std::vector<t> const& v, t num)
{
	std::vector<bool> result(v.size());
	for (int i = 0; i < v.size(); i++)
	{
		result[i] = v[i] >= num;
	}
	return result;
}

template<typename t>
std::vector<bool> operator<=(std::vector<t> const& v, t num)
{
	std::vector<bool> result(v.size());
	for (int i = 0; i < v.size(); i++)
	{
		result[i] = v[i] <= num;
	}
	return result;
}

namespace std
{
	template<class InputIt1, class InputIt2, class InputIt3, // similar to https://en.cppreference.com/w/cpp/algorithm/transform.html "possible implementation"
		class OutputIt, class TrinaryOp>
	constexpr
		OutputIt transform3(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt3 first3,
			OutputIt d_first, TrinaryOp trinary_op)
	{
		for (; first1 != last1; ++d_first, ++first1, ++first2, ++first3)
			*d_first = trinary_op(*first1, *first2, *first3);
		return d_first;
	}
}

struct idx
{
	idx(int _sidelength, int _period) : sidelength(_sidelength), period(_period) {}
	int sidelength;
	int period;
	int counter = 0;
	int idx_counter = sidelength - 1;

	int Step()
	{
		if (counter % period == 0)
		{
			next();
		}

		counter += 1;

		return idx_counter;
	}

	void next()
	{
		if (idx_counter < sidelength - 1)
		{
			idx_counter += 1;
		}
		else
		{
			idx_counter = 0;
		}
	}
};

template<typename t>
std::vector<int> searchsorted(std::vector<t> const& a, std::vector<t> const& v, std::string const& side)
{
	for (int i = 0; i < a.size() - 1; i++)
	{
		if (a[i] > a[i + 1])
		{
			throw std::invalid_argument("Not ascending!");
		}
	}
	std::vector<int> result(v.size());
	for (int i = 0; i < v.size();i++)
	{
		if (side == "left")
		{
			if (a[a.size() - 1] < v[i])
			{
				result[i] = a.size();
			}
			else if (v[i] <= a[0])
			{
				result[i] = 0;
			}
			else
			{
				for (int j = 1; j < a.size();j++)
				{
					if (a[j - 1] < v[i] && v[i] <= a[j])
					{
						result[i] = j;
						break;
					}
				}
			}
		}
		else
		{
			if (a[a.size() - 1] <= v[i])
			{
				result[i] = a.size();
			}
			else if (v[i] < a[0])
			{
				result[i] = 0;
			}
			else
			{
				for (int j = a.size() - 1; j > 0; j--)
				{
					if (a[j - 1] <= v[i] && v[i] < a[j])
					{
						result[i] = j;
						break;
					}
				}
			}
		}
	}
	return result;
}

template<typename t>
std::vector<int> digitize(std::vector<t> const& x, std::vector<t> const& bins)
{
	return searchsorted(bins, x, "right");
}


std::vector<double> fftfreq(int n, double d)
{
	std::vector<double> result(n);
	double factor = 1. / (d * n);
	if (n % 2 == 0)
	{
		for (int i = 0; i < n / 2; i++)
		{
			result[i] = (double)i;
		}
		for (int i = n / 2; i < n; i++)
		{
			result[i] = (double)i - n;
		}
	}
	else
	{
		for (int i = 0; i < (n + 1) / 2; i++)
		{
			result[i] = (double)i;
		}
		for (int i = (n + 1) / 2; i < n; i++)
		{
			result[i] = (double)i - n;
		}
	}
	std::transform(result.begin(), result.end(), result.begin(), [&](double r) {return factor * r;});
	return result;
}

std::vector<double> fourier_grid_box_3d(int sidelength, double side, int dim)
{
	double dcell = side / sidelength;
	std::vector<double> k(fftfreq(sidelength, dcell));
	std::transform(k.begin(), k.end(), k.begin(), [&](double k) {return 2 * pi.real() * k;});
	//int s3 = sidelength * sidelength * sidelength;
	int sn = sidelength;
	std::vector<int> sidelengthPowers(dim,1);
	for (int i = dim - 2; i > -1;i--)
	{
		sidelengthPowers[i] = sidelengthPowers[i + 1] * sidelength;
		sn *= sidelength;
	}

	std::vector<double> kvec(sn * dim);

	for (int j = 0; j < dim; j++)
	{
		idx Idx(sidelength, sidelengthPowers[j]);
		for (int i = j * sn; i < (j+1) * sn; i++)
		{
			kvec[i] = k[Idx.Step()];
		}
	}
	//idx Idx1(sidelength, sidelength * sidelength);
	//idx Idx2(sidelength, sidelength);
	//idx Idx3(sidelength, 1);

	//for (int i = 0; i < s3; i++)
	//{
	//	kvec[i] = k[Idx1.Step()];
	//}
	//for (int i = s3; i < 2 * s3; i++)
	//{
	//	kvec[i] = k[Idx2.Step()];
	//}
	//for (int i = 2 * s3; i < 3 * s3; i++)
	//{
	//	kvec[i] = k[Idx3.Step()];
	//}
	std::vector<double> result(sn);
	//std::transform3(kvec.begin(), kvec.begin() + s3, kvec.begin() + s3, kvec.begin() + 2 * s3, result.begin(), [](double k1, double k2, double k3) {return std::sqrt(k1 * k1 + k2 * k2 + k3 * k3);});
	for (int j = 0; j < dim; j++)
	{
		std::transform(result.begin(), result.end(), kvec.begin() + j * sn, result.begin(), [](double r, double k) {return r + k * k;});
	}
	std::transform(result.begin(), result.end(), result.begin(), [](double r) {return std::sqrt(r);});
	return result;
}

std::vector<double> build_k_bin_edges(int sidelength, double side)
{
	double dcell = side / sidelength;
	double kmin = 2 * pi.real() / side;
	double kmax = pi.real() / dcell;
	double dk_bin = kmin;
	std::vector<double> result;
	for (int i = 0; kmin + dk_bin * (i - 0.5) < kmax + dk_bin; i++)
	{
		result.push_back(kmin + dk_bin * (i - 0.5));
	}
	return result;
}

std::vector<double> bincount(std::vector<int> const& x, std::vector<double> const& weights, int minlength = 0)
{
	int maximum = std::accumulate(x.begin(), x.end(), 0, [](int a, int b) {return std::max(a, b);});
	std::vector<double> result;
	if (maximum + 1 > minlength)
	{
		result.resize(maximum + 1);
	}
	else
	{
		result.resize(minlength);
	}
	for (int i = 0; i < std::min((int)result.size(),minlength);i++)
	{
		int k = -1;
		result[i] = std::accumulate(x.begin(), x.end(), 0., [&](double a, int b) {++k; return b == i ? a + weights[k] : a;});
	}
	return result;
}

std::tuple<std::vector<double>, std::vector<double>> bin_isotropic_modes(std::vector<double> values, std::vector<double> kmod, int sidelength, double side)
{
	double dcell = side / sidelength;
	std::vector<double> bin_edges(build_k_bin_edges(sidelength, side));
	int n_bins = bin_edges.size() - 1;
	std::vector<int> bin_idx = digitize(kmod, bin_edges) - 1;
	std::vector<bool> valid = (bin_idx >= 0) && (bin_idx < n_bins) && (kmod > 0.);
	std::vector<double> weight(kmod.size(), 1.);

	std::vector<int> idx_valid = boolidx(bin_idx, valid);
	std::vector<double> k_valid = boolidx(kmod, valid);
	std::vector<double> values_valid = boolidx(values, valid);
	std::vector<double> w_valid = boolidx(weight, valid);

	std::vector<double> nmodes(bincount(idx_valid, w_valid, n_bins));
	std::vector<double> k_sum(bincount(idx_valid, w_valid * k_valid, n_bins));
	std::vector<double> values_sum(bincount(idx_valid, w_valid * values_valid, n_bins));

	std::vector<bool> populated = nmodes > 0.;

	std::vector<double> k_centres_populated = boolidx(k_sum, populated) / boolidx(nmodes, populated);
	std::vector<double> values_binned_populated = boolidx(values_sum, populated) / boolidx(nmodes, populated);

	return std::make_tuple(k_centres_populated, values_binned_populated);
}

