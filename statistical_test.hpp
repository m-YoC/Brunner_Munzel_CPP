#pragma once

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

namespace test {

	class function {
	public:
		inline static double gamma(double x) {
			return std::tgamma(x);
		}

		inline static double log_gamma(double x) {
			return std::lgamma(x);
		}

		inline static double beta(double a, double b) {
			return std::exp(log_beta(a, b));
		}

		inline static double log_beta(double a, double b) {
			return log_gamma(a) + log_gamma(b) - log_gamma(a + b);
		}

		template<unsigned int LOOP_MAXCOUNT = 100>
		inline static double regularized_incomplete_beta(double x, double a, double b) {

			/* 定義域 0 ≦ Re(x) ≦ 1 */
			if (x < 0. || 1. < x) return std::nan("1");

			/* ibeta関数は対称型．収束が速い方を用いる */
			if (x * (a + b + 2.) > (a + 1.)) {
				return 1. - regularized_incomplete_beta<LOOP_MAXCOUNT>(1. - x, b, a);
			}

			/* ribeta(x, a, b) = res1 * res2 */

			/* 非反復部分の計算 発散を避けるために対数で計算を行う
			   res1 = x^a (1-x)^b / [a B(a,b)]
			*/
			const double res1 = exp(a * log(x) + b * log(1. - x) - log_beta(a, b)) / a;

			/* 反復部分の計算
			   a_(2m+1) = -(a+m)(a+b+m)x / [(a+2m)(a+2m+1)]
			   a_(2m  ) = m(b-m)x / [(a+2m-1)(a+2m)]
			   b = 1

			   P_-1 = 1, P_0 = b, Q_-1 = 0, Q_0 = 1
			   P_n = b P_(n-1) + a_n P_(n-2)
			   Q_n = b Q_(n-1) + a_n Q_(n-2)

			   f_n = P_n / Q_n

			   res2 = 1 / f_n
			*/

			std::array<double, 4> P = {0, 0, 1, 1};
			std::array<double, 4> Q = {0, 0, 0, 1};
			const static int n = 2;
			for (int i = 0; i < LOOP_MAXCOUNT; ++i) {
				/*shift*/
				P[n - 2] = P[n]; P[n - 1] = P[n + 1];
				Q[n - 2] = Q[n]; Q[n - 1] = Q[n + 1];

				const double m1 = i;
				const double a1 = -(a + m1) * (a + b + m1) * x / ((a + 2. * m1 ) * (a + 2. * m1 + 1.));

				/*漸化式*/
				P[n] = P[n - 1] + a1 * P[n - 2];
				Q[n] = Q[n - 1] + a1 * Q[n - 2];

				
				const double m2 = i + 1;
				const double a2 = m2 * (b - m2) * x / ((a + 2. * m2 - 1.) * (a + 2. * m2));

				/*漸化式*/
				P[n + 1] = P[n] + a2 * P[n - 1];
				Q[n + 1] = Q[n] + a2 * Q[n - 1];

			}

			const double res2 = Q[n + 1] / P[n + 1];
			
			return res1 * res2;
		}
	};

	class cdf {
	public:
		/*正規分布　累積分布関数*/
		static double n_cdf(double t) {
			return 0.5 * (1.0 + std::erf(t * 0.70710678118654752440084));
		}

		/*t分布　累積分布関数*/
		static double t_cdf(double t, double v) {
			const double y = sqrt(t * t + v);
			const double x = (t + y) / (2.0 * y);
			return function::regularized_incomplete_beta(x, 0.5 * v, 0.5 * v);
		}
	};


	typedef enum {
		normal_dist,
		t_dist
	} dist_type;

	template<typename scalar>
	class test_base {	
	public:
		using data = std::vector<scalar>;

		virtual scalar pvalue_greater() = 0;
		virtual scalar pvalue_less() = 0;
		virtual scalar pvalue_both() = 0;
	protected:
		scalar p;

		scalar getMean(const data& x) {
			scalar sum = .0;
			for (int i = 0; i < x.size(); ++i) {
				sum += x[i];
			}
			return sum / (scalar)x.size();
		}

		scalar getVariance(const data& x, scalar mean, int ddof = 0) {
			scalar sum = .0;
			for (int i = 0; i < x.size(); ++i) {
				sum += (x[i] - mean) * (x[i] - mean);
			}
			return sum / (scalar)(x.size() - ddof);
		}

		/*vectorのmergeメソッド*/
		template<typename type>
		static std::vector<type> merge(std::vector<std::vector<type>*> elems) {

			std::vector<type> res;

			size_t count = 0;
			for (auto v : elems) count += v->size();

			res.reserve(count);
			for (int i = 0; i < elems.size(); ++i) {
				std::copy(elems[i]->begin(), elems[i]->end(), std::back_inserter(res));
			}

			return std::move(res);
		}
	};


	/* ddof = 0 [標本分散] / 1 [不偏分散] */
	template<typename scalar>
	class welch_t : public test_base<scalar> {
	public:
		welch_t(data& x, data& y, const int ddof = 1) {
			if (x.size() <= 1 || y.size() <= 1) {
				p = -1;
				return;
			}

			auto meanx = getMean(x);
			auto meany = getMean(y);

			auto S2x = getVariance(x, meanx, ddof);
			auto S2y = getVariance(y, meany, ddof);

			scalar nx = (scalar)x.size();
			scalar ny = (scalar)y.size();

			scalar S2nx = S2x / nx;
			scalar S2ny = S2y / ny;

			scalar t = (meanx - meany) / sqrt(S2nx + S2ny);

			scalar v = pow(S2nx + S2ny, 2) / (S2nx * S2nx / (nx - static_cast<scalar>(1.0)) + S2ny * S2ny / (ny - static_cast<scalar>(1.0)));

			p = cdf::t_cdf(t, v);

		}

		virtual scalar pvalue_greater() override {
			return static_cast<scalar>(1.0) - p;
		}
		virtual scalar pvalue_less() override {
			return p;
		}
		virtual scalar pvalue_both() override {
			return static_cast<scalar>(2.0) * min(p, static_cast<scalar>(1.0) - p);
		}

	private:

	};

	template<typename scalar>
	class brunner_munzel : public test_base<scalar> {
		using rankdata = std::vector<scalar>;
		using tag = unsigned int;
		using data_tag = std::vector<std::pair<scalar, tag>>;

	public:
		brunner_munzel(data& x, data& y, dist_type dtype = t_dist) {
			if (x.size() <= 1 || y.size() <= 1) {
				p = -1;
				return;
			}

			auto rankcxy = createRank(merge<scalar>({ &x, &y }));
			auto rankcx_mean = getMean(rankcxy, 0, x.size());
			auto rankcy_mean = getMean(rankcxy, x.size(), y.size());

			auto rankx = createRank(x);
			auto ranky = createRank(y);
			auto rankx_mean = ((scalar)(x.size() + 1)) * static_cast<scalar>(0.5);
			auto ranky_mean = ((scalar)(y.size() + 1)) * static_cast<scalar>(0.5);

			scalar S2x = 0, S2y = 0;
			for (int i = 0; i < x.size(); ++i) {
				S2x += std::pow((rankcxy[i] - rankcx_mean) - (rankx[i] - rankx_mean), 2);
			}
			S2x /= (scalar)(x.size() - 1);

			for (int i = 0; i < y.size(); ++i) {
				S2y += std::pow((rankcxy[i + x.size()] - rankcy_mean) - (ranky[i] - ranky_mean), 2);
			}
			S2y /= (scalar)(y.size() - 1);

			/*2つの分布が完全分離している*/
			if (S2x == 0 && S2y == 0) {
				/*互いの最初の値の大小で評価してしまう*/
				if (x[0] < y[0]) {
					p = static_cast<scalar>(1.0);
				}
				else if (x[0] > y[0]) {
					p = static_cast<scalar>(0.0);
				}
				else {
					p = static_cast<scalar>(0.5);
				}
				return;
			}

			scalar wbfn = (scalar)x.size() * (scalar)y.size() * (rankcy_mean - rankcx_mean);
			wbfn /= ((scalar)x.size() + (scalar)y.size()) * std::sqrt((scalar)x.size() * S2x + (scalar)y.size() * S2y);

			if (dtype == t_dist) {
				p = calcP_t(wbfn, (scalar)x.size(), S2x, (scalar)y.size(), S2y);
			}
			else if (dtype == normal_dist) {
				p = calcP_normal(wbfn);
			}
			else {
				p = -1;
			}
		}

		virtual scalar pvalue_greater() override {
			return p;
		}
		virtual scalar pvalue_less() override {
			return static_cast<scalar>(1.0) - p;
		}
		virtual scalar pvalue_both() override {
			return static_cast<scalar>(2.0) * min(p, static_cast<scalar>(1.0) - p);
		}

	private:

		data_tag copyDataWithTag(const data& x) {
			data_tag res(x.size());
			for (int i = 0; i < x.size(); ++i) {
				res[i].first = x[i];
				res[i].second = i;
			}

			return res;
		}

		rankdata createRank(const data& x) {
			rankdata res(x.size());
			data_tag _x = copyDataWithTag(x);
			std::sort(_x.begin(), _x.end(), [](auto l, auto r) {
				return l.first < r.first;
			});

			for (int i = 0; i < _x.size();) {
				/*同順位のもので平均を取る*/
				int c = 0;
				/*以降のデータで同値の数を数える*/
				for (int j = i; j < _x.size(); ++j) {
					/*自身は必ず含まれる*/
					if (_x[j].first == _x[i].first) {
						c++;
					}
					else {
						break;
					}
				}
				/*
				{ Σ(n + n+1 + n+2 + ... + n+c-1) } / c
				= n + Σ(0 + 1 + 2 + 3 + ... + c-1) / c
				= n + c(c-1) / 2c
				= n + (c-1) / 2
				*/
				const scalar ave = (scalar)(i + 1) + (scalar)(c - 1) * static_cast<scalar>(0.5);
				for (int j = i; j < i + c; ++j) {
					res[_x[j].second] = ave;
				}
				i += c;
			}

			return res;
		}

		scalar getMean(const rankdata& x, unsigned begin, unsigned size) {
			scalar res = 0;
			for (int i = begin; i < begin + size; ++i) {
				res += (scalar)x[i];
			}
			return res / (scalar)size;
		}

		/*greaterのときのp値*/
		scalar calcP_normal(scalar wbfn) {
			return cdf::n_cdf(wbfn);
		}

		/*greaterのときのp値*/
		scalar calcP_t(scalar wbfn, scalar nx, scalar S2x, scalar ny, scalar S2y) {
			scalar df = std::pow(nx * S2x + ny * S2y, 2) / (std::pow(nx * S2x, 2) / (nx - static_cast<scalar>(1.0)) + std::pow(ny * S2y, 2) / (ny - static_cast<scalar>(1.0)));
			return cdf::t_cdf(wbfn, df);
		}

	};

}