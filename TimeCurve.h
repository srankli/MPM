#ifndef __TIMECURVE_H__
#define __TIMECURVE_H__

#define PI 3.14159265358979323846264338327950288419716939937510

enum class TimeCurveType : unsigned int
{
	Invalid = 0,
	Zero = 1,
	Constant = 2,
	Linear = 3,
	CubicSmooth = 4
};

// No polymorphism
struct TimeCurve_inline
{
protected:
	TimeCurveType type;
public:
	TimeCurve_inline(TimeCurveType tp = TimeCurveType::Invalid) : type(tp) {}
	inline TimeCurveType getType() noexcept { return type; }
};

struct TimeCurveZero_inline : public TimeCurve_inline
{
	TimeCurveZero_inline() : TimeCurve_inline(TimeCurveType::Zero) {}
	inline double x(double t) noexcept { return 0.0; }
	inline double dx_dt(double t) noexcept { return 0.0; }
};

struct TimeCurveConstant_inline : public TimeCurve_inline
{
	double x_const;
	TimeCurveConstant_inline() : TimeCurve_inline(TimeCurveType::Constant) {}
	inline double x(double t) noexcept { return x_const; }
	inline double dx_dt(double t) noexcept { return 0.0; }
};

struct TimeCurveLinear_inline : public TimeCurve_inline
{
	double x_begin;
	double x_end;
	double t_len;
protected:
	double gradient;
public:
	TimeCurveLinear_inline() : TimeCurve_inline(TimeCurveType::Linear) {}
	inline void init(void) noexcept
	{
		gradient = (x_end - x_begin) / t_len;
	}
	inline double x(double t) noexcept
	{
		if (t < 0.0)
			return x_begin;
		else if (t > t_len)
			return x_end;
		else
			return x_begin + (x_end - x_begin) / t_len * t;
	}
	inline double dx_dt(double t) noexcept
	{
		if (t < 0.0 || t > t_len)
			return 0.0;
		else
			return gradient;
	}
};

struct TimeCurveCubicSmooth_inline : public TimeCurve_inline
{
	double x_begin;
	double x_end;
	double t_len;
protected:
	double half_t_len;
	// for x()
	double M;
	double M_div_3tlen;
	double M_div_2;
	double minus_M_pro_tlen2_div_6_plus_xend;
	// for dx_dt()
	double M_div_tlen;
public:
	TimeCurveCubicSmooth_inline() :
		TimeCurve_inline(TimeCurveType::CubicSmooth) {}
	inline void init(void) noexcept
	{
		half_t_len = t_len / 2.0;
		// for x()
		M = 6.0 * (x_end - x_begin) / (t_len * t_len);
		M_div_3tlen = M / (3.0 * t_len);
		M_div_2 = M / 2.0;
		minus_M_pro_tlen2_div_6_plus_xend = -M * t_len * t_len / 6.0 + x_end;
		// for dx_dt()
		M_div_tlen = M / t_len;
	}
	inline double x(double t) noexcept
	{
		double t2 = t * t;
		double t3 = t2 * t;
		if (t < 0.0)
			return x_begin;
		else if (t < half_t_len)
			return -M_div_3tlen * t3 + M_div_2 * t2 + x_begin;
		//return -M / (3.0*t_len) * t3 + M / 2.0 * t2 + x_begin;
		else if (t < t_len)
			return -M_div_3tlen * t3 + M_div_2 * t2 + minus_M_pro_tlen2_div_6_plus_xend;
			//return -M / (3.0*t_len) * t3 + M / 2.0 * t2 - M / 6.0 * t_len2 + x_end;
		else
			return x_end;
	}
	inline double dx_dt(double t) noexcept
	{
		if (t < 0.0 || t > t_len)
			return 0.0;
		else
			return (-M_div_tlen * t + M) * t;
	}
};


struct TimeCurveSine_inline : public TimeCurve_inline
{
	double amplitude;
	double period;
	/*
	    Initial phase, can adjust initial phase
	to make the time curve become cosine.
	*/
	double fai;
protected:
	double w; // 2 * pi / period

public:
	TimeCurveSine_inline() : fai(0.0) {}
	inline void init(void) noexcept
	{
		w = 2.0 * PI / period;
	}
	inline double x(double t) noexcept
	{
		return amplitude * sin(w * t + fai);
	}
	inline double dx_dt(double t) noexcept
	{
		return w * amplitude * cos(w * t + fai);
	}
};


#endif