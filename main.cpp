#include <iostream>


typedef double Time;
typedef double Rate;

/********************************************************************************************
General diffusion process classes
This class describes a stochastic process governed by dx(t) = mu (t, x(t))dt + sigma(t, x(t))dz(t).
********************************************************************************************/

class DiffusionProcess
{
public: 
	DiffusionProcess(double x0) : x0_(x0) {}
	virtual ~DiffusionProcess() {}

	double x0() const { return x0_; }

	// return the drift part of the equation, i.e. mu(t, x_t)
	virtual double drift(Time t, double x) const = 0;

	//returns the diffusion part of the equation, i.e. sigma(t, x_t)
	virtual double diffusion(Time t, double x) const = 0;

	// returns the expectation of the process after a time interval
	// returns E(x_{t_0 + delta t} | x_{t_0} = x_0) since it is Markov.
	// By default, it returns the Euler approximation defined by 
	// x_0 + mu(t_0, x_0) delta t.
	virtual double expectation(Time t0, double x0, Time dt) const {
		return x0 + drift(t0,x0) * dt;
	}

	// returns the variance of the process after a time reversal
	// returns Var(x_{t_0 + Delta t} | x_{t_0} = x_0)
	// By default, it returns the Euler approximation defined by
	// sigma(t_0, x_0)^2 \Delta t
	virtual double variance(Time t0, double x0, Time dt) const {
		double sigma = diffusion(t0, x0);
		return sigma * sigma * dt;
	}
private:
	double x0_;


};


/****************************************************************
Black - Scholes diffusion process class
THis class describes the stochastic process goerned by dS = (r - 0.5 * {sigma ^2}) dt + sigmadz (t)
*****************************************************************/

class BlackScholesProcess : public DiffusionProcess
{
public: 
	BlackScholesProcess (Rate rate, double volatility, double s0 = 0.0) : DiffusionProcess(s0), r_(rate), sigma_(volatility) {}

	double drift(Time t, double x) const {
		return r_ - 0.5*sigma_*sigma_;
	}
	double diffusion(Time t, double x) const {
		return sigma_;
	}
private:
	double r_, sigma_;
};

bool BlackScholesProcess_expected_price_after_time_0_should_be_the_starting_price()
{
	const double START_PRICE = 12.1;
	const Rate RATE = 0.05;
	const double VOLATILITY = 0.06;
	BlackScholesProcess process(RATE, VOLATILITY);
	return START_PRICE == process.expectation(0, START_PRICE, 0);
}

int main() {
	std::cout << "BlackScholesProcess_expected_price_after_time_0_should_be_the_starting_price: " 
		      << BlackScholesProcess_expected_price_after_time_0_should_be_the_starting_price() << std::endl;
	return 0;
}




