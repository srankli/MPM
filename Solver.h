#ifndef __DATAANDSOLVER_H__
#define __DATAANDSOLVER_H__

class Output;
class OutputRequest;

// --------------------------------- SolverError ---------------------------------
// Error class of solver
class SolverError
{
protected:
	size_t code;
	std::string message;

	static const std::string errorCodeToExplanation[];
	static const size_t errorExplanationNum;

public:
	SolverError(size_t cd = 0, const char *msg = "");
	~SolverError() {}
	void setCode(unsigned int cd) { code = cd; }
	size_t getCode(void) { return code; }
	void setMessage(char *msg) { message = msg; }
	const std::string &getMessage(void) { return message; }
	const std::string &getExplaination(void);
};

// ----------------------------------- Solver -----------------------------------
// Each solver calculate one step of analysis
class Solver
{
	friend Output;
	friend OutputRequest;
protected:
	// The step from with this solver start (start from 1)
	size_t index;
	std::string name;

	// Start time of this step
	const double t_start;
	// Time length of this step
	const double t_step;
	// Time already calcuated (start from 0.0)
	double t_cal;
	// cureent time; (== t_start + t_cur)
	double t_cur;

	// Size of time increment
	double t_increment;
	/*
	t_increment_a = t_increment/2 at first interation, 
	t_increment_a = t_increment at other interations.
	*/
	double t_increment_a;
	unsigned long long iteration_index;
	
	// output class
	OutputRequest &output;
	
public:
	Solver(double time_step, OutputRequest &out, const char *na = "");
	Solver(double time_step, Solver &prev_solver, const char *na = "");
	~Solver() {}

	inline void setTimeIncrement(double time_increment) noexcept
	{ t_increment = time_increment; } // use CFL in the future

	/*
	Initialize Calculation
	Note that behavior of init() function need to change according to step_index
	*/
	virtual int init(void) { return 0; }
	// complete one iteration
	virtual int iteration(void) = 0;
	// Solve the whole problem
	// CFL condition may be used in the future
	virtual int solve(double time_increment);

/*
// ----------- NodeVar list ----------
protected:
	NodeVar *stack_top;
	size_t stack_num;
public:
	inline void push_stack(NodeVar *item) noexcept
	{
		item->next_solver = stack_top;
		stack_top = item;
		++stack_num;
	}
	inline size_t get_num_stack() noexcept { return stack_num; }
	inline void reset_stack(void) noexcept
	{
		stack_top = nullptr;
		stack_num = 0;
	}
	inline NodeVar *get_top_stack(void) noexcept { return stack_top; }
	inline void start_stack(NodeVar **cur) noexcept { *cur = stack_top; }
	inline void next_stack(NodeVar **cur) noexcept { *cur = (*cur)->next_solver; }
*/
};

#endif