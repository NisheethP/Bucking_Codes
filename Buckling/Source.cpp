#include <boost\numeric\odeint.hpp>
#include <SFML\Graphics.hpp>
#include <iostream>
#include <vector>
#include <cmath>

using ldouble = long double;
using uint = unsigned short;
using state_type = std::vector<ldouble>;

const ldouble PI = acos(-1.0L);
const double stepSize = 0.01;
const ldouble precision = 1e-6;
ldouble crackMoment = 0;


struct push_back_state_and_pos
{
	std::vector< state_type > &states;
	std::vector< double > &pos;

	push_back_state_and_pos(std::vector< state_type > &states, std::vector< double > &times)
		: states(states),
		pos(times)
	{ }

	void operator()(const state_type &x, double t)
	{
		states.push_back(x);
		pos.push_back(t);
	}

	void clearThis()
	{
		states.clear();
		pos.clear();
	}
};

struct properties
{
	double L;
	ldouble B;
	ldouble w;
	ldouble E;
	ldouble I;
	ldouble P;
	ldouble a;

	properties(double pL, ldouble pB, ldouble pW, ldouble pE, ldouble pA = 0) :
		L(pL), 
		B(pB),
		w(pW),
		E(pE),
		a(pA*w),
		I(B * w*w*w / 12L),
		P((PI*PI*E*I) / (4 * L*L))
	{}

	ldouble K1C(ldouble moment)
	{
		ldouble temp = 6*moment/B;
		temp /= (w*w);
		temp *= sqrt(PI*a);
		return temp;
	}

};

struct Coord
{
	float x;
	float y;
	
	Coord(float pX = 0, float pY = 0):
		x(pX),
		y(pY)
	{}

	Coord(sf::Vector2f v) :
		x(v.x),
		y(v.y)
	{}

	sf::Vector2f getVector();
	sf::Vector2f mapScreen(Coord origin, Coord scale = { 1,1 });
};

//properties param1(10.0,5L,0.05L,0.1L,200e9);
properties param1(10.0, 0.05L, 0.1L, 200e9, 0.5);
ldouble C = -(param1.P) / (param1.E*param1.I);

const ldouble s0 = 5;

void ColumnEquation(const state_type &theta, state_type &dqds, const double x);
void printData(push_back_state_and_pos Obs, uint steps1);

ldouble simpleBuckling(state_type BC1, push_back_state_and_pos *Observer = nullptr, properties param = param1);
ldouble crackedColumn(state_type BC1, push_back_state_and_pos *Observer = nullptr, properties param = param1);

void solveColumn(state_type BC1, ldouble loadFactor, push_back_state_and_pos *Observer, size_t *steps, properties param = param1);
void plot(std::vector<Coord> &points, sf::RenderWindow &window, Coord origin, Coord scale = {1,1});
void getState(ldouble pos, push_back_state_and_pos &Data, state_type &state);

template <typename T = ldouble>
ldouble Exp(T x, int n);


int main()
{
	state_type BC1(2);
	const ldouble init_Angle = 45L;
	const ldouble init_Curve = 0L;
	ldouble loadRatio = 1;

	sf::RenderWindow window(sf::VideoMode(1280,900), "Column");

	BC1[0] = init_Angle * PI / 180L;
	BC1[1] = init_Curve;
	
	std::vector<state_type> finalState;
	std::vector<double> finalPos;
	push_back_state_and_pos finalObserver(finalState, finalPos);

	loadRatio = simpleBuckling(BC1, &finalObserver);

	crackedColumn(BC1, &finalObserver);

	std::cout << loadRatio;
	
	std::vector<Coord> points;
	{
		std::vector<state_type>::reverse_iterator stateIter = finalState.rbegin();
		//std::vector<double>::reverse_iterator posIter = finalPos.rbegin();
		Coord prevCoord(0,0), tempCoord;
		ldouble angle = 0, dl, totL = 0;
		while (stateIter != finalState.rend())
		{
			angle = (*stateIter)[0];
			//std::cout << angle << '\n';
			tempCoord.x = prevCoord.x + stepSize * sin(angle);
			tempCoord.y = prevCoord.y + stepSize * cos(angle);

			points.push_back(Coord(tempCoord.x, tempCoord.y));
			ldouble dx2 = (tempCoord.x - prevCoord.x);
			ldouble dy2 = (tempCoord.y - prevCoord.y);
			dl = std::sqrt(dx2*dx2 + dy2*dy2);
			totL += dl;

			++stateIter;
			//++posIter;
			prevCoord = tempCoord;
			
		}
		std::cout << '\n' << param1.L << '\t' << totL;
		//state_type tempState;
		//getState(10, finalObserver, tempState);
		//std::cout << '\n' << tempState[0] << '\t' << tempState[1];
		
	}
	Coord origin(150,800);

	while (window.isOpen())
	{
		sf::Event winEvent;

		while (window.pollEvent(winEvent))
		{
			if (winEvent.type == sf::Event::Closed)
				window.close();
		}

		window.clear(sf::Color::Black);

		plot(points, window, origin, {50,50});

		window.display();
	}

	return 0;
}

ldouble simpleBuckling(state_type BC1, push_back_state_and_pos *Observer, properties param)
{
	state_type BC(2);

	const int NO_OF_OBSERVERS = 3;
	std::vector<std::vector<state_type>> state(NO_OF_OBSERVERS);
	std::vector<std::vector<double>> pos(NO_OF_OBSERVERS);

	std::vector<push_back_state_and_pos> Obs;
	for (int i = 0; i < NO_OF_OBSERVERS; ++i)
		Obs.push_back(push_back_state_and_pos(state[i],pos[i]));
	
	bool inLoop = true;
	size_t steps1, steps2, steps;
	ldouble k, k1, k2;
	k1 = 1L;
	k2 = 10L;

	while(inLoop)
	{
		for (std::vector<push_back_state_and_pos>::iterator iter = Obs.begin(); iter != Obs.end(); ++iter)
			iter->clearThis();
	
		solveColumn(BC1,k1,&Obs[0],&steps1,param);
		solveColumn(BC1, k2, &Obs[1], &steps2, param);
		
		k = (k1 + k2) / 2;
		solveColumn(BC1, k, &Obs[2], &steps, param);

		bool k1Negative = (Obs[0].states[steps1-1][0] < 0);
		
		if (abs(Obs[0].states[steps1 - 1][0]) < precision)
			inLoop = false;

		if (k1Negative)
		{
			if (Obs[2].states[steps - 1][0] < 0)
				k1 = k;
			else if (Obs[2].states[steps - 1][0] > 0)
				k2 = k;
		}

		else
		{
			if (Obs[2].states[steps - 1][0] < 0)
				k2 = k;
			else if (Obs[2].states[steps - 1][0] > 0)
				k1 = k;
		}
	}

	//std::cout << "The ratio P/Pcr is: " << std::setprecision(11) << k;	
	if (Observer != nullptr)
	{
		for (std::vector<state_type>::iterator stateIter = Obs[2].states.begin(); stateIter != Obs[2].states.end(); ++stateIter)
			Observer->states.push_back(*stateIter);
		for (std::vector<double>::iterator posIter = Obs[2].pos.begin(); posIter != Obs[2].pos.end(); ++posIter)
			Observer->pos.push_back(*posIter);
	}
	return k;

}

ldouble crackedColumn(state_type BC1, push_back_state_and_pos *Observer, properties param)
{
	state_type BC(2);

	const int NO_OF_OBSERVERS = 6;
	std::vector<std::vector<state_type>> state(NO_OF_OBSERVERS);
	std::vector<std::vector<double>> pos(NO_OF_OBSERVERS);

	std::vector<push_back_state_and_pos> Obs;
	for (int i = 0; i < NO_OF_OBSERVERS; ++i)
		Obs.push_back(push_back_state_and_pos(state[i], pos[i]));

	bool inLoop = true;
	size_t steps1, steps2, steps;
	ldouble k, k1, k2;
	k1 = 1L;
	k2 = 10L;

	ldouble Compliance = 7000*Exp<>(param.B,4);
	Compliance -= 6540 * Exp<>(param.B,3)*param.w;
	Compliance += 3665 * param.B*param.B*param.w*param.w;
	Compliance -= 700 * param.B*Exp<>(param.w, 3);
	Compliance += 561 * Exp<>(param.w, 4);
	Compliance *= Compliance * PI * param.a*9;
	Compliance /= (31250* param.B*Exp<>(param.w,12)*param.E);
	Compliance += param.L / (param.E*param.I);
	
	while (inLoop)
	{

		state_type tempState;

		for (std::vector<push_back_state_and_pos>::iterator iter = Obs.begin(); iter != Obs.end(); ++iter)
			iter->clearThis();

		//Solving for case 1
		solveColumn(BC1, k1, &Obs[0], &steps1, param);
		getState(s0, Obs[0], tempState);
		ldouble pMoment = param.E*param.I*tempState[1];
		ldouble angleKink = pMoment*Compliance;
		ldouble lowerAngle = tempState[0]+angleKink;
		state_type lowerBC = {lowerAngle, tempState[1]};

		solveColumn(BC1, k2, &Obs[1], &steps2, param);
		getState(s0, Obs[1], tempState);

		k = (k1 + k2) / 2;
		solveColumn(BC1, k, &Obs[2], &steps, param);

		bool k1Negative = (Obs[0].states[steps1 - 1][0] < 0);

		if (abs(Obs[0].states[steps1 - 1][0]) < precision)
			inLoop = false;

		if (k1Negative)
		{
			if (Obs[2].states[steps - 1][0] < 0)
				k1 = k;
			else if (Obs[2].states[steps - 1][0] > 0)
				k2 = k;
		}

		else
		{
			if (Obs[2].states[steps - 1][0] < 0)
				k2 = k;
			else if (Obs[2].states[steps - 1][0] > 0)
				k1 = k;
		}
	}
	return k;
}

void solveColumn(state_type BC1, ldouble loadFactor, push_back_state_and_pos * Observer, size_t *steps, properties param)
{
	boost::numeric::odeint::runge_kutta4< state_type > stepper;
	state_type BC(2);

	BC = BC1;
	C = -(loadFactor* param1.P) / (param1.E*param1.I);
	
	*steps = boost::numeric::odeint::integrate_const(stepper, ColumnEquation, BC, 0.0, param1.L + stepSize, stepSize, *Observer);
}

void ColumnEquation(const state_type &theta, state_type &dqds, const double t)
{	
	dqds[0] = theta[1];
	//double temp = theta[0];
	dqds[1] = C*sin(theta[0]);
}

void printData(push_back_state_and_pos Obs, uint steps1)
{
	for (uint i = 0; i < steps1; ++i)
	{
		std::cout << Obs.pos[i] << '\t' << Obs.states[i][0] << '\t' << Obs.states[i][1] << '\n';
	}
}

void plot(std::vector<Coord> &points, sf::RenderWindow &window, Coord origin, Coord scale )
{
	sf::VertexArray Curve(sf::LineStrip);

	for (std::vector<Coord>::iterator iter = points.begin(); iter != points.end(); ++iter)
		Curve.append(sf::Vertex(iter->mapScreen(origin, scale)));

	window.draw(Curve);
}

void getState(ldouble pos, push_back_state_and_pos &Data, state_type &state)
{
	state_type tempState;
	
	std::vector<state_type>::iterator stateIter = Data.states.begin();
	std::vector<double>::iterator posIter = Data.pos.begin();

	while ((posIter != Data.pos.end()) && (stateIter != Data.states.end()) && (*posIter != pos))
	{
		state = *stateIter;
		++posIter;
		++stateIter;
	}

}

sf::Vector2f Coord::getVector()
{
	return sf::Vector2f(x,y);
}

sf::Vector2f Coord::mapScreen(Coord origin, Coord scale)
{
	sf::Vector2f temp;
	temp.x = origin.x + this->x*scale.x;
	temp.y = origin.y - this->y*scale.y;

	return temp;
}

template <typename T>
ldouble Exp(T x, int n)
{
	ldouble temp = 1;
	for (int i = 0; i < n; ++i)
		temp *= x;

	return temp;
}