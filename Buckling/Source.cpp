#include <boost\numeric\odeint.hpp>
#include <SFML\Graphics.hpp>
#include <iostream>
#include <vector>
#include <cmath>

using ldouble = long double;
using uint = unsigned short;
using state_type = std::vector<ldouble>;

const ldouble PI = acos(-1.0L);
const double stepSize = 0.001;
const ldouble precision = 1e-10;
ldouble crackMoment = 0;
const float scaleFactor = 75;


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
ldouble crackedColumn(state_type BC1, push_back_state_and_pos *Observer = nullptr, push_back_state_and_pos *BottomObserver = nullptr, properties param = param1);

void solveColumn(state_type BC1, ldouble loadFactor, push_back_state_and_pos *Observer, size_t *steps, properties param = param1, ldouble initPos = 0);
void plot(std::vector<Coord> &points, sf::RenderWindow &window, Coord origin, Coord scale = {1,1}, sf::Color pColor = sf::Color::White);
void getState(ldouble pos, push_back_state_and_pos &Data, state_type &state);

template <typename T = ldouble>
ldouble Exp(T x, int n);

void printCrackDebug(state_type initState, ldouble pMoment, ldouble angleKink, ldouble lowerAngle, state_type lowerBC, push_back_state_and_pos out);

int main()
{
	state_type BC1(2);
	const ldouble init_Angle = 60L;
	const ldouble init_Curve = 0L;
	ldouble loadRatio = 1;

	sf::RenderWindow window(sf::VideoMode(1280,900), "Column");

	BC1[0] = init_Angle * PI / 180L;
	BC1[1] = init_Curve;
	
	std::vector<state_type> finalState;
	std::vector<double> finalPos;
	push_back_state_and_pos finalObserver(finalState, finalPos);

	loadRatio = simpleBuckling(BC1, &finalObserver);
	
	{
		state_type testState;
		std::vector<state_type> testState1;
		std::vector<double> testPos1;
		push_back_state_and_pos testObs(testState1, testPos1);
		size_t testSteps;

		solveColumn(BC1, 1, &testObs, &testSteps, param1);

		//getState(s0, testObs, testState);
		//std::cout << '\n' << testState[0] << '\t' << testState[1];
	}

	std::vector<state_type> crackFinalStateTop, crackFinalStateBottom;
	std::vector<double> crackFinalPosTop, crackFinalPosBottom;
	push_back_state_and_pos crackFinalObserverTop(crackFinalStateTop, crackFinalPosTop);
	push_back_state_and_pos crackFinalObserverBottom(crackFinalStateBottom, crackFinalPosBottom);
	ldouble crackLoadRatio = crackedColumn(BC1, &crackFinalObserverTop, &crackFinalObserverBottom);

//	std::cout << std::endl << crackFinalStateTop.back()[0]*180/PI << '\t' << crackFinalStateBottom.front()[0] * 180 / PI << '\n';
//	std::cout << 180*(crackFinalStateTop.back()[0] - crackFinalStateBottom.front()[0]) << std::endl;

	std::cout << std::setprecision(10) << "Load Ratio without Crack: " << loadRatio;
	std::cout << std::setprecision(10) << "\nLoad Ration with Crack: " << crackLoadRatio;
	//getState();


	std::vector<Coord> points, crackPointsTop, crackPointsBottom;
	
	Coord origin(150, 800);
	Coord crackOrigin(150, 800);
	
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

			++stateIter;
			prevCoord = tempCoord;			
		}
	}

	{
		std::vector<state_type>::reverse_iterator stateIter = crackFinalStateBottom.rbegin();
		//std::vector<double>::reverse_iterator posIter = finalPos.rbegin();
		Coord prevCoord(0, 0), tempCoord;
		ldouble angle = 0, dl, totL = 0;
		while (stateIter != crackFinalStateBottom.rend())
		{
			angle = (*stateIter)[0];
			//std::cout << angle << '\n';
			tempCoord.x = prevCoord.x + stepSize * sin(angle);
			tempCoord.y = prevCoord.y + stepSize * cos(angle);

			crackPointsBottom.push_back(Coord(tempCoord.x, tempCoord.y));

			++stateIter;
			prevCoord = tempCoord;
		}

		stateIter = crackFinalStateTop.rbegin();
		
		while (stateIter != crackFinalStateTop.rend())
		{
			angle = (*stateIter)[0];
			//std::cout << angle << '\n';
			tempCoord.x = prevCoord.x + stepSize * sin(angle);
			tempCoord.y = prevCoord.y + stepSize * cos(angle);

			crackPointsTop.push_back(Coord(tempCoord.x, tempCoord.y));

			++stateIter;
			prevCoord = tempCoord;
		}
	}

	

	while (window.isOpen())
	{
		sf::Event winEvent;

		while (window.pollEvent(winEvent))
		{
			if (winEvent.type == sf::Event::Closed)
				window.close();
		}

		window.clear(sf::Color::Black);

		plot(points, window, origin, { scaleFactor,scaleFactor });
		plot(crackPointsBottom, window, crackOrigin, { scaleFactor,scaleFactor }, sf::Color::Green);
		plot(crackPointsTop, window, crackOrigin, { scaleFactor,scaleFactor },sf::Color::Red);

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

	while (inLoop)
	{
		for (std::vector<push_back_state_and_pos>::iterator iter = Obs.begin(); iter != Obs.end(); ++iter)
			iter->clearThis();

		solveColumn(BC1, k1, &Obs[0], &steps1, param);
		solveColumn(BC1, k2, &Obs[1], &steps2, param);

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

ldouble crackedColumn(state_type BC1, push_back_state_and_pos *TopObserver, push_back_state_and_pos *BottomObserver, properties param)
{
	state_type BC(2);

	const int NO_OF_OBSERVERS = 6;
	std::vector<std::vector<state_type>> state(NO_OF_OBSERVERS);
	std::vector<std::vector<double>> pos(NO_OF_OBSERVERS);

	std::vector<push_back_state_and_pos> Obs;
	for (int i = 0; i < NO_OF_OBSERVERS; ++i)
		Obs.push_back(push_back_state_and_pos(state[i], pos[i]));

	bool inLoop = true;
	std::vector<size_t> steps(NO_OF_OBSERVERS);
	ldouble k, k1, k2;
	k1 = 0L;
	k2 = 10L;

	ldouble Compliance = 7000*Exp<>(param.B,4);
	Compliance -= 6540 * Exp<>(param.B,3)*param.w;
	Compliance += 3665 * param.B*param.B*param.w*param.w;
	Compliance -= 700 * param.B*Exp<>(param.w, 3);
	Compliance += 561 * Exp<>(param.w, 4);
	Compliance *= Compliance * PI * param.a*9;
	Compliance /= (31250* param.B*Exp<>(param.w,12)*param.E);
	Compliance *= param.B;
	//Compliance += param.L / (param.E*param.I);

	std::cout << "\nCompliance: " << Compliance;

	while (inLoop)
	{
		state_type tempState = { 0,0 };
		state_type lowerBC = { 0,0 };
		ldouble pMoment = 0, angleKink = 0, lowerAngle = 0;
		properties upperParam = param;
		upperParam.L = s0;

		for (std::vector<push_back_state_and_pos>::iterator iter = Obs.begin(); iter != Obs.end(); ++iter)
			iter->clearThis();

		//Solving for case 1
		solveColumn(BC1, k1, &Obs[0], &steps[0], upperParam);
		getState(s0, Obs[0], tempState);

		pMoment = param.E*param.I*tempState[1];
		angleKink = pMoment*Compliance;
		lowerAngle = tempState[0]+angleKink;
		lowerBC = {lowerAngle, tempState[1]};
		
		solveColumn(lowerBC, k1, &Obs[1], &steps[1], param, s0);

		/*
		Coord tempCoord = { 0,0 };
		{
			std::vector<state_type>::reverse_iterator stateIter = Obs[0].states.rbegin();
			ldouble angle = 0, dl, totL = 0;

			while (stateIter != Obs[0].states.rend())
			{
				angle = (*stateIter)[0];
				//std::cout << angle << '\n';
				tempCoord.x += stepSize * sin(angle);
				tempCoord.y += stepSize * cos(angle);

				++stateIter;
			}
		}
		std::cout << k1*param.P*tempCoord.x << '\t' << tempCoord.y;
		*/

		printCrackDebug(tempState,pMoment,angleKink,lowerAngle,lowerBC, Obs[1]);

		//Solving for case 2
		solveColumn(BC1, k2, &Obs[2], &steps[2], upperParam);
		getState(s0, Obs[2], tempState);
		
		pMoment = param.E*param.I*tempState[1];
		angleKink = pMoment * Compliance;
		lowerAngle = tempState[0] + angleKink;
		lowerBC = { lowerAngle, tempState[1] };
		
		solveColumn(lowerBC, k2, &Obs[3], &steps[3], param, s0);
		printCrackDebug(tempState, pMoment, angleKink, lowerAngle, lowerBC, Obs[3]);

		//Solving for Midpoint
		k = (k1 + k2) / 2;
		solveColumn(BC1, k, &Obs[4], &steps[4], upperParam);
		getState(s0, Obs[4], tempState);

		pMoment = param.E*param.I*tempState[1];
		angleKink = pMoment * Compliance;
		lowerAngle = tempState[0] + angleKink;
		lowerBC = { lowerAngle, tempState[1] };

		solveColumn(lowerBC, k, &Obs[5], &steps[5], param, s0);
		printCrackDebug(tempState, pMoment, angleKink, lowerAngle, lowerBC, Obs[5]);

		//Biseciton Check
		bool k1Negative = (Obs[1].states.back()[0] < 0);

		if (abs(Obs[5].states.back()[0]) < precision)
			inLoop = false;

		if (k1Negative)
		{
			if (Obs[5].states.back()[0] < 0)
				k1 = k;
			else if (Obs[5].states.back()[0] > 0)
				k2 = k;
		}

		else
		{
			if (Obs[5].states.back()[0] < 0)
				k2 = k;
			else if (Obs[5].states.back()[0] > 0)
				k1 = k;
		}

		std::cout << "\n~~~~~~~~~~~~~~";
		std::cin.get();

	}

	const int top = 4, bottom = 5;
	if (TopObserver != nullptr)
	{
		TopObserver->clearThis();
		//upper part of column
		for (std::vector<state_type>::iterator stateIter = Obs[top].states.begin(); stateIter != Obs[top].states.end(); ++stateIter)
			TopObserver->states.push_back(*stateIter);
		for (std::vector<double>::iterator posIter = Obs[top].pos.begin(); posIter != Obs[top].pos.end(); ++posIter)
			TopObserver->pos.push_back(*posIter);

	}
	if (BottomObserver != nullptr)
	{
		BottomObserver->clearThis();
		//Lower part of column
		for (std::vector<state_type>::iterator stateIter = Obs[bottom].states.begin(); stateIter != Obs[bottom].states.end(); ++stateIter)
			BottomObserver->states.push_back(*stateIter);
		for (std::vector<double>::iterator posIter = Obs[bottom].pos.begin(); posIter != Obs[bottom].pos.end(); ++posIter)
			BottomObserver->pos.push_back(*posIter);
	}



	return k;
}

void solveColumn(state_type BC1, ldouble loadFactor, push_back_state_and_pos * Observer, size_t *steps, properties param, ldouble initPos)
{
	boost::numeric::odeint::runge_kutta4< state_type > stepper;
	state_type BC(2);

	BC = BC1;
	C = -(loadFactor* param.P) / (param.E*param.I);
	
	*steps = boost::numeric::odeint::integrate_const(stepper, ColumnEquation, BC, static_cast<double>(initPos), param.L + stepSize, stepSize, *Observer);
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

void plot(std::vector<Coord> &points, sf::RenderWindow &window, Coord origin, Coord scale, sf::Color pColor)
{
	sf::VertexArray Curve(sf::LineStrip);

	for (std::vector<Coord>::iterator iter = points.begin(); iter != points.end(); ++iter)
		Curve.append(sf::Vertex(iter->mapScreen(origin, scale)));
	for (int i = 0; i < Curve.getVertexCount(); ++i)
		Curve[i].color = pColor;
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

void printCrackDebug(state_type initState, ldouble pMoment, ldouble angleKink, ldouble lowerAngle, state_type lowerBC, push_back_state_and_pos out)
{
	std::cout << "\nInitial State - Q = " << initState[0] << " Q' = " << initState[1];
	std::cout << "\nCompliance = " << angleKink/pMoment;
	std::cout << "\nMoment = " << pMoment;
	std::cout << "\nangleKink = " << angleKink << '\t' << angleKink*180L/PI;
	std::cout << "\nlowerAngle = " << lowerAngle;
	std::cout << "\nBC at Crack - Q = " << lowerBC[0] << " Q' = " << lowerBC[1];
	std::cout << "\nEnd state - Q = " << out.states.back()[0] << " Q' = " << out.states.back()[1];
	std::cout << "\nFinal Pos = " << out.pos.back();
	std::cout << "\n=====================";

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