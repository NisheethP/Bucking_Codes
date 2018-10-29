#include <boost\numeric\odeint.hpp>
#include <SFML\Graphics.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

using ldouble = long double;
using uint = unsigned short;
using state_type = std::vector<ldouble>;

const ldouble PI = acos(-1.0L);
const double stepSize = 0.01;
const ldouble precision = 1e-7;
ldouble crackMoment = 0;
const float scaleFactor = 60;

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
properties param1(10.0L, 0.05L, 0.1L, 200e9, 0.75L);
properties CompParam(param1);

ldouble C = -(param1.P) / (param1.E*param1.I);

ldouble s0 = 5;

void ColumnEquation(const state_type &theta, state_type &dqds, const double x);
void ComplianceEquation(const state_type &Comp, state_type &dcda, const double a);

void printData(push_back_state_and_pos Obs, uint steps1);

ldouble simpleBuckling(state_type BC1, push_back_state_and_pos *Observer = nullptr, properties param = param1);
ldouble crackedColumn(state_type BC1, push_back_state_and_pos *Observer = nullptr, push_back_state_and_pos *BottomObserver = nullptr, properties param = param1, ldouble loadLimitLower = 0.0L, ldouble loadLimitUpper = 10.0L);

void solveColumn(state_type BC1, ldouble loadFactor, push_back_state_and_pos *Observer, size_t *steps, properties param = param1, ldouble initPos = 0);
void plot(std::vector<Coord> &points, sf::RenderWindow &window, Coord origin, Coord scale = {1,1}, sf::Color pColor = sf::Color::White);
void getState(ldouble pos, push_back_state_and_pos &Data, state_type &state);
ldouble SIF_Factor(ldouble x);		//x = a/w
ldouble SIF_Factor(properties prop);

template <typename T = ldouble>
ldouble Exp(T x, int n);

void printCrackDebug(state_type initState, ldouble pMoment, ldouble angleKink, ldouble lowerAngle, state_type lowerBC, push_back_state_and_pos out);

int AnimateMain();
Coord getPoints(std::vector<state_type> state, std::vector<Coord> &points,Coord initCoord = { 0,0 });
void outputToFile(push_back_state_and_pos &Observer, std::string fileName);

ldouble percentChange(ldouble init, ldouble fin);
int crackLength();

int main()
{
	AnimateMain();
	return 0;
}

int crackLength()
{
	std::fstream file("F:\\C++\\C++ Programs\\Fracture\\Crack Length Output\\output.txt", std::ios::out);
	if (!file.is_open())
	{
		std::cout << "File not opened";
		std::cin.get();
		return -1;
	}

	file << "Load Drop\t\t\tCurve Legnth\t\t\ta/w\t\t\tError" << std::endl;
	properties paramNeg(10.0L, 0.05L, 0.1L, 200e9, 0.0L);
	properties paramPos(10.0L, 0.05L, 0.1L, 200e9, 1.0L);
	properties paramTemp(paramNeg);

	const ldouble init_Angle = 80.0L;
	const ldouble init_Curve = 0.0L;
	state_type BC1(2);
	BC1[0] = init_Angle * PI / 180.0L;
	BC1[1] = init_Curve;

	ldouble loadRatio = 1, loadRatioElastic = 1;
	loadRatioElastic = simpleBuckling(BC1, nullptr, paramNeg);
	
	bool inLoop = true;
	s0 = 2;

	ldouble percDropNeg, percDropPos, percDropTemp;

	ldouble target = 10.0L;
	for (int i = 1; i <= 100; ++i)
	{
		std::cout << "\nLoad Drop = " << i << " started";
		target = i/100.0L;
		target *= 25.0L / 100.0L;
		//target *= 25.0L / 100.0L;
		paramNeg.a = paramTemp.a;
		paramPos.a = 1.0L*param1.w;
		paramTemp.a = 0;
		inLoop = true;

		int iterCount = 0;
		ldouble error = 1;
		while (inLoop)
		{
			
			loadRatio = crackedColumn(BC1, nullptr, nullptr, paramNeg, 0, loadRatioElastic);
			percDropNeg = percentChange(loadRatioElastic, loadRatio);

			loadRatio = crackedColumn(BC1, nullptr, nullptr, paramPos, 0, loadRatioElastic);
			percDropPos = percentChange(loadRatioElastic, loadRatio);

			paramTemp.a = (paramNeg.a + paramPos.a) / 2.0L;

			loadRatio = crackedColumn(BC1, nullptr, nullptr, paramTemp, 0, loadRatioElastic);
			percDropTemp = percentChange(loadRatioElastic, loadRatio);

			error = (percDropTemp - target);

			if (0 > (error))
				paramNeg = paramTemp;
			if (0 < (error))
				paramPos = paramTemp;
			if (abs(error) <= 1e-4/2)
				inLoop = false;
			//std::cout << "\nError: " << error;
			//std::cout << std::setprecision(7) << "\n" << percDropNeg- target << '\t' << percDropPos- target << '\t' << percDropTemp << '\t' << paramTemp.a;
			//std::cin.get();

			if (iterCount > 1e4)
			{
				std::cout << "\nIteratons over " << iterCount << " for target = " << i << "; Final Error = " << error;
				break;
			}

			++iterCount;
			//std::cout << "\t" << iterCount;
		}
		//std::cout << "a = " << paramTemp.a;
		//std::cin.get();
		file << std::setprecision(8) << target * 100 << "\t\t\t\t" << paramTemp.a << "\t\t\t\t" << paramTemp.a / paramTemp.w << "\t\t\t" << error << "\t\t\t\t" << std::endl;
	}
	
	return 0;
}

int AnimateMain()
{

	s0 = 5;
	state_type BC1(2);
	const ldouble init_Angle = 80L;
	const ldouble init_Curve = 0L;
	ldouble loadRatio = 1;

	BC1[0] = init_Angle * PI / 180L;
	BC1[1] = init_Curve;
	
	std::vector<state_type> finalState;
	std::vector<double> finalPos;
	push_back_state_and_pos finalObserver(finalState, finalPos);

	loadRatio = simpleBuckling(BC1, &finalObserver);
	
	/*
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
	*/
	
	std::vector<state_type> crackFinalStateTop, crackFinalStateBottom;
	std::vector<double> crackFinalPosTop, crackFinalPosBottom;
	push_back_state_and_pos crackFinalObserverTop(crackFinalStateTop, crackFinalPosTop);
	push_back_state_and_pos crackFinalObserverBottom(crackFinalStateBottom, crackFinalPosBottom);
	
	ldouble crackLoadRatio = crackedColumn(BC1, &crackFinalObserverTop, &crackFinalObserverBottom, param1, 0, loadRatio);

	//std::cout << std::endl << crackFinalStateTop.back()[0]*180/PI << '\t' << crackFinalStateBottom.front()[0] * 180 / PI << '\n';
	//std::cout << 180*(crackFinalStateTop.back()[0] - crackFinalStateBottom.front()[0]) << std::endl;

	std::cout << std::setprecision(10) << "\nLoad Ratio without Crack: " << loadRatio;
	std::cout << std::setprecision(10) << "\nLoad Ration with Crack: " << crackLoadRatio;

	std::vector<Coord> points, crackPointsTop, crackPointsBottom;
	 
	std::vector<std::vector<Coord>> pointsFrame, crackPointsFrameTop, crackPointsFrameBottom;
	const int numFrames = 100;
	const ldouble minAngle = 0L;
	
	//Simple Column
	for (int i = 0; i < numFrames; i++)
	{
		state_type frameBC(BC1);
		frameBC[0] = minAngle + i*((init_Angle-minAngle)/ numFrames);
		frameBC[0] *= PI / 180L;
		
		std::vector<state_type> frameState;
		std::vector<double> framePos;
		push_back_state_and_pos frameObserver(frameState, framePos);

		simpleBuckling(frameBC, &frameObserver);
		std::vector<Coord> tempPoints;
		
		getPoints(frameState, tempPoints);

		pointsFrame.push_back(tempPoints);
		//std::cout << "\n Iteration " << i << "done";
	}

	//Cracked Column
	for (int i = 0; i < numFrames; i++)
	{
		state_type frameBC(BC1);
		frameBC[0] = minAngle + i * ((init_Angle - minAngle) / numFrames);
		frameBC[0] *= PI / 180L;

		std::vector<state_type> frameStateTop;
		std::vector<double> framePosTop;
		push_back_state_and_pos frameObserverTop(frameStateTop, framePosTop);

		std::vector<state_type> frameStateBottom;
		std::vector<double> framePosBottom;
		push_back_state_and_pos frameObserverBottom(frameStateBottom, framePosBottom);

		crackedColumn(frameBC, &frameObserverTop, &frameObserverBottom);

		std::vector<Coord> tempPoints;
		Coord crackCoord(0,0);

		crackCoord = getPoints(frameStateBottom, tempPoints);
		crackPointsFrameBottom.push_back(tempPoints);

		getPoints(frameStateTop, tempPoints, crackCoord);
		crackPointsFrameTop.push_back(tempPoints);
	}

	Coord origin(150, 850);
	Coord crackOrigin(350, 850);
	
	//Simple Column
	getPoints(finalState, points);

	//Cracked Column
	Coord tempCrackCoords = getPoints(crackFinalStateBottom, crackPointsBottom);
	getPoints(crackFinalStateTop, crackPointsTop, tempCrackCoords);

	
	pointsFrame.push_back(points);
	crackPointsFrameTop.push_back(crackPointsTop);
	crackPointsFrameBottom.push_back(crackPointsBottom);
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//GRAPHICS SECTION
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	
	sf::RenderWindow window(sf::VideoMode(1280, 900), "Column");
	auto iter = pointsFrame.begin();
	auto cIterTop = crackPointsFrameTop.begin();
	auto cIterBottom = crackPointsFrameBottom.begin();
	while (window.isOpen())
	{
		if (iter == pointsFrame.end())
		{
			iter = pointsFrame.begin();
			sf::sleep(sf::milliseconds(500));
		}

		if (cIterTop == crackPointsFrameTop.end())
		{
			cIterTop = crackPointsFrameTop.begin();
			sf::sleep(sf::milliseconds(500));
		}

		if (cIterBottom == crackPointsFrameBottom.end())
		{
			cIterBottom = crackPointsFrameBottom.begin();
			sf::sleep(sf::milliseconds(500));
		}

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
		
		if (iter != pointsFrame.end())
			plot(*iter, window, origin, { scaleFactor,scaleFactor });

		if (cIterTop != crackPointsFrameTop.end())
			plot(*cIterTop, window, crackOrigin, { scaleFactor,scaleFactor }, sf::Color::Red);
		
		if (cIterBottom != crackPointsFrameBottom.end())
			plot(*cIterBottom, window, crackOrigin, { scaleFactor,scaleFactor }, sf::Color::Green);
		
		window.display();
		
		sf::sleep(sf::milliseconds(10000/numFrames));
		++iter;
		++cIterTop;
		++cIterBottom;
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
	int i = 0;

	while (inLoop)
	{
		++i;
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

	//std::cout << "\nColumn solven in: " << i << " iterations.";
	return k;

}

ldouble crackedColumn(state_type BC1, push_back_state_and_pos *TopObserver, push_back_state_and_pos *BottomObserver, properties param, ldouble loadLimitLower, ldouble loadLimitUpper )
{
	bool crackDebugOutput = false;

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
	k1 = loadLimitLower;
	k2 = loadLimitUpper;

	ldouble Compliance = 0;
	{
		std::vector<state_type> compType;
		std::vector<double> compPos;
		push_back_state_and_pos compObs(compType, compPos);

		boost::numeric::odeint::runge_kutta_fehlberg78< state_type > stepper;

		state_type compBC(1);
		compBC[0] = 0;

		//stepper, ComplianceEquation, BC, 0, s0, 0.0001, compObs
		boost::numeric::odeint::integrate_const(stepper, ComplianceEquation, compBC, 0.0, static_cast<double>(param.a), 0.00001, compObs);
		Compliance = compType.back()[0];
	}

	//std::cout << "\nCompliance: " << Compliance << std::endl;

	int i = 0;

	while (inLoop)
	{
		++i;
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

		if (crackDebugOutput)
			printCrackDebug(tempState,pMoment,angleKink,lowerAngle,lowerBC, Obs[1]);

		//Solving for case 2
		solveColumn(BC1, k2, &Obs[2], &steps[2], upperParam);
		getState(s0, Obs[2], tempState);
		
		pMoment = param.E*param.I*tempState[1];
		angleKink = pMoment * Compliance;
		lowerAngle = tempState[0] + angleKink;
		lowerBC = { lowerAngle, tempState[1] };
		
		solveColumn(lowerBC, k2, &Obs[3], &steps[3], param, s0);
		if (crackDebugOutput)
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
		if (crackDebugOutput)
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

		if(crackDebugOutput)
		{
			std::cout << "\n~~~~~~~~~~~~~~";
			std::cin.get();
		}
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

	//std::cout << "\nCracked Column solven in: " << i << " iterations.";
	return k;
}

void solveColumn(state_type BC1, ldouble loadFactor, push_back_state_and_pos * Observer, size_t *steps, properties param, ldouble initPos)
{
	//runge_kutta_fehlberg78
	boost::numeric::odeint::runge_kutta4 < state_type > stepper;
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

void ComplianceEquation(const state_type &Comp, state_type &dcda, const double a)
{
	//ldouble K1c = CompParam.K1C();
	ldouble temp = 72 * PI*a;
	temp /= (CompParam.E*CompParam.B*Exp<>(CompParam.w,4));
	ldouble F = SIF_Factor(a / CompParam.w);
	temp *= F*F;
	dcda[0] = temp;

	//std::cout << "\na: " << a << "\tComp: " << Comp[0] << "\tdcda " << dcda[0] << "\tF: " << F;
	//std::cin.get();
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

ldouble SIF_Factor(ldouble x)
{
	ldouble F = 1.222 - 1.4*x + 7.33*x*x - 13.08*x*x*x + 14 * x*x*x*x;
	return F;
}

ldouble SIF_Factor(properties prop)
{
	return SIF_Factor(prop.a/prop.w);
}

void outputToFile(push_back_state_and_pos &Observer, std::string fileName)
{
	std::ofstream file(fileName);

	for (int i = 0; i < Observer.pos.size(); ++i)
	{
		file << Observer.pos[i] << '\t' << Observer.states[0][i] << '\t' << Observer.states[1][i] << std::endl;
	}
}

ldouble percentChange(ldouble init, ldouble fin)
{
	return ((1-fin/init));
}

Coord getPoints(std::vector<state_type> state, std::vector<Coord> &points, Coord initCoord)
{
	std::vector<state_type>::reverse_iterator stateIter = state.rbegin();
	
	Coord prevCoord(initCoord), tempCoord;
	ldouble angle = 0;
	
	while (stateIter != state.rend())
	{
		angle = (*stateIter)[0];
		//std::cout << angle << '\n';
		tempCoord.x = prevCoord.x + stepSize * sin(angle);
		tempCoord.y = prevCoord.y + stepSize * cos(angle);

		points.push_back(Coord(tempCoord.x, tempCoord.y));

		++stateIter;
		prevCoord = tempCoord;
	}

	return prevCoord;
}