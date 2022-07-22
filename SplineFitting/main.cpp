
#include "spline_curve_fitting.h"
#include "read_write_asc.h"

#include <iostream>

//int main()
//{
//
//	string inFileName = "generated.txt";
//	string outFileName1 = inFileName + "_controls.txt";
//	string outFileName2 = inFileName + "_spline.txt";
//
//	CubicBSplineCurve curve(0.002);
//	SplineCurveFitting scf;
//
//	std::vector<Vector2d> points = generate_points(100);
//
//	scf.apply(points, curve, 28, 50, 0.005, 0.005, 0.0001, SPHERE_INIT);
//
//		CReadWriteAsc::writeAsc( inFileName, points);
//
//	CReadWriteAsc::writeAsc(outFileName1, curve.getControls());
//	CReadWriteAsc::writeAsc(outFileName2, curve.getSamples());
//	CReadWriteAsc::writeAsc(inFileName, points);
//}


int main(int argc, char* argv[])
{
	//char inpf[200], * input;
	//argc--; argv++;					//Skip program name arg

	//if (argc < 1)
	//{
	//	cout << "Input file:" << endl;
	//	cin >> inpf;
	//	input = inpf;
	//}
	//else input = argv[0];

	string inFileName("01.txt");
	string outFileName1 = inFileName + "_controls.txt";
	string outFileName2 = inFileName + "_spline.txt";

	CubicBSplineCurve curve(0.002);
	SplineCurveFitting scf;


	std::vector<Vector2d> points;
	CReadWriteAsc::readAsc(inFileName, points);

	scf.apply(points, curve, 28, 50, 0.005, 0.005, 0.0001, SPHERE_INIT);

	//	CReadWriteAsc::writeAsc( inFileName, points);
	CReadWriteAsc::writeAsc(outFileName1, curve.getControls());
	CReadWriteAsc::writeAsc(outFileName2, curve.getSamples());
}
