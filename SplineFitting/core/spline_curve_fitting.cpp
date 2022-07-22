#include "spline_curve_fitting.h"
#include <eigen/SVD>


#define TEMP_PI       3.14159265358979323846


void SplineCurveFitting::initControlPoint(const vector<Vector2d>& datapoints,
	vector<Vector2d>& controlPs,
	int controlNum,
	EInitType initType)
{
	// compute the initial 12 control points
	controlPs.clear();
	int perNum = controlNum / 4;

	// ���������һif��������ʲô��������Ϊ��ʵ�ֵĺܺã�������Ĭ��Ҳ�õĲ������ʵ�֣�
	if (initType == RECT_INIT)
	{
		Vector2d v1 = datapoints[0];
		Vector2d v2 = datapoints[0];

		// compute bounding box for datapoints
		for (unsigned int i = 0; i != datapoints.size(); ++i) {
			Vector2d v = datapoints[i];
			if (v1.x() > v.x())  v1.x() = v.x();
			if (v1.y() > v.y())  v1.y() = v.y();
			if (v2.x() < v.x())  v2.x() = v.x();
			if (v2.y() < v.y())  v2.y() = v.y();
		}


		Vector2d dir = (v2 - v1) * 0.5;

		// center of the bounding box
		Vector2d cen = v1 + dir;

		// extend the bounding box a bit further to potentially aid floating point inprecise computation
		v1 = cen - 1.05 * dir;
		v2 = cen + 1.05 * dir;

		vector<Vector2d> rets;
		rets.push_back(v1);
		rets.push_back(Vector2d(v1.x(), v2.y()));
		rets.push_back(v2);
		rets.push_back(Vector2d(v2.x(), v1.y()));
		rets.push_back(v1);


		for (int i = 0; i < 4; i++)
		{
			Vector2d p1 = rets[i];
			Vector2d p2 = rets[i + 1];
			for (int j = 0; j < perNum; j++) {
				controlPs.push_back(p1 + (p2 - p1) * j / (double)(perNum));
			}
		}
	}
	else
	{
		// control points��һ��bounding circle��ȡ��
		// ��bounding circle�ϲ���ָ����Ŀ�ĵ�����Ϊcontrol points.
		// ���⣺bounding circle�Ǳպϵģ��Ƿ�Ҳʹ�����õ���bsplineҲ�Ǳպϵģ�which is not sth i want��
		// compute bounding circle.
		Vector2d cen(0, 0);
		for (size_t i = 0; i != datapoints.size(); ++i)
		{
			cen += datapoints[i];
		}
		cen /= datapoints.size();

		double radius = 0;
		for (size_t i = 0; i != datapoints.size(); ++i)
		{
			double len = (datapoints[i] - cen).norm();
			if (radius < len)
				radius = len;
		}

		double theta = (2 * TEMP_PI) / controlNum;
		for (size_t i = 0; i != controlNum; ++i)
		{
			Vector2d pos = cen + radius * Vector2d(std::cos(theta * i), std::sin(theta * i));
			controlPs.push_back(pos);
		}
	}

}

double SplineCurveFitting::apply(
	const vector<Vector2d>& datapoints,
	CubicBSplineCurve& curve,
	int controlNum /* = 28 */,
	int maxIterNum  /*= 30 */,
	double alpha /* = 0.002 */,
	double gama /* = 0.002 */,
	double eplison /* = 0.0001 */,
	EInitType initType /* =SPHERE_INIT */)
{
	// ��control point �ĸ������¿���4�ı���
	controlNum = controlNum / 4 * 4;

	// initialize the cube B-spline
	CubicBSplineCurve* spline = &curve;
	vector<Vector2d> controlPs;

	// ok
	initControlPoint(datapoints, controlPs, controlNum, initType);


	spline->setNewControl(controlPs);

	// update the control point
	// compute P"(t)
	MatrixXd pm = spline->getSIntegralSq(); // ����������������ĺ�������ʲô
	MatrixXd sm = spline->getFIntegralSq();
	// end test

	// find the (approximate, not exact) foot prints (for all datapoints), and return the total corresponding error
	std::vector< std::pair<int, double> > parameters;
	double fsd = spline->findFootPrint(datapoints, parameters);
	int iterNum = 0;
	while (fsd > eplison && iterNum < maxIterNum)
	{
		MatrixXd ehm(2 * controlNum, 2 * controlNum);
		VectorXd ehv(2 * controlNum);

		ehm.setZero();
		ehv.setZero();

		// compute h(D)
		// ForEach: �Ծ�ÿ��datapoint�������sample point��
		for (int i = 0; i < (int)parameters.size(); i++)
		{
			// 	if (labels[i] == false)
			// 		continue;

			// compute d, rho, Tkv, Nkv
			double curvature = spline->getCurvature(parameters[i]);
			// ���ʰ뾶����Ϊ����curvature��0���԰ѳ�ʼ������������
			double rho = 10e+6;

			// TODO: ��Ū�������getPos����ô����
			// �䷵��parameters[i] ��Ӧ�� b-spline �ϵĲ����㣨��Ȼ��û�з���positions_��������������һ�飩
			Vector2d neip = spline->getPos(parameters[i]);

			Vector2d Tkv = spline->getTangent(parameters[i]);
			Vector2d Nkv = spline->getNormal(parameters[i]);

			// datapoint �� sample point ֮��ľ���
			double d = (datapoints[i] - neip).norm();
			Vector2d curvature_center(0.0, 0.0);
			bool is_datapoints_neip_on_the_same_side = true;


			if (curvature != 0.0f)
			{
				rho = 1 / curvature;
				curvature_center = spline->getCurvCenter(parameters[i]);
				double ddd = (curvature_center - neip).norm();
				is_datapoints_neip_on_the_same_side = spline->checkSameSide(curvature_center, datapoints[i], neip);
			}

			// �����е���ء���
			// ���մ����￪ʼŪ

			VectorXd coffv = spline->getCoffe(parameters[i]);
			MatrixXd tempcoffm1(controlNum, 1);
			for (int ij = 0; ij < controlNum; ij++)
				tempcoffm1(ij, 0) = coffv[ij];
			MatrixXd tempcoffm = tempcoffm1 * (tempcoffm1.transpose());

			// update the matrix
			double fxx = Tkv.x() * Tkv.x();
			double fyy = Tkv.y() * Tkv.y();
			double fxy = Tkv.x() * Tkv.y();

			Vector2d oldp = neip - datapoints[i];

			if (!is_datapoints_neip_on_the_same_side)
			{
				d = -d;
				VectorXd tempv1 = (coffv) * (d / (d - rho)) * (fxx * datapoints[i].x() + fxy * datapoints[i].y());
				VectorXd tempv2 = (coffv) * (d / (d - rho)) * (fyy * datapoints[i].y() + fxy * datapoints[i].x());
				for (int i2 = 0; i2 < controlNum; i2++)
				{
					for (int j = 0; j < controlNum; j++)
					{
						double fp = (d / (d - rho)) * tempcoffm(i2, j);
						ehm(i2, j) += fxx * fp;
						ehm(i2, j + controlNum) += fxy * fp;
						ehm(i2 + controlNum, j) += fxy * fp;
						ehm(i2 + controlNum, j + controlNum) += fyy * fp;
					}
					ehv[i2] += tempv1[i2];
					ehv[i2 + controlNum] += tempv2[i2];
				}
			}
			fxx = Nkv.x() * Nkv.x();
			fyy = Nkv.y() * Nkv.y();
			fxy = Nkv.x() * Nkv.y();
			VectorXd tempv1 = (coffv) * (fxx * datapoints[i].x() + fxy * datapoints[i].y());
			VectorXd tempv2 = (coffv) * (fyy * datapoints[i].y() + fxy * datapoints[i].x());
			for (int i2 = 0; i2 < controlNum; i2++)
			{
				for (int j = 0; j < controlNum; j++)
				{
					double fp = tempcoffm(i2, j);
					ehm(i2, j) += fxx * fp;
					ehm(i2, j + controlNum) += fxy * fp;
					ehm(i2 + controlNum, j) += fxy * fp;
					ehm(i2 + controlNum, j + controlNum) += fyy * fp;
				}
				ehv[i2] += tempv1[i2];
				ehv[i2 + controlNum] += tempv2[i2];
			}
		}


		// check if ehm, ehv right
		//solve the function
		MatrixXd fm = ehm * 0.5 + pm * alpha + sm * gama;
		VectorXd ehv2 = ehv * 0.5;
		JacobiSVD<MatrixXd> svd(fm, ComputeThinU | ComputeThinV);
		VectorXd resultxy = svd.solve(ehv2);

		// update the curve
		for (int i = 0; i < controlNum; i++)
			controlPs[i] = Vector2d(resultxy[i], resultxy[i + controlNum]);
		spline->setNewControl(controlPs);
		++iterNum;

		fsd = spline->findFootPrint(datapoints, parameters);
	}

	return fsd;
}