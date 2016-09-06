#define EIGEN_NO_DEBUG
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include "spline.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
using namespace std;
using namespace spline;
using namespace Eigen;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::SparseMatrix<double> SpMat;
MatrixXd globalExtraction(double* knot, int p, int m) {
	vector<double> knots;
	MatrixXd CC = MatrixXd::Identity(m - p, m - p);
	for (int i = 0; i < m + 1; i++)
		knots.push_back(knot[i]);
	int delay = 0, cont = 0, conts = 0;
	for (int i = 0; i < m + 1 - 2 * (p + 1); i++) {
		if (knot[i + p + 1] == knot[i + p + 2]) {
			cont++;
			delay++;
			continue;
		}
		int k = (i - delay + 1) * (p + 1);
		for (int j = 0; j < p - cont; j++) {
			MatrixXd C = MatrixXd::Zero(m - p + (i - delay) * p - conts + j + 1,
			                            m - p + (i - delay) * p - conts + j);
			for (int ii = 0; ii < k - p + 1; ii++)
				C(ii, ii) = 1;
			for (int ii = k - p + 1; ii <= k; ii++) {
				C(ii, ii) = (knot[i + p + 1] - knots[ii]) / (knots[ii + p] - knots[ii]);
				C(ii, ii - 1) = 1 - (knot[i + p + 1] - knots[ii]) / (knots[ii + p] - knots[ii]);
			}
			for (int ii = k + 1; ii < m - p + (i - delay) * p - conts + j + 1; ii++)
				C(ii, ii - 1) = 1;
			CC *= C.transpose();
			knots.insert(knots.begin() + (i - delay) * (p + 1) + p + 1, knot[i + p + 1]);
		}
		conts = delay;
		cont = 0;
	}
	return CC;
}
void GenerateKnot(int order, int refine, int insert, double* insert_knot,
                  double* & knot, int& m);
void CombineKnot(double* knot1, int m1, double* knot2, int m2, double* & knot,
                 int& m);
void Geometry(double xi, double eta, double& pxpxi, double& pxpeta,
              double& pypxi, double& pypeta);
void CompToPhy(double xi, double eta, double& x, double& y);
void AnalyticalDeformation(double x, double y, double L, double D, double I,
                           double P, double E,
                           double nu,
                           Vector2d& deformation);
void AnalyticalStress(double x, double y, double L, double D, double I,
                      double P, Vector3d& stress);
int main() {
	double kkkkk[] = {0, 0, 0, 0, 0, .5, 1, 1, 1, 1, 1};
	cout << globalExtraction(kkkkk, 4, 10) << endl;
	double P = -1000;
	double E = 1000;
	double nu = 0.3;
	double L = 4;
	double D = 2;
	double I = D * D * D / 12;
	double* gaussian = x6;
	double* weight = w6;
	int gaussian_points = 6;
	MatrixXd DD(3, 3);
	DD << 1 - nu, nu, 0, nu, 1 - nu, 0, 0, 0, (1.0 - 2 * nu) / 2;
	DD *= E / (1 + nu) / (1 - 2 * nu);
	int order = 4;
	int refine = 4;
	int m_x_patch1, m_y_patch1;
	int m_x_patch2, m_y_patch2;
	int p_x = order, p_y = order;
	double* knots_x_patch1, *knots_y_patch1;
	double* knots_x_patch2, *knots_y_patch2;
	double insertion_patch1[] = {.5};
	double insertion_patch2[] = {1.0/4, .5, 3.0/4, };
	GenerateKnot(p_x, refine, 1, insertion_patch1, knots_x_patch1, m_x_patch1);
	GenerateKnot(p_y, refine, 3, insertion_patch2, knots_y_patch1, m_y_patch1);
	GenerateKnot(p_x, refine, 3, insertion_patch2, knots_x_patch2, m_x_patch2);
	GenerateKnot(p_y, refine, 3, insertion_patch2, knots_y_patch2, m_y_patch2);
	double* knots_y_coupling;
	int m_y_coupling;
	CombineKnot(knots_y_patch1, m_y_patch1, knots_y_patch2, m_y_patch2,
	            knots_y_coupling,
	            m_y_coupling);
	int dof_x_patch1 = m_x_patch1 - p_x, dof_y_patch1 = m_y_patch1 - p_y,
	    dof_x_patch2 = m_x_patch2 - p_x,
	    dof_y_patch2 = m_y_patch2 - p_y;
	int dof_patch1 = dof_x_patch1 * dof_y_patch1,
	    dof_patch2 = dof_x_patch2 * dof_y_patch2;
	int elements_x_patch1 = m_x_patch1 - 2 * p_x,
	    elements_y_patch1 = m_y_patch1 - 2 * p_y,
	    elements_x_patch2 = m_x_patch2 - 2 * p_x,
	    elements_y_patch2 = m_y_patch2 - 2 * p_y;
	int elements_y_coupling = m_y_coupling - 2 * p_y;
	MatrixXd M = MatrixXd::Zero(dof_y_patch2, dof_y_patch2),
	         N2N1 = MatrixXd::Zero(dof_y_patch2, dof_y_patch1);
	for (int ii_y = 0; ii_y < elements_y_coupling; ii_y++) {
		double J_y = (knots_y_coupling[ii_y + p_y + 1] - knots_y_coupling[ii_y + p_y]) /
		             2;
		double Middle_y = (knots_y_coupling[ii_y + p_y + 1] + knots_y_coupling[ii_y +
		                   p_y]) / 2;
		int i_y_patch1 = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
		int i_y_patch2 = Findspan(m_y_patch2, p_y, knots_y_patch2, Middle_y);
		for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
			double eta = Middle_y + J_y * gaussian[jj_y];
			double** ders_y_patch1, **ders_y_patch2;
			DersBasisFuns(i_y_patch1, eta, p_y, knots_y_patch1, 0, ders_y_patch1);
			DersBasisFuns(i_y_patch2, eta, p_y, knots_y_patch2, 0, ders_y_patch2);
			VectorXd Neta_patch1 = VectorXd::Zero(dof_y_patch1),
			         Neta_patch2 = VectorXd::Zero(dof_y_patch2);
			for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
				Neta_patch1(i_y_patch1 - p_y + kk_y) = ders_y_patch1[0][kk_y];
				Neta_patch2(i_y_patch2 - p_y + kk_y) = ders_y_patch2[0][kk_y];
			}
			for (int k = 0; k < 1; k++)
				delete ders_y_patch1[k];
			delete[] ders_y_patch1;
			for (int k = 0; k < 1; k++)
				delete ders_y_patch2[k];
			delete[] ders_y_patch2;
			M += weight[jj_y] * Neta_patch2 * Neta_patch2.transpose() * J_y;;
			N2N1 += weight[jj_y] * Neta_patch2 * Neta_patch1.transpose() * J_y;
		}
	}
	MatrixXd MN2N1 = M.fullPivHouseholderQr().solve(N2N1);
	MatrixXd Iden = MatrixXd::Identity(2, 2);
	MatrixXd couplingMatrix = MatrixXd::Zero(2 * dof_patch2,
	                          2 * (dof_patch2 - dof_y_patch2 + dof_y_patch1));
	couplingMatrix.block(0, 0, 2 * dof_y_patch2,
	                     2 * dof_y_patch1) = kroneckerProduct(MN2N1, Iden);
	couplingMatrix.block(2 * dof_y_patch2, 2 * dof_y_patch1,
	                     2 * (dof_patch2 - dof_y_patch2),
	                     2 * (dof_patch2 - dof_y_patch2)) = MatrixXd::Identity(2 *
	                             (dof_patch2 - dof_y_patch2), 2 * (dof_patch2 - dof_y_patch2));
	MatrixXd K_patch1 = MatrixXd::Zero(2 * dof_patch1, 2 * dof_patch1);
	for (int ii_x = 0; ii_x < elements_x_patch1; ii_x++) {
		double J_x = (knots_x_patch1[ii_x + p_x + 1] - knots_x_patch1[ii_x + p_x]) / 2;
		double Middle_x = (knots_x_patch1[ii_x + p_x + 1] + knots_x_patch1[ii_x +
		                   p_x]) / 2;
		int i_x = Findspan(m_x_patch1, p_x, knots_x_patch1, Middle_x);
		for (int ii_y = 0; ii_y < elements_y_patch1; ii_y++) {
			double J_y = (knots_y_patch1[ii_y + p_y + 1] - knots_y_patch1[ii_y + p_y]) / 2;
			double Middle_y = (knots_y_patch1[ii_y + p_y + 1] + knots_y_patch1[ii_y +
			                   p_y]) / 2;
			int i_y = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
			for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
				for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
					double xi = Middle_x + J_x * gaussian[jj_x];
					double eta = Middle_y + J_y * gaussian[jj_y];
					double** ders_x, **ders_y;
					DersBasisFuns(i_x, xi, p_x, knots_x_patch1, 1, ders_x);
					DersBasisFuns(i_y, eta, p_y, knots_y_patch1, 1, ders_y);
					VectorXd Nxi(p_x + 1), Nxi_xi(p_x + 1), Neta(p_y + 1), Neta_eta(p_y + 1);
					for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
						Nxi(kk_x) = ders_x[0][kk_x];
						Nxi_xi(kk_x) = ders_x[1][kk_x];
					}
					for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
						Neta(kk_y) = ders_y[0][kk_y];
						Neta_eta(kk_y) = ders_y[1][kk_y];
					}
					for (int k = 0; k < 2; k++)
						delete ders_x[k];
					delete[] ders_x;
					for (int k = 0; k < 2; k++)
						delete ders_y[k];
					delete[] ders_y;
					VectorXd Nxi_xiNeta, NxiNeta_eta;
					Nxi_xiNeta = kroneckerProduct(Nxi_xi, Neta);
					NxiNeta_eta = kroneckerProduct(Nxi, Neta_eta);
					double pxpxi, pxpeta, pypxi, pypeta;
					Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
					double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
					VectorXd Nx_xNy, NxNy_y;
					Nx_xNy = 1.0 / Jacobian * (Nxi_xiNeta * pypeta - NxiNeta_eta * pypxi);
					NxNy_y = 1.0 / Jacobian * (-Nxi_xiNeta * pxpeta + NxiNeta_eta * pxpxi);
					for (int kkx = 0; kkx < (p_x + 1) * (p_y + 1); kkx++) {
						for (int kky = 0; kky < (p_x + 1) * (p_y + 1); kky++) {
							MatrixXd Bx(3, 2);
							MatrixXd By(3, 2);
							Bx << Nx_xNy(kkx), 0, 0, NxNy_y(kkx), NxNy_y(kkx), Nx_xNy(kkx);
							By << Nx_xNy(kky), 0, 0, NxNy_y(kky), NxNy_y(kky), Nx_xNy(kky);
							K_patch1.block(2 * ((m_y_patch1 - p_y) * (kkx / (p_y + 1) + i_x - p_x) + kkx %
							                    (p_y + 1) + i_y - p_y), 2 * ((m_y_patch1 - p_y) * (kky /
							                            (p_y + 1) + i_x - p_x) + kky
							                            % (p_y + 1) + i_y - p_y), 2,
							               2) += weight[jj_x] * weight[jj_y] * Jacobian * Bx.transpose() * DD * By * J_x *
							                     J_y;
						}
					}
				}
			}
		}
	}
	MatrixXd K_patch2 = MatrixXd::Zero(2 * dof_patch2, 2 * dof_patch2);
	for (int ii_x = 0; ii_x < elements_x_patch2; ii_x++) {
		double J_x = (knots_x_patch2[ii_x + p_x + 1] - knots_x_patch2[ii_x + p_x]) / 2;
		double Middle_x = (knots_x_patch2[ii_x + p_x + 1] + knots_x_patch2[ii_x +
		                   p_x]) / 2;
		int i_x = Findspan(m_x_patch2, p_x, knots_x_patch2, Middle_x);
		for (int ii_y = 0; ii_y < elements_y_patch2; ii_y++) {
			double J_y = (knots_y_patch2[ii_y + p_y + 1] - knots_y_patch2[ii_y + p_y]) / 2;
			double Middle_y = (knots_y_patch2[ii_y + p_y + 1] + knots_y_patch2[ii_y +
			                   p_y]) / 2;
			int i_y = Findspan(m_y_patch2, p_y, knots_y_patch2, Middle_y);
			for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
				for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
					double xi = Middle_x + J_x * gaussian[jj_x];
					double eta = Middle_y + J_y * gaussian[jj_y];
					double** ders_x, **ders_y;
					DersBasisFuns(i_x, xi, p_x, knots_x_patch2, 1, ders_x);
					DersBasisFuns(i_y, eta, p_y, knots_y_patch2, 1, ders_y);
					VectorXd Nxi(p_x + 1), Nxi_xi(p_x + 1), Neta(p_y + 1), Neta_eta(p_y + 1);
					for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
						Nxi(kk_x) = ders_x[0][kk_x];
						Nxi_xi(kk_x) = ders_x[1][kk_x];
					}
					for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
						Neta(kk_y) = ders_y[0][kk_y];
						Neta_eta(kk_y) = ders_y[1][kk_y];
					}
					for (int k = 0; k < 2; k++)
						delete ders_x[k];
					delete[] ders_x;
					for (int k = 0; k < 2; k++)
						delete ders_y[k];
					delete[] ders_y;
					VectorXd Nxi_xiNeta, NxiNeta_eta;
					Nxi_xiNeta = kroneckerProduct(Nxi_xi, Neta);
					NxiNeta_eta = kroneckerProduct(Nxi, Neta_eta);
					double pxpxi, pxpeta, pypxi, pypeta;
					Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
					double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
					VectorXd Nx_xNy, NxNy_y;
					Nx_xNy = 1.0 / Jacobian * (Nxi_xiNeta * pypeta - NxiNeta_eta * pypxi);
					NxNy_y = 1.0 / Jacobian * (-Nxi_xiNeta * pxpeta + NxiNeta_eta * pxpxi);
					for (int kkx = 0; kkx < (p_x + 1) * (p_y + 1); kkx++) {
						for (int kky = 0; kky < (p_x + 1) * (p_y + 1); kky++) {
							MatrixXd Bx(3, 2);
							MatrixXd By(3, 2);
							Bx << Nx_xNy(kkx), 0, 0, NxNy_y(kkx), NxNy_y(kkx), Nx_xNy(kkx);
							By << Nx_xNy(kky), 0, 0, NxNy_y(kky), NxNy_y(kky), Nx_xNy(kky);
							K_patch2.block(2 * ((m_y_patch2 - p_y) * (kkx / (p_y + 1) + i_x - p_x) + kkx %
							                    (p_y + 1) + i_y - p_y), 2 * ((m_y_patch2 - p_y) * (kky /
							                            (p_y + 1) + i_x - p_x) + kky
							                            % (p_y + 1) + i_y - p_y), 2,
							               2) += weight[jj_x] * weight[jj_y] * Jacobian * Bx.transpose() * DD * By * J_x *
							                     J_y;
						}
					}
				}
			}
		}
	}
	MatrixXd K = MatrixXd::Zero(2 * (dof_patch1 + dof_patch2 -
	                                 dof_y_patch2), 2 * (dof_patch1 + dof_patch2 - dof_y_patch2));
	K.block(0, 0, 2 * dof_patch1, 2 * dof_patch1) = K_patch1;
	K.block(2 * (dof_patch1 - dof_y_patch1), 2 * (dof_patch1 - dof_y_patch1),
	        2 * (dof_patch2 - dof_y_patch2 + dof_y_patch1),
	        2 * (dof_patch2 - dof_y_patch2 + dof_y_patch1)) += couplingMatrix.transpose()
	                * K_patch2 * couplingMatrix;


	VectorXd FT = VectorXd::Zero(2 * (dof_patch1 + dof_patch2 -
	                                  dof_y_patch2));
	for (int ii_y = 0; ii_y < elements_y_patch2; ii_y++) {
		double J_y = (knots_y_patch2[ii_y + p_y + 1] - knots_y_patch2[ii_y + p_y]) / 2;
		double Middle_y = (knots_y_patch2[ii_y + p_y + 1] + knots_y_patch2[ii_y +
		                   p_y]) / 2;
		int i_y = Findspan(m_y_patch2, p_y, knots_y_patch2, Middle_y);
		for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
			double xi = 0.9999999999999999999999;
			double eta = Middle_y + J_y * gaussian[jj_y];
			double** ders_y;
			DersBasisFuns(i_y, eta, p_y, knots_y_patch2, 0, ders_y);
			VectorXd Nxi = VectorXd::Zero(dof_x_patch2),
			         Neta = VectorXd::Zero(dof_y_patch2);
			for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
				Neta(i_y - p_y + kk_y) = ders_y[0][kk_y];
			}
			Nxi(dof_x_patch2 - 1) = 1;
			for (int k = 0; k < 1; k++)
				delete ders_y[k];
			delete[] ders_y;
			VectorXd NxiNetaY, NxiNeta;
			VectorXd Y(2);
			Y(1) = 1, Y(0) = 0;
			NxiNeta = kroneckerProduct(Nxi, Neta);
			NxiNetaY = kroneckerProduct(NxiNeta, Y);
			double pxpxi, pxpeta, pypxi, pypeta;
			Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
			double x, y;
			CompToPhy(xi, eta, x, y);
			Vector3d stress;
			AnalyticalStress(x, y, L, D, I, P, stress);
			FT.segment(2 * (dof_patch1 - dof_y_patch1),
			           2 * (dof_patch2 - dof_y_patch2 + dof_y_patch1)) += weight[jj_y] * (stress(
			                       2) * NxiNetaY) * J_y * pypeta;
		}
	}

	
	MatrixXd K_compute = K.block(2 * dof_y_patch1, 2 * dof_y_patch1,
	                             2 * (dof_patch1 + dof_patch2 -
	                                  dof_y_patch2 - dof_y_patch1), 2 * (dof_patch1 + dof_patch2 -
	                                          dof_y_patch2 - dof_y_patch1));
	VectorXd F_compute = FT.segment(2 * dof_y_patch1, 2 * (dof_patch1 + dof_patch2 -
	                                dof_y_patch2 - dof_y_patch1));
	VectorXd U_result = K_compute.partialPivLu().solve(F_compute);
	cout<<U_result<<endl;
	return 0;
}
void GenerateKnot(int order, int refine, int insert, double* insert_knot,
                  double* & knot, int& m) {
	vector<double> Knot;
	for (int i = 0; i <= order; i++) {
		Knot.push_back(0);
	}
	for (int i = 0; i < insert; i++) {
		Knot.push_back(insert_knot[i]);
	}
	for (int i = 0; i <= order; i++) {
		Knot.push_back(1);
	}
	for (int i = 0; i < refine; i++) {
		for (int j = 0; j < Knot.end() - Knot.begin() - 1; j++) {
			double insertion;
			if (Knot[j] != Knot[j + 1]) {
				insertion = (Knot[j] + Knot[j + 1]) / 2;
				Knot.insert(Knot.begin() + j + 1, insertion);
				j++;
			}
		}
	}
	m = Knot.end() - Knot.begin() - 1;
	knot = new double[m + 1];
	for (int i = 0; i < m + 1; i++)
		knot[i] = Knot[i];
}
void CombineKnot(double* knot1, int m1, double* knot2, int m2, double* & knot,
                 int& m) {
	vector<double> Knot;
	int i1 = 0, i2 = 0;
	while (i1 <= m1) {
		if (knot1[i1] == knot2[i2]) {
			Knot.push_back(knot1[i1]);
			i1++, i2++;
		} else if (knot1[i1] < knot2[i2]) {
			Knot.push_back(knot1[i1]);
			i1++;
		} else {
			Knot.push_back(knot2[i2]);
			i2++;
		}
	}

	m = Knot.end() - Knot.begin() - 1;
	knot = new double[m + 1];
	for (int i = 0; i < m + 1; i++)
		knot[i] = Knot[i];
}
void AnalyticalDeformation(double x, double y, double L, double D, double I,
                           double P, double E,
                           double nu,
                           Vector2d& deformation) {
	double E_bar = E / (1 - pow(nu, 2));
	double nu_bar = nu / (1 - nu);
	deformation(0) = -P * y / 6 / E_bar / I * ((6 * L - 3 * x) * x +
	                 (2 + nu_bar) * (pow(y,
	                                     2) - pow(D / 2, 2)));
	deformation(1) = P / 6 / E_bar / I * (3 * nu_bar * pow(y,
	                                      2) * (L - x) + (4 + 5 * nu_bar) * D * D * x / 4 + (3 * L - x) * x * x);
}
void AnalyticalStress(double x, double y, double L, double D, double I,
                      double P,
                      Vector3d& stress) {
	stress(0) = -P * (L - x) * y / I;
	stress(1) = 0;
	stress(2) = P / 2 / I * (pow(D / 2, 2) - pow(y, 2));
}
void Geometry(double xi, double eta, double& pxpxi, double& pxpeta,
              double& pypxi, double& pypeta) {
	pxpxi = 2;
	pxpeta = 0;
	pypxi = 0;
	pypeta = 2;
}
void CompToPhy(double xi, double eta, double& x, double& y) {
	y = (eta - 0.5) * 2;
	x = xi * 4;
}
