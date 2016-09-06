#include <cmath>
#include <stdlib.h>

namespace spline
{
	void DersBasisFuns(int i, double u, int p, double* U, int n, double** & ders) {	//checked
		ders = new double*[n + 1];
		for (int m = 0; m < n + 1; m++)
			ders[m] = new double[p + 1];
		double** ndu;
		ndu = new double*[p + 1];
		for (int m = 0; m < p + 1; m++)
			ndu[m] = new double[p + 1];
		double** a;
		a = new double*[2];
		for (int m = 0; m < 2; m++)
			a[m] = new double[p + 1];
		ndu[0][0] = 1;
		double* left;
		left = new double[p + 1];
		double* right;
		right = new double[p + 1];
		for (int j = 1; j <= p; j++) {
			left[j] = u - U[i + 1 - j];
			right[j] = U[i + j] - u;
			double saved = 0;
			for (int r = 0; r < j; r++) {
				ndu[j][r] = right[r + 1] + left[j - r];
				double temp = ndu[r][j - 1] / ndu[j][r];
				ndu[r][j] = saved + right[r + 1] * temp;
				saved = left[j - r] * temp;
			}
			ndu[j][j] = saved;
		}
		for (int j = 0; j <= p; j++)
			ders[0][j] = ndu[j][p];
		for (int r = 0; r <= p; r++) {
			int s1 = 0, s2 = 1;
			a[0][0] = 1;
			for (int k = 1; k <= n; k++) {
				double d = 0;
				int rk = r - k, pk = p - k;
				if (r >= k) {
					a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
					d = a[s2][0] * ndu[rk][pk];
				}
				int j1, j2;
				if (rk >= -1)
					j1 = 1;
				else
					j1 = -rk;
				if (r - 1 <= pk)
					j2 = k - 1;
				else
					j2 = p - r;
				for (int j = j1; j <= j2; j++) {
					a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
					d += a[s2][j] * ndu[rk + j][pk];
				}
				if (r <= pk) {
					a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
					d += a[s2][k] * ndu[r][pk];
				}
				ders[k][r] = d;
				int j = s1;
				s1 = s2;
				s2 = j;
			}
		}
		int r = p;
		for (int k = 1; k <= n; k++) {
			for (int j = 0; j <= p; j++)
				ders[k][j] *= r;
			r *= (p - k);
		}
		delete[] left;
		delete[] right;
		for (int k = 0; k < 2; k++)
			delete a[k];
		delete[] a;
		for (int k = 0; k < p + 1; k++)
			delete ndu[k];
		delete[] ndu;
	}
	int Findspan(int m, int p, double* U, double u) {					//checked
																		/* This function determines the knot span.*/
		int n = m - p - 1;
		if (u >= U[n + 1])
			return n;
		if (u <= U[p])
			return p;
		int low = p, high = n + 1;
		int mid = (low + high) / 2;
		while (u < U[mid] || u >= U[mid + 1]) {
			if (u < U[mid])
				high = mid;
			else
				low = mid;
			mid = (low + high) / 2;
		}
		return mid;
		/*Test knot={ 0,0,0,1,2,4,4,5,6,6,6 }, u=0 return 2, 2.5 return 4, 4 return 6,3.9999 return 4,
		6 return 7.
		*/
	}
}
