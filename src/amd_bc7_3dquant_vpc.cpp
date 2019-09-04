//===============================================================================
// Copyright (c) 2007-2016  Advanced Micro Devices, Inc. All rights reserved.
// Copyright (c) 2004-2006 ATI Technologies Inc.
//===============================================================================
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

#include "al2o3_platform/platform.h"
#include "amd_bc7_3dquant_vpc.hpp"

#define EPSILON        0.000001
#undef MAX_TRY
#define MAX_TRY        20
#undef DBL_MAX_EXP
#define DBL_MAX_EXP 1024
#undef  TRACE

#define MAX_TRACE    250000

struct TRACE
{
	int  k;
	float d;
};

static int trcnts[MAX_CLUSTERS][MAX_ENTRIES_QUANT_TRACE];

#define USE_TRACE_WITH_DYNAMIC_MEM

#ifdef USE_TRACE_WITH_DYNAMIC_MEM
int*        amd_codes[MAX_CLUSTERS][MAX_ENTRIES_QUANT_TRACE] = {};
TRACE*      amd_trs[MAX_CLUSTERS][MAX_ENTRIES_QUANT_TRACE]   = {};
#else
int         amd_codes[MAX_CLUSTERS][MAX_ENTRIES_QUANT_TRACE][MAX_TRACE];
    TRACE       amd_trs[MAX_CLUSTERS][MAX_ENTRIES_QUANT_TRACE][MAX_TRACE];
#endif

static int g_Quant_init = 0;
void traceBuilder (int numEntries, int numClusters,struct TRACE tr [], int code[], int *trcnt );

void Quant_Init(void)
{
	if (g_Quant_init > 0)
	{
		g_Quant_init++;
		return;
	}
	if (amd_codes[0][0])  return;

	for ( int numClusters = 0; numClusters < MAX_CLUSTERS; numClusters++ )
	{
		for ( int numEntries = 0; numEntries < MAX_ENTRIES_QUANT_TRACE; numEntries++ )
		{
#ifdef USE_TRACE_WITH_DYNAMIC_MEM
			amd_codes[ numClusters][ numEntries ]    = new int[ MAX_TRACE ];
			amd_trs[ numClusters ][ numEntries ]    = new TRACE[ MAX_TRACE ];

			ASSERT(amd_codes[ numClusters][ numEntries ]);
			ASSERT(amd_trs[ numClusters ][ numEntries ]);
#endif

			traceBuilder (  numEntries+1,
											numClusters+1,
											amd_trs[numClusters][numEntries],
											amd_codes[numClusters][numEntries],
											trcnts[numClusters]+(numEntries));
		}
	}

	g_Quant_init++;
}

void Quant_DeInit(void)
{
	g_Quant_init--;
	if (g_Quant_init > 1)
	{
		return;
	}
	else
	{
		g_Quant_init = 0; // Reset in case user called Quant_DeInit too many times without matching Quant_Init
		if (amd_codes[0][0] == nullptr)  return;

#ifdef USE_TRACE_WITH_DYNAMIC_MEM

		for (int i = 0; i < MAX_CLUSTERS; i++)
		{
			for (int j = 0; j < MAX_ENTRIES_QUANT_TRACE; j++)
			{
				if (amd_codes[i][j])
				{
					delete[] amd_codes[i][j];
					amd_codes[i][j] = nullptr;
				}
				if (amd_trs[i][j])
				{
					delete[] amd_trs[i][j];
					amd_trs[i][j] = nullptr;
				}
			}
		}

#endif
	}

}

//=========================================================================================

inline int a_compare( const void *arg1, const void *arg2 )
{
	if (((a* )arg1)->d-((a* )arg2)->d > 0 ) return 1;
	if (((a* )arg1)->d-((a* )arg2)->d < 0 ) return -1;
	return 0;
};

//
// We ignore the issue of ordering equal elements here, though it can affect results abit
//
void sortProjection(float projection[MAX_ENTRIES], int order[MAX_ENTRIES], int numEntries)
{
	int i;
	a what[MAX_ENTRIES+MAX_PARTITIONS_TABLE];

	for (i=0; i < numEntries;i++)
		what[what[i].i=i].d = projection[i];

	qsort((void*)&what, numEntries, sizeof(a),a_compare);

	for (i=0; i < numEntries;i++)
		order[i]=what[i].i;
};

void covariance(float data[][DIMENSION], int numEntries, float cov[DIMENSION][DIMENSION])
{
	int i,j,k;

	for(i=0;i<DIMENSION;i++)
		for(j=0;j<=i;j++) {
			cov[i][j]=0;
			for(k=0;k<numEntries;k++)
				cov[i][j]+=data[k][i]*data[k][j];
		}

	for(i=0;i<DIMENSION;i++)
		for(j=i+1;j<DIMENSION;j++)
			cov[i][j] = cov[j][i];
}

void covariance_d(float data[][MAX_DIMENSION_BIG], int numEntries, float cov[MAX_DIMENSION_BIG][MAX_DIMENSION_BIG], int dimension)
{
	int i,j,k;

	for(i=0;i<dimension;i++)
		for(j=0;j<=i;j++)
		{
			cov[i][j]=0;
			for(k=0;k<numEntries;k++)
				cov[i][j]+=data[k][i]*data[k][j];
		}

	for(i=0;i<dimension;i++)
		for(j=i+1;j<dimension;j++)
			cov[i][j] = cov[j][i];
}

void centerInPlace(float data[][DIMENSION], int numEntries, float mean[DIMENSION])
{
	int i,k;

	for(i=0;i<DIMENSION;i++) {
		mean[i]=0;
		for(k=0;k<numEntries;k++)
			mean[i]+=data[k][i];
	}

	if (!numEntries)
		return;

	for(i=0;i<DIMENSION;i++) {
		mean[i]/=(float) numEntries;
		for(k=0;k<numEntries;k++)
			data[k][i]-=mean[i];
	}
}

void centerInPlace_d(float data[][MAX_DIMENSION_BIG], int numEntries, float mean[MAX_DIMENSION_BIG], int dimension)
{
	int i,k;

	for(i=0;i<dimension;i++)
	{
		mean[i]=0;
		for(k=0;k<numEntries;k++)
			mean[i]+=data[k][i];
	}

	if (!numEntries)
		return;

	for(i=0;i<dimension;i++)
	{
		mean[i]/=(float) numEntries;
		for(k=0;k<numEntries;k++)
			data[k][i]-=mean[i];
	}
}

void project(float data[][DIMENSION], int numEntries, float vector[DIMENSION], float projection[MAX_ENTRIES])
{
	// assume that vector is normalized already
	int i,k;

	for(k=0;k<numEntries;k++) {
		projection[k]=0;
		for(i=0;i<DIMENSION;i++) {
			projection[k]+=data[k][i]*vector[i];
		}
	}
}

void project_d(float data[][MAX_DIMENSION_BIG], int numEntries, float vector[MAX_DIMENSION_BIG], float projection[MAX_ENTRIES], int dimension)
{
	// assume that vector is normalized already
	int i,k;

	for(k=0;k<numEntries;k++)
	{
		projection[k]=0;
		for(i=0;i<dimension;i++)
		{
			projection[k]+=data[k][i]*vector[i];
		}
	}
}

void eigenVector(float cov[DIMENSION][DIMENSION], float vector[DIMENSION])
{
	// calculate an eigenvecto corresponding to a biggest eigenvalue
	// will work for non-zero non-negative matricies only

#define EV_ITERATION_NUMBER 20
#define EV_SLACK            2        /* additive for exp base 2)*/


	int i,j,k,l, m, n,p,q;
	float c[2][DIMENSION][DIMENSION];
	float maxDiag;

	for(i=0;i<DIMENSION;i++)
		for(j=0;j<DIMENSION;j++)
			c[0][i][j] =cov[i][j];

	p = (int) floor(log( (DBL_MAX_EXP - EV_SLACK) / ceil (log((float)DIMENSION)/log(2.)) )/log(2.));

	ASSERT(p>0);

	p = p >0 ? p : 1;

	q =  (EV_ITERATION_NUMBER+p-1) / p;


	l=0;

	for(n=0;n<q; n++) {

		maxDiag = 0;

		for(i=0;i<DIMENSION;i++)
			maxDiag = c[l][i][i] > maxDiag ? c[l][i][i] : maxDiag;

		if (maxDiag<=0)
		{
			return;
		}
		ASSERT(maxDiag > 0);

		for(i=0;i<DIMENSION;i++)
			for(j=0;j<DIMENSION;j++)
				c[l][i][j] /=maxDiag;

		for(m=0;m<p;m++) {
			for(i=0;i<DIMENSION;i++)
				for(j=0;j<DIMENSION;j++) {
					c[1-l][i][j]=0;
					for(k=0;k<DIMENSION;k++)
						c[1-l][i][j]+=c[l][i][k]*c[l][k][j];
				}
			l=1-l;
		}
	}

	maxDiag = 0;
	k =0;

	for(i=0;i<DIMENSION;i++) {
		k = c[l][i][i] > maxDiag ? i : k;
		maxDiag = c[l][i][i] > maxDiag ? c[l][i][i] : maxDiag;
	}
	float t;
	t=0;
	for(i=0;i<DIMENSION;i++) {
		t+=c[l][k][i]*c[l][k][i];
		vector[i]=c[l][k][i];
	}
	// normalization is really optional
	t= sqrt(t);
	ASSERT(t>0);
	if (t<=0)
	{
		return;
	}

	for(i=0;i<DIMENSION;i++)
		vector[i]/=t;
}

void eigenVector_d(float cov[MAX_DIMENSION_BIG][MAX_DIMENSION_BIG], float vector[MAX_DIMENSION_BIG], int dimension)
{
	// calculate an eigenvecto corresponding to a biggest eigenvalue
	// will work for non-zero non-negative matricies only

#define EV_ITERATION_NUMBER 20
#define EV_SLACK            2        /* additive for exp base 2)*/


	int i,j,k,l, m, n,p,q;
	float c[2][MAX_DIMENSION_BIG][MAX_DIMENSION_BIG];
	float maxDiag;

	for(i=0;i<dimension;i++)
		for(j=0;j<dimension;j++)
			c[0][i][j] =cov[i][j];

	p = (int) floor(log( (DBL_MAX_EXP - EV_SLACK) / ceil (log((float)dimension)/log(2.)) )/log(2.));

	ASSERT(p>0);

	p = p >0 ? p : 1;

	q =  (EV_ITERATION_NUMBER+p-1) / p;

	l=0;

	for(n=0;n<q; n++)
	{
		maxDiag = 0;

		for(i=0;i<dimension;i++)
			maxDiag = c[l][i][i] > maxDiag ? c[l][i][i] : maxDiag;

		if (maxDiag<=0)
		{
			return;
		}
		ASSERT(maxDiag >0);

		for(i=0;i<dimension;i++)
			for(j=0;j<dimension;j++)
				c[l][i][j] /=maxDiag;

		for(m=0;m<p;m++) {
			for(i=0;i<dimension;i++)
				for(j=0;j<dimension;j++) {
					float temp=0;
					for(k=0;k<dimension;k++)
					{
						// Notes:
						// This is the most consuming portion of the code and needs optimizing for perfromance
						temp += c[l][i][k]*c[l][k][j];
					}
					c[1-l][i][j]=temp;
				}
			l=1-l;
		}
	}

	maxDiag = 0;
	k =0;

	for(i=0;i<dimension;i++)
	{
		k = c[l][i][i] > maxDiag ? i : k;
		maxDiag = c[l][i][i] > maxDiag ? c[l][i][i] : maxDiag;
	}
	float t;
	t=0;
	for(i=0;i<dimension;i++)
	{
		t+=c[l][k][i]*c[l][k][i];
		vector[i]=c[l][k][i];
	}
	// normalization is really optional
	t= sqrt(t);
	ASSERT(t>0);
	if (t<=0)
	{
		return;
	}
	for(i=0;i<dimension;i++)
		vector[i]/=t;
}

float partition2(float data[][DIMENSION], int numEntries,int index[])
{
	int i,j,k;
	float cov[2][DIMENSION][DIMENSION];
	float center[2][DIMENSION];
	float cnt[2] ={0,0};
	float vector[2][DIMENSION];
	float acc=0;

	for(k=0;k<numEntries;k++)
		cnt[index[k]]++;


	for(i=0;i<DIMENSION;i++) {
		center[0][i]=center[1][i]=0;
		for(k=0;k<numEntries;k++)
			center[index[k]][i]+=data[k][i];
	}


	for(i=0;i<DIMENSION;i++)
		for(j=0;j<=i;j++) {
			cov[0][i][j]=cov[1][i][j]=0;
			for(k=0;k<numEntries;k++)
				cov[index[k]][i][j]+=data[k][i]*data[k][j];
		}

	for(i=0;i<DIMENSION;i++)
		for(j=0;j<=i;j++)
			for (k=0;k<2;k++)
				if (cnt[k]!=0)
					cov[k][i][j] -=center[k][i]*center[k][j]/(float)cnt[k];


	for(i=0;i<DIMENSION;i++)
		for(j=i+1;j<DIMENSION;j++)
			for(k=0;k<2;k++)
				cov[k][i][j] = cov[k][j][i];


	for(k=0;k<2;k++)
		eigenVector(cov[k], vector[k]); // assume the returned vector is nomalized



	for(i=0;i<DIMENSION;i++)
		for(k=0;k<2;k++)
			acc+=cov[k][i][i];

	for(i=0;i<DIMENSION;i++)
		for(j=0;j<DIMENSION;j++)
			for(k=0;k<2;k++)
				acc-=cov[k][i][j]*vector[k][i]*vector[k][j];

	return(acc);
}

void quantEven(float data[MAX_ENTRIES][DIMENSION],int numEntries, int numClusters, int index[MAX_ENTRIES])
{
	// Data should be centered, otherwise will not work
	// The running time (number of iteration of the external loop) is
	//   binomial(numEntries+numClusters-2, numClusters-1)
	// First cluster is always used, (without loss of generality)
	// Ramp should be shifted such, that the first element ramp[0] is 0
	int i,k;

	int level;

	float t,s;
	int c =1;

	int cluster[MAX_CLUSTERS];

	int bestCluster[MAX_CLUSTERS];
	// stores the las index for the cluster


	float  dpAcc       [MAX_CLUSTERS][DIMENSION];
	float  index2Acc   [MAX_CLUSTERS];     // for backtraking
	float  indexAcc    [MAX_CLUSTERS];

	float dRamp2[MAX_CLUSTERS];    // first differenses of the (shifted) ramp squared

	float S;

	float nErrorNum=0;   // not the actual error, but some (decreasing) linear functional of it represented
	// as numerator and denominator
	float nErrorDen=1;

	level=1;

	bestCluster[0]=cluster[0]=cluster[1]=numEntries;


	indexAcc[0]=index2Acc[0]=indexAcc[1]=index2Acc[1]=0;

	for(i=0;i<DIMENSION;i++)
		dpAcc[0][i]=dpAcc[1][i]=0;

	S =  1/sqrt((float) numEntries);

	for(i=1;i<MAX_CLUSTERS;i++) {
		dRamp2[i] = 2*i-1;
	}

	level=1;

	do {
		k = --cluster[level-1];

		indexAcc    [level] += S;
		index2Acc   [level] += dRamp2 [level];

		t=0;
		for(i=0;i<DIMENSION;i++) {
			// using scaled ramp instead of non-scaled here effectively scales the data, so
			// the resulting quantisation will be the same, but the error metric value  will be different

			dpAcc [level][i] += data[k][i];
			t += dpAcc[level][i] * dpAcc[level][i];
		}

		if ((cluster[level]!= numEntries || cluster[level-1]!=0) &&
				nErrorNum * (s=index2Acc[level]-indexAcc[level] * indexAcc[level]) < nErrorDen * t) {
			nErrorNum=t;
			nErrorDen=s;
			for(i=0;i<=level;i++)
				bestCluster[i]=cluster[i];
		}
		c++;

		if (level < numClusters - 1 ) {
			// go up
			level++;
			indexAcc [level]=indexAcc [level-1];
			index2Acc[level]=index2Acc[level-1];

			for(i=0;i<DIMENSION;i++)
				dpAcc [level][i]=dpAcc    [level-1][i];

		}
		else while ((level-1) && cluster[level-1]==cluster[level-2])
				level--;

		cluster[level]=numEntries;

	} while (level != 1 || cluster[level-1] != 0);

	for (level=i=0;i< numEntries;i++) {
		while (i==bestCluster[level])
			level++;
		index[i]=level;
	}
}

void quantLineConstr(float data[][DIMENSION], int order[MAX_ENTRIES],int numEntries, int numClusters, int index[MAX_ENTRIES])
{
	// Data should be centered, otherwise will not work
	// The running time (number of iteration of the external loop) is
	//   binomial(numEntries+numClusters-2, numClusters-1)
	// Index just defines which points should be combined in a cluster
	int i,j,k;

	int level;

	float t,s;

	// We need paddingof 0 on -1 index
	int  cluster_[MAX_CLUSTERS+1]={0};
	int *cluster = cluster_+1;

	int bestCluster[MAX_CLUSTERS];
	// stores the las index for the cluster

	float cov[DIMENSION][DIMENSION];
	float dir[DIMENSION];


	float gcAcc[MAX_CLUSTERS][DIMENSION];// Clusters' graviti centers
	float gcSAcc[MAX_CLUSTERS][DIMENSION];// Clusters' graviti centers


	float nError=0;   // not the actual error, but some (decreasing) linear functional of it represented
	// as numerator and denominator

	level=1;


	bestCluster[0]=cluster[0]=cluster[1]=numEntries;


	for(i=0;i<DIMENSION;i++)
		gcAcc[0][i]=gcAcc[1][i]=0;


	level=1;

	do {
		ASSERT(level >0);

		k = order[--cluster[level-1]];

		s=(cluster[level-1]-cluster[level-2]) == 0 ? 0: 1/sqrt( (float) (cluster[level-1]-cluster[level-2])); // see cluster_ decl for
		// cluster[-1] value
		t=1/sqrt((float) (numEntries-cluster[level-1]));

		for(i=0;i<DIMENSION;i++) {
			gcAcc[level  ][i] += data[k][i];
			gcAcc[level-1][i] -= data[k][i];

			gcSAcc[level-1][i] = gcAcc[level-1][i] * s;
			gcSAcc[level  ][i] = gcAcc[level  ][i] * t;
		}
		covariance(gcSAcc, level+1,  cov);
		eigenVector(cov, dir);
		// assume the vector is normalized here
		t=0;
		for(i=0;i<DIMENSION;i++)
			for(j=0;j<DIMENSION;j++)
				t+= cov[i][j]*dir[i]*dir[j];

		if (t>nError) {

			nError=t;

			for(i=0;i<=level;i++)
				bestCluster[i]=cluster[i];
		}

		if (level < numClusters - 1 ) {
			// go up
			level++;
			for(i=0;i<DIMENSION;i++)
				gcAcc [level][i]=0;

		}
		else while ((level-1) && cluster[level-1]==cluster[level-2]) {
				level--;
				for(i=0;i<DIMENSION;i++)
					gcAcc [level][i]+=gcAcc [level+1][i];
			}


		cluster[level]=numEntries;

	} while (level != 1 || cluster[level-1] != 1);

	for (level=i=0;i< numEntries;i++) {
		while (i==bestCluster[level])
			level++;
		index[order[i]]=level;
	}
}

float totalError(float data[MAX_ENTRIES][DIMENSION],float data2[MAX_ENTRIES][DIMENSION],int numEntries)
{
	int i,j;
	float t=0;
	for (i=0;i<numEntries;i++)
		for (j=0;j<DIMENSION;j++)
			t+= (data[i][j]-data2[i][j])*(data[i][j]-data2[i][j]);
	return t;
};

float totalError_d(float data[MAX_ENTRIES][MAX_DIMENSION_BIG],float data2[MAX_ENTRIES][MAX_DIMENSION_BIG],int numEntries, int dimension)
{
	int i,j;
	float t=0;
	for (i=0;i<numEntries;i++)
		for (j=0;j<dimension;j++)
			t+= (data[i][j]-data2[i][j])*(data[i][j]-data2[i][j]);

	return t;
};

float optQuantEven(
		float data[MAX_ENTRIES][DIMENSION],
		int numEntries, int numClusters, int index[MAX_ENTRIES],
		float out[MAX_ENTRIES][DIMENSION],
		float direction [DIMENSION],float *step
)
{
	int maxTry=MAX_TRY;
	int i,j,k;
	float t,s;
	float centered[MAX_ENTRIES][DIMENSION];
	float ordered[MAX_ENTRIES][DIMENSION];
	float mean[DIMENSION];
	float cov[DIMENSION][DIMENSION];
	float projected[MAX_ENTRIES];

	int order[MAX_ENTRIES];

	for (i=0;i<numEntries;i++)
		for (j=0;j<DIMENSION;j++)
			centered[i][j]=data[i][j];

	centerInPlace(centered, numEntries, mean);
	covariance(centered, numEntries, cov);

	// check if they all are the same

	t=0;
	for (j=0;j<DIMENSION;j++)
		t+= cov[j][j];

	if (t==0 || numEntries==0) {
		for (i=0;i<numEntries;i++) {
			index[i]=0;
			for (j=0;j<DIMENSION;j++)
				out[i][j]=mean[j];
		}
		return 0.;
	}

	eigenVector(cov, direction);
	project(centered, numEntries, direction, projected);

	for (i=0;i<maxTry;i++) {

		if (i) {
			t=0;
			for (j=0;j<DIMENSION;j++) {
				direction[j]=0;
				for (k=0;k<numEntries;k++)
					direction[j]+=ordered[k][j]*index[k];
				t+=direction[j]*direction[j];
			}

			// Actually we don't need to normailize direction here, as the
			// optimal quntization (index) is invariant of the scale.
			// Hence we don't care about possible degenration of the <direction> either
			// though normally it should not happen

			// However, the EPSILON should be scaled, otherwise is does not make sense

			t = sqrt(t)*EPSILON;

			project(centered, numEntries, direction, projected);

			for (j=1; j < numEntries;j++)
				if (projected[order[j]] < projected[order[j-1]]-t /*EPSILON*/)
					break;

			if (j >= numEntries)
				break;
		}

		sortProjection(projected, order, numEntries);

		for (k=0;k<numEntries;k++)
			for (j=0;j<DIMENSION;j++)
				ordered[k][j]=centered[order[k]][j];

		quantEven(ordered, numEntries, numClusters, index);

	}
	s=t=0;

	float q=0;

	for (k=0;k<numEntries;k++) {
		s+= index[k];
		t+= index[k]*index[k];
	}

	for (j=0;j<DIMENSION;j++) {
		direction[j]=0;
		for (k=0;k<numEntries;k++)
			direction[j]+=ordered[k][j]*index[k];
		q+= direction[j]* direction[j];

	}



	s /= (float) numEntries;

	t = t - s * s * (float) numEntries;

	ASSERT(t !=0);


	t = (t == 0 ? 0. : 1/t);

	for (i=0;i<numEntries;i++)
		for (j=0;j<DIMENSION;j++)
			out[order[i]][j]=mean[j]+direction[j]*t*(index[i]-s);

	// normalize direction for output

	q=sqrt(q);
	*step=t*q;
	for (j=0;j<DIMENSION;j++)
		direction[j]/=q;
	return totalError(data,out,numEntries);
};


int requantize(float data[MAX_ENTRIES][DIMENSION],
							 float centers[MAX_CLUSTERS][DIMENSION], int numEntries, int numClusters,int index[MAX_ENTRIES] )
{
	int i,j,k;
	float p,q;
	int cnt[MAX_CLUSTERS];
	int change =0;

	for (i=0;i<numEntries;i++) {
		p=0;
		index[i]=0;
		for(k=0;k<DIMENSION;k++)
			p+=(data[i][k]-centers[index[i]][k])*(data[i][k]-centers[index[i]][k]);

		for(j=0;j<numClusters;j++) {
			q=0;
			for(k=0;k<DIMENSION;k++)
				q+=(data[i][k]-centers[j][k])*(data[i][k]-centers[j][k]);

			change |= q < p ? (j!= index[i]) : 0;
			index[i]= q < p ? j : index[i];
			p    = q < p ? q : p;
		}
	}

	for(j=0;j<numClusters;j++)
		cnt[j]=0;

	for(j=0;j<numClusters;j++)
		for(k=0;k<DIMENSION;k++)
			centers[j][k]=0;

	for (i=0;i<numEntries;i++) {
		cnt[index[i]]++;
		for(k=0;k<DIMENSION;k++)
			centers[index[i]][k]+=data[i][k];
	}
	for(j=0;j<numClusters;j++)
		for(k=0;k<DIMENSION;k++)
			centers[j][k]/=(float) cnt[j];

	return(change);
}

float optQuantLineConstr(
		float data[MAX_ENTRIES][DIMENSION],
		int numEntries, int numClusters, int index[MAX_ENTRIES],
		float out[MAX_ENTRIES][DIMENSION]
)
{

	int maxTry=MAX_TRY;

	int i,j,k;
	float t;

	float centered[MAX_ENTRIES][DIMENSION];

	float mean[DIMENSION];

	float cov[DIMENSION][DIMENSION];

	float projected[MAX_ENTRIES];

	float direction [DIMENSION];

	int order[MAX_ENTRIES];

	for (i=0;i<numEntries;i++)
		for (j=0;j<DIMENSION;j++)
			centered[i][j]=data[i][j];

	centerInPlace(centered, numEntries, mean);
	covariance(centered, numEntries, cov);

	// check if they all are the same

	t=0;
	for (j=0;j<DIMENSION;j++)
		t+= cov[j][j];

	if (t==0 || numEntries==0) {
		for (i=0;i<numEntries;i++) {
			index[i]=0;
			for (j=0;j<DIMENSION;j++)
				out[i][j]=mean[j];
		}
		return 0.;
	}

	eigenVector(cov, direction);
	project(centered, numEntries, direction, projected);

	for (i=0;i<maxTry;i++) {

		if (i) {
			t=0;
			for (j=0;j<DIMENSION;j++) {
				direction[j]=0;
				for (k=0;k<numEntries;k++)
					direction[j]+=centered[k][j]*index[k];
				t=direction[j]*direction[j];
			}

			// Actually we don't need to normailize direction here, as the
			// optimal quntization (index) is invariant of the scale.
			// Hence we don't care about possible degenration of the <direction> either
			// though normally it should not happen

			// However, the EPSILON should be scaled, otherwise is does not make sense

			t = sqrt(t)*EPSILON;

			project(centered, numEntries, direction, projected);

			for (j=1; j < numEntries;j++)
				if (projected[order[j]] < projected[order[j-1]]-t /*EPSILON*/)
					break;

			if (j >= numEntries)
				break;
		}

		sortProjection(projected, order, numEntries);

		quantLineConstr(centered, order, numEntries, numClusters, index);

	}
	float gcAcc[MAX_CLUSTERS][DIMENSION];
	float gcSAcc[MAX_CLUSTERS][DIMENSION];
	float gcS[MAX_CLUSTERS];


	for(i=0;i<MAX_CLUSTERS;i++) {
		gcS[i]=0;
		for(j=0;j<DIMENSION;j++)
			gcAcc[i][j]=0;
	}

	for (k=0;k<numEntries;k++) {
		gcS[index[k]]+=1;
		for (j=0;j<DIMENSION;j++)
			gcAcc[index[k]][j]+=centered[k][j];

	}

	for(i=0;i<numClusters;i++)
		for (j=0;j<DIMENSION;j++)
			if (gcS[i]!=0) {
				gcSAcc[i][j] = gcAcc[i][j]/sqrt((float)gcS[i]);
				gcAcc[i][j] /= ((float)gcS[i]);
			}
			else
				gcSAcc[i][j] = 0;


	covariance(gcSAcc, numClusters,  cov);
	eigenVector(cov, direction);
	// assume the vector is normalized here

	for(i=0;i<numClusters;i++) {
		gcS[i]=0;
		for (j=0;j<DIMENSION;j++)
			gcS[i]+=direction[j]*gcAcc[i][j];
	}

	for (i=0;i<numEntries;i++)
		for (j=0;j<DIMENSION;j++)
			out[i][j]=mean[j]+direction[j]*gcS[index[i]];

	return totalError(data,out,numEntries);
};

void quantTrace(float data[MAX_ENTRIES_QUANT_TRACE][DIMENSION],int numEntries, int numClusters, int index[MAX_ENTRIES_QUANT_TRACE]) {
	// Data should be centered, otherwise will not work
	int i,j,k;
	float sdata[2*MAX_ENTRIES][DIMENSION];
	float  dpAcc [DIMENSION];
	float M =0;
	struct TRACE  *tr ;

	tr=amd_trs[numClusters-1][numEntries-1];

	int trcnt =trcnts[numClusters-1][numEntries-1];

	int *code;
	code=amd_codes[numClusters-1][numEntries-1];

	for (i=0;i<numEntries;i++)
		for (j=0;j<DIMENSION;j++) {
			sdata[2*i][j]= data[i][j];
			sdata[2*i+1][j]=-data[i][j];
		}

	for (j=0;j<DIMENSION;j++)
		dpAcc[j]=0;

	k=-1;

#define UROLL_STEP(i) \
    dpAcc[0]+=sdata[tr[i].k][0];\
    dpAcc[1]+=sdata[tr[i].k][1];\
    dpAcc[2]+=sdata[tr[i].k][2];\
    { float c; \
    c = (dpAcc[0]*dpAcc[0]+dpAcc[1]*dpAcc[1]+dpAcc[2]*dpAcc[2])*tr[i].d;\
    if (c > M) {k=i;M=c;};};

	for (i=0;i+15<trcnt;i+=16) {
		UROLL_STEP(i)
		UROLL_STEP(i+1)
		UROLL_STEP(i+2)
		UROLL_STEP(i+3)
		UROLL_STEP(i+4)
		UROLL_STEP(i+5)
		UROLL_STEP(i+6)
		UROLL_STEP(i+7)
		UROLL_STEP(i+8)
		UROLL_STEP(i+9)
		UROLL_STEP(i+10)
		UROLL_STEP(i+11)
		UROLL_STEP(i+12)
		UROLL_STEP(i+13)
		UROLL_STEP(i+14)
		UROLL_STEP(i+15)
	}

	for (;i<trcnt;i++) {
		UROLL_STEP(i)
	}

	if ((k<0)||(k >=MAX_TRACE)) {  // NP
		return;
	}

	k = code[k];
	i=0;
	for (j=0;j<numEntries;j++) {
		while ((k & 1) ==0) {
			i++;
			k>>=1;
		}
		index[j]=i;
		k>>=1;
	}
}

void quantTrace_d(float data[MAX_ENTRIES_QUANT_TRACE][MAX_DIMENSION_BIG],int numEntries, int numClusters, int index[MAX_ENTRIES_QUANT_TRACE],int dimension)
{
	// Data should be centered, otherwise will not work

	int i,j,k;

	float sdata[2*MAX_ENTRIES][MAX_DIMENSION_BIG];

	float  dpAcc [MAX_DIMENSION_BIG];

	float M =0;

	struct TRACE  *tr ;
	tr=amd_trs[numClusters-1][numEntries-1];

	int trcnt =trcnts[numClusters-1][numEntries-1];

	int *code;
	code=amd_codes[numClusters-1][numEntries-1];

	for (i=0;i<numEntries;i++)
		for (j=0;j<dimension;j++)
		{
			sdata[2*i][j]= data[i][j];
			sdata[2*i+1][j]=-data[i][j];
		}

	for (j=0;j<dimension;j++)
		dpAcc[j]=0;

	k=-1;

#define UROLL_STEP_1(i) \
    dpAcc[0]+=sdata[tr[i].k][0];\
    {\
        float c; \
        c = (dpAcc[0]*dpAcc[0])*tr[i].d;\
        if (c > M) {k=i;M=c;};\
    };

#define UROLL_STEP_2(i) \
    dpAcc[0]+=sdata[tr[i].k][0];\
    dpAcc[1]+=sdata[tr[i].k][1];\
    { float c; \
    c = (dpAcc[0]*dpAcc[0]+dpAcc[1]*dpAcc[1])*tr[i].d;\
    if (c > M) {k=i;M=c;};};

#define UROLL_STEP_3(i) \
    dpAcc[0]+=sdata[tr[i].k][0];\
    dpAcc[1]+=sdata[tr[i].k][1];\
    dpAcc[2]+=sdata[tr[i].k][2];\
    { float c; \
    c = (dpAcc[0]*dpAcc[0]+dpAcc[1]*dpAcc[1]+dpAcc[2]*dpAcc[2])*tr[i].d;\
    if (c > M) {k=i;M=c;};};

#define UROLL_STEP_4(i) \
    dpAcc[0]+=sdata[tr[i].k][0];\
    dpAcc[1]+=sdata[tr[i].k][1];\
    dpAcc[2]+=sdata[tr[i].k][2];\
    dpAcc[3]+=sdata[tr[i].k][3];\
    { float c; \
    c = (dpAcc[0]*dpAcc[0]+dpAcc[1]*dpAcc[1]+dpAcc[2]*dpAcc[2]+dpAcc[3]*dpAcc[3])*tr[i].d;\
    if (c > M) {k=i;M=c;};};

#undef UROLL_STEP

#define UROLL_MACRO(UROLL_STEP){\
\
\
    for (i=0;i+15<trcnt;i+=16)\
    {\
        UROLL_STEP(i)\
        UROLL_STEP(i+1)\
        UROLL_STEP(i+2)\
        UROLL_STEP(i+3)\
        UROLL_STEP(i+4)\
        UROLL_STEP(i+5)\
        UROLL_STEP(i+6)\
        UROLL_STEP(i+7)\
        UROLL_STEP(i+8)\
        UROLL_STEP(i+9)\
        UROLL_STEP(i+10)\
        UROLL_STEP(i+11)\
        UROLL_STEP(i+12)\
        UROLL_STEP(i+13)\
        UROLL_STEP(i+14)\
        UROLL_STEP(i+15)\
    }\
\
    for (;i<trcnt;i++) {\
        UROLL_STEP(i)\
    }};

	switch(dimension)
	{
	case    1:
	UROLL_MACRO(UROLL_STEP_1);
		break;
	case    2:
	UROLL_MACRO(UROLL_STEP_2);
		break;
	case    3:
	UROLL_MACRO(UROLL_STEP_3);
		break;
	case    4:
	UROLL_MACRO(UROLL_STEP_4);
		break;
	default:
		return;
		break;
	}


	if (k<0)
	{
		return;
	}

	k = code[k];
	i=0;
	for (j=0;j<numEntries;j++)
	{
		while ((k & 1) ==0)
		{
			i++;
			k>>=1;
		}
		index[j]=i;

		k>>=1;
	}
}

void quant_AnD_Shell(float* v_, int k, int n, int *idx) {
	// input:
	//
	// v_  points, might be uncentered
	// k - number of points in the ramp
	// n - number of points in v_
	//
	// output:
	//
	// index, uncentered, in the range 0..k-1
	//
#define MAX_BLOCK MAX_ENTRIES
	int i,j;
	float v[MAX_BLOCK];
	float z[MAX_BLOCK];
	a d[MAX_BLOCK];
	float l;
	float mm;
	float r=0;
	int mi;

	ASSERT((v_ != NULL) && (n>1) && (k>1));

	float m, M, s, dm=0.;
	m=M=v_[0];

	for (i=1; i < n;i++) {
		m = m < v_[i] ? m : v_[i];
		M = M > v_[i] ? M : v_[i];
	}
	if (M==m) {
		for (i=0; i < n;i++)
			idx[i]=0;
		return;
	}

	ASSERT(M-m >0);
	s = (k-1)/(M-m);
	for (i=0; i < n;i++) {
		v[i] = v_[i]*s;

		idx[i]=(int)(z[i] = floor(v[i] +0.5 /* stabilizer*/ - m *s));

		d[i].d = v[i]-z[i]- m *s;
		d[i].i = i;
		dm+= d[i].d;
		r += d[i].d*d[i].d;
	}
	if (n*r- dm*dm >= (float)(n-1)/4 /*slack*/ /2) {

		dm /= (float)n;

		for (i=0; i < n;i++)
			d[i].d -= dm;

		qsort((void*)&d, n, sizeof(a),a_compare);

		// got into fundamental simplex
		// move coordinate system origin to its center
		for (i=0; i < n;i++)
			d[i].d -= (2.*(float)i+1-(float)n)/2./(float)n;

		mm=l=0.;
		j=-1;
		for (i=0; i < n;i++) {
			l+=d[i].d;
			if (l < mm) {
				mm =l;
				j=i;
			}
		}

		// position which should be in 0
		j = ++j % n;

		for (i=j; i < n;i++)
			idx[d[i].i]++;
	}
	// get rid of an offset in idx
	mi=idx[0];
	for (i=1; i < n;i++)
		mi = mi < idx[i]? mi :idx[i];

	for (i=0; i < n;i++)
		idx[i]-=mi;
}

float optQuantTrace(
		float data[MAX_ENTRIES][DIMENSION],
		int numEntries, int numClusters, int index_[MAX_ENTRIES],
		float out[MAX_ENTRIES][DIMENSION],
		float direction [DIMENSION],float *step
)
{
	
	int index[MAX_ENTRIES];

	int maxTry=MAX_TRY;

	int i,j,k;
	float t,s;

	float centered[MAX_ENTRIES][DIMENSION];

	float ordered[MAX_ENTRIES][DIMENSION];

	float mean[DIMENSION];

	float cov[DIMENSION][DIMENSION];

	float projected[MAX_ENTRIES];

	int order[MAX_ENTRIES];


	for (i=0;i<numEntries;i++)
		for (j=0;j<DIMENSION;j++)
			centered[i][j]=data[i][j];

	centerInPlace(centered, numEntries, mean);
	covariance(centered, numEntries, cov);

	// check if they all are the same

	t=0;
	for (j=0;j<DIMENSION;j++)
		t+= cov[j][j];

	if (t==0 || numEntries==0) {
		for (i=0;i<numEntries;i++) {
			index_[i]=0;
			for (j=0;j<DIMENSION;j++)
				out[i][j]=mean[j];
		}
		return 0.;
	}


	eigenVector(cov, direction);
	project(centered, numEntries, direction, projected);


	for (i=0;i<maxTry;i++) {

		if (i) {
			t=0;
			for (j=0;j<DIMENSION;j++) {
				direction[j]=0;
				for (k=0;k<numEntries;k++)
					direction[j]+=ordered[k][j]*index[k];
				t+=direction[j]*direction[j];
			}

			// Actually we don't need to normailize direction here, as the
			// optimal quntization (index) is invariant of the scale.
			// Hence we don't care about possible degenration of the <direction> either
			// though normally it should not happen

			// However, the EPSILON should be scaled, otherwise is does not make sense

			t = sqrt(t)*EPSILON;

			project(centered, numEntries, direction, projected);

			for (j=1; j < numEntries;j++)
				if (projected[order[j]] < projected[order[j-1]]-t /*EPSILON*/)
					break;

			if (j >= numEntries)
				break;
		}

		sortProjection(projected, order, numEntries);

		for (k=0;k<numEntries;k++)
			for (j=0;j<DIMENSION;j++)
				ordered[k][j]=centered[order[k]][j];

		quantTrace(ordered, numEntries, numClusters, index);
	}
	s=t=0;

	float q=0;

	for (k=0;k<numEntries;k++) {
		s+= index[k];
		t+= index[k]*index[k];
	}

	for (j=0;j<DIMENSION;j++) {
		direction[j]=0;
		for (k=0;k<numEntries;k++)
			direction[j]+=ordered[k][j]*index[k];
		q+= direction[j]* direction[j];

	}


	s /= (float) numEntries;

	t = t - s * s * (float) numEntries;

	ASSERT(t !=0);


	t = (t == 0 ? 0. : 1/t);

	for (i=0;i<numEntries;i++) {
		for (j=0;j<DIMENSION;j++)
			out[order[i]][j]=mean[j]+direction[j]*t*(index[i]-s);
		index_[order[i]]=index[i];
	}


	// normalize direction for output

	q=sqrt(q);
	*step=t*q;
	for (j=0;j<DIMENSION;j++)
		direction[j]/=q;

	return totalError(data,out,numEntries);
}

float optQuantTrace_d(
		float data[MAX_ENTRIES][MAX_DIMENSION_BIG],
		int numEntries, int numClusters, int index_[MAX_ENTRIES],
		float out[MAX_ENTRIES][MAX_DIMENSION_BIG],
		float direction [MAX_DIMENSION_BIG],float *step,
		int dimension
)
{
	int index[MAX_ENTRIES];
	int maxTry=MAX_TRY;
	int i,j,k;
	float t,s;
	float centered[MAX_ENTRIES][MAX_DIMENSION_BIG];
	float ordered[MAX_ENTRIES][MAX_DIMENSION_BIG];
	float mean[MAX_DIMENSION_BIG];
	float cov[DIMENSION][MAX_DIMENSION_BIG];
	float projected[MAX_ENTRIES];
	int order[MAX_ENTRIES];

	for (i=0;i<numEntries;i++)
		for (j=0;j<dimension;j++)
			centered[i][j]=data[i][j];

	centerInPlace_d(centered, numEntries, mean, dimension);
	covariance_d(centered, numEntries, cov, dimension);

	// check if they all are the same

	t=0;
	for (j=0;j<dimension;j++)
		t+= cov[j][j];

	if (t<EPSILON || numEntries==0) {
		for (i=0;i<numEntries;i++) {
			index_[i]=0;
			for (j=0;j<dimension;j++)
				out[i][j]=mean[j];
		}
		return 0.;
	}


	eigenVector_d(cov, direction, dimension);
	project_d(centered, numEntries, direction, projected, dimension);

	for (i=0;i<maxTry;i++)
	{
		if (i)
		{
			t=0;
			for (j=0;j<dimension;j++)
			{
				direction[j]=0;
				for (k=0;k<numEntries;k++)
					direction[j]+=ordered[k][j]*index[k];
				t+=direction[j]*direction[j];
			}

			// Actually we don't need to normailize direction here, as the
			// optimal quntization (index) is invariant of the scale.
			// Hence we don't care about possible degenration of the <direction> either
			// though normally it should not happen

			// However, the EPSILON should be scaled, otherwise is does not make sense

			t = sqrt(t)*EPSILON;

			project_d(centered, numEntries, direction, projected, dimension);

			for (j=1; j < numEntries;j++)
				if (projected[order[j]] < projected[order[j-1]]-t /*EPSILON*/)
					break;

			if (j >= numEntries)
				break;
		}

		sortProjection(projected, order, numEntries);

		for (k=0;k<numEntries;k++)
			for (j=0;j<dimension;j++)
				ordered[k][j]=centered[order[k]][j];

		quantTrace_d(ordered, numEntries, numClusters, index, dimension);
	}

	s=t=0;

	float q=0;

	for (k=0;k<numEntries;k++)
	{
		s+= index[k];
		t+= index[k]*index[k];
	}

	for (j=0;j<dimension;j++)
	{
		direction[j]=0;
		for (k=0;k<numEntries;k++)
			direction[j]+=ordered[k][j]*index[k];
		q+= direction[j]* direction[j];

	}

	s /= (float) numEntries;

	t = t - s * s * (float) numEntries;

	ASSERT(t !=0);

	t = (t == 0 ? 0. : 1/t);

	for (i=0;i<numEntries;i++)
	{
		for (j=0;j<dimension;j++)
			out[order[i]][j]=mean[j]+direction[j]*t*(index[i]-s);
		index_[order[i]]=index[i];
	}

	// normalize direction for output

	q=sqrt(q);
	*step=t*q;

	for (j=0;j<dimension;j++)
		direction[j]/=q;

	return totalError_d(data,out,numEntries, dimension);
}


void traceBuilder (int numEntries, int numClusters,struct TRACE tr [], int code[], int *trcnt )
{
	//=================
#define DIG(J_IN,I,N,J_OUT,DIR,NC,NCC)            \
        for (I=J_IN;I<N || NC < NCC ;I++) {                   \
            J_OUT = ((((J_IN) & 0x1)==DIR) ? I : N-1-(I-(J_IN)));    \

	//=================

	int i[7];
	int j[7];
	int k[7]={0,1,2,3,4,5,6};
	int n;
	int c  =0;
	int c0 =0;
	int p;
	int h[8]={0,0,0,0, 0,0,0,0};

	if (numClusters == 1)
	{
		tr[c].k=0;
		tr[c].d=0;
		code[c]=0;
		*trcnt=0;
		return;
	}

	h[numClusters-1]=numEntries;

	int q  = numEntries*(numClusters-1);
	int q2 = numEntries*(numClusters-1)*(numClusters-1);

	n = numEntries + numClusters -2; // higest delimiter postion; all points start in highest cluster

	int cd = -(1<< (numClusters-1));

	DIG(     0,i[0],n,j[0],0,numClusters,2)
		DIG(j[0]+1,i[1],n,j[1],1,numClusters,3)
			DIG(j[1]+1,i[2],n,j[2],0,numClusters,4)
				DIG(j[2]+1,i[3],n,j[3],1,numClusters,5)
					DIG(j[3]+1,i[4],n,j[4],0,numClusters,6)
						DIG(j[4]+1,i[5],n,j[5],1,numClusters,7)
							DIG(j[5]+1,i[6],n,j[6],0,numClusters,8)

								int rescan;
								do {
									rescan=0;
									for (p=0;p<numClusters-1;p++) {

										if (abs(j[p]-k[p]) >1 )  {
											return;
										}

										else if (j[p]-k[p]== 1 ) {
											int ci= k[p]-p;                    // move it one cluster down "-"
											int cn=p+1;

											h[cn]--;
											h[cn-1]++;

											if (h[cn] < 0 || h[cn-1]>= numEntries) {
												rescan =1;
												h[cn]++;
												h[cn-1]--;

											}

											else {

												q2+= -2*cn+1;
												q--;

												{
													int i1,cc=0; for(i1=0;i1<numClusters;i1++) cc += i1*i1*h[i1];
												};
												
												cd |=  (1<<k[p]);
												cd &= ~(1<<j[p]);

												if (c < MAX_TRACE) // NP
												{
													tr[c].k=2*ci+1;
													tr[c].d=1./((float) q2 - (float) q*(float) q /(float) (numEntries));
													code[c]=cd;
													c++;
												}
												else
												{
													// What to do here?
													tr[c].k=0;
													tr[c].d=0;
													code[c]=0;
													*trcnt=0;
													return;
												}
												k[p]=j[p];
											}
										}
										else if (j[p]-k[p]==-1 )
										{
											int ci=j[p]-p;                    // move it up
											int cn =p;

											h[cn]--;
											h[cn+1]++;

											if (h[cn] < 0 || h[cn+1]>= numEntries) {
												rescan =1;
												h[cn]++;
												h[cn+1]--;
											}
											else {


												q2+= 2*cn+1;
												q++;
												{ 
													int i1,cc=0; for(i1=0;i1<numClusters;i1++) cc += i1*i1*h[i1];
												};

												cd |=  (1<<k[p]);
												cd &= ~(1<<j[p]);

												if (c < MAX_TRACE) // NP
												{
													tr[c].k=2*ci;
													tr[c].d=1./((float) q2 - (float) q*(float) q /(float) (numEntries));
													code[c]=cd;
													c++;
												}
												else
												{
													// What to do here?
													tr[c].k=0;
													tr[c].d=0;
													code[c]=0;
													*trcnt=0;
													return;
												}
												k[p]=j[p];
											}
										}
									}
								}
								while (rescan);
								c0++;
								if (numClusters < 8) break; }
							if (numClusters < 7) break; }
						if (numClusters < 6) break; }
					if (numClusters < 5) break; }
				if (numClusters < 4) break; }
			if (numClusters < 3) break; }
		if (numClusters < 2) break; }

	*trcnt=c;
}

float optQuantAnD(
		float data[MAX_ENTRIES][DIMENSION],
		int numEntries, int numClusters, int index[MAX_ENTRIES],
		float out[MAX_ENTRIES][DIMENSION],
		float direction [DIMENSION],float *step
)
{
	int index_[MAX_ENTRIES];
	int maxTry=MAX_TRY*10;
	int try_two=50;
	int i,j,k;
	float t,s;
	float centered[MAX_ENTRIES][DIMENSION];
	float mean[DIMENSION];
	float cov[DIMENSION][DIMENSION];
	float projected[MAX_ENTRIES];

	int order_[MAX_ENTRIES];


	for (i=0;i<numEntries;i++)
		for (j=0;j<DIMENSION;j++)
			centered[i][j]=data[i][j];

	centerInPlace(centered, numEntries, mean);
	covariance(centered, numEntries, cov);

	// check if they all are the same

	t=0;
	for (j=0;j<DIMENSION;j++)
		t+= cov[j][j];

	if (t==0 || numEntries==0) {
		for (i=0;i<numEntries;i++) {
			index[i]=0;
			for (j=0;j<DIMENSION;j++)
				out[i][j]=mean[j];
		}
		return 0.;
	}


	eigenVector(cov, direction);
	project(centered, numEntries, direction, projected);



	for (i=0;i<maxTry;i++) {
		int done =0;

		if (i) {
			do {
				float q;
				q=s=t=0;

				for (k=0;k<numEntries;k++) {
					s+= index[k];
					t+= index[k]*index[k];
				}

				for (j=0;j<DIMENSION;j++) {
					direction[j]=0;
					for (k=0;k<numEntries;k++)
						direction[j]+=centered[k][j]*index[k];
					q+= direction[j]* direction[j];

				}

				s /= (float) numEntries;
				t = t - s * s * (float) numEntries;
				ASSERT(t !=0);
				t = (t == 0 ? 0. : 1/t);
				// We need to requantize

				q = sqrt(q);
				t *=q;

				if (q !=0)
					for (j=0;j<DIMENSION;j++)
						direction[j]/=q;

				// direction normalized

				project(centered, numEntries, direction, projected);
				sortProjection(projected, order_, numEntries);

				int index__[MAX_ENTRIES];

				// it's projected and centered; cluster centers are (index[i]-s)*t (*dir)
				k=0;
				for (j=0; j < numEntries;j++) {
					while (projected[order_[j]] > (k+0.5 -s)*t  && k < numClusters-1)
						k++;
					index__[order_[j]]=k;
				}
				done =1;
				for (j=0; j < numEntries;j++) {
					done = (done && (index__[j]==index[j]));
					index[j]=index__[j];
				}
			} while (! done && try_two--);
			if (i==1)
				for (j=0; j < numEntries;j++)
					index_[j]=index[j];
			else {
				done =1;
				for (j=0; j < numEntries;j++) {
					done = (done && (index_[j]==index[j]));
					index_[j]=index_[j];
				}
				if (done)
					break;

			}
		}

		quant_AnD_Shell(projected,  numClusters,numEntries, index);

	}
	s=t=0;

	float q=0;

	for (k=0;k<numEntries;k++) {
		s+= index[k];
		t+= index[k]*index[k];
	}

	for (j=0;j<DIMENSION;j++) {
		direction[j]=0;
		for (k=0;k<numEntries;k++)
			direction[j]+=centered[k][j]*index[k];
		q+= direction[j]* direction[j];

	}

	s /= (float) numEntries;

	t = t - s * s * (float) numEntries;
	ASSERT(t !=0);


	t = (t == 0 ? 0. : 1/t);

	for (i=0;i<numEntries;i++)
		for (j=0;j<DIMENSION;j++)
			out[i][j]=mean[j]+direction[j]*t*(index[i]-s);

	// normalize direction for output

	q=sqrt(q);
	*step=t*q;
	for (j=0;j<DIMENSION;j++)
		direction[j]/=q;


	return totalError(data,out,numEntries);
}

float optQuantAnD_d(
		float data[MAX_ENTRIES][MAX_DIMENSION_BIG],
		int numEntries, int numClusters, int index[MAX_ENTRIES],
		float out[MAX_ENTRIES][MAX_DIMENSION_BIG],
		float direction [MAX_DIMENSION_BIG],float *step,
		int dimension
)
{
	int index_[MAX_ENTRIES];

	int maxTry=MAX_TRY*10;
	int try_two=50;

	int i,j,k;
	float t,s;

	float centered[MAX_ENTRIES][MAX_DIMENSION_BIG];

	float mean[MAX_DIMENSION_BIG];

	float cov[MAX_DIMENSION_BIG][MAX_DIMENSION_BIG];

	float projected[MAX_ENTRIES];

	int order_[MAX_ENTRIES];


	for (i=0;i<numEntries;i++)
		for (j=0;j<dimension;j++)
			centered[i][j]=data[i][j];

	centerInPlace_d(centered, numEntries, mean, dimension);
	covariance_d(centered, numEntries, cov, dimension);

	// check if they all are the same

	t=0;
	for (j=0;j<dimension;j++)
		t+= cov[j][j];

	if (t<(1./256.) || numEntries==0) {
		for (i=0;i<numEntries;i++) {
			index[i]=0;
			for (j=0;j<dimension;j++)
				out[i][j]=mean[j];
		}
		return 0.;
	}

	eigenVector_d(cov, direction, dimension);
	project_d(centered, numEntries, direction, projected, dimension);

	for (i=0;i<maxTry;i++)
	{
		int done =0;

		if (i)
		{
			do
			{
				float q;
				q=s=t=0;

				for (k=0;k<numEntries;k++)
				{
					s+= index[k];
					t+= index[k]*index[k];
				}

				for (j=0;j<dimension;j++)
				{
					direction[j]=0;
					for (k=0;k<numEntries;k++)
						direction[j]+=centered[k][j]*index[k];
					q+= direction[j]* direction[j];

				}

				s /= (float) numEntries;
				t = t - s * s * (float) numEntries;
				ASSERT(t !=0);
				t = (t == 0 ? 0. : 1/t);
				// We need to requantize

				q = sqrt(q);
				t *=q;

				if (q !=0)
					for (j=0;j<dimension;j++)
						direction[j]/=q;

				// direction normalized

				project_d(centered, numEntries, direction, projected, dimension);
				sortProjection(projected, order_, numEntries);

				int index__[MAX_ENTRIES];

				// it's projected and centered; cluster centers are (index[i]-s)*t (*dir)
				k=0;
				for (j=0; j < numEntries;j++)
				{
					while (projected[order_[j]] > (k+0.5 -s)*t  && k < numClusters-1)
						k++;
					index__[order_[j]]=k;
				}
				done =1;
				for (j=0; j < numEntries;j++)
				{
					done = (done && (index__[j]==index[j]));
					index[j]=index__[j];
				}
			} while (! done && try_two--);

			if (i==1)
				for (j=0; j < numEntries;j++)
					index_[j]=index[j];
			else
			{
				done =1;
				for (j=0; j < numEntries;j++)
				{
					done = (done && (index_[j]==index[j]));
					index_[j]=index_[j];
				}
				if (done)
					break;

			}
		}

		quant_AnD_Shell(projected,  numClusters,numEntries, index);
	}
	s=t=0;

	float q=0;

	for (k=0;k<numEntries;k++)
	{
		s+= index[k];
		t+= index[k]*index[k];
	}

	for (j=0;j<dimension;j++)
	{
		direction[j]=0;
		for (k=0;k<numEntries;k++)
			direction[j]+=centered[k][j]*index[k];
		q+= direction[j]* direction[j];
	}

	s /= (float) numEntries;

	t = t - s * s * (float) numEntries;

	ASSERT(t !=0);

	t = (t == 0 ? 0. : 1/t);

	for (i=0;i<numEntries;i++)
		for (j=0;j<dimension;j++)
			out[i][j]=mean[j]+direction[j]*t*(index[i]-s);

	// normalize direction for output

	q=sqrt(q);
	*step=t*q;
	for (j=0;j<dimension;j++)
		direction[j]/=q;

	return totalError_d(data,out,numEntries, dimension);
}

