#include "generateEllipseCandidates.h"
#include "candidatesFunctions.cpp"
#include "gradientFunctions.cpp"
#include "lsdFunctions.cpp"
#include "lsdInterface.cpp"
#include "mexFunction.cpp"
#include "miscFunctions.cpp"
#include "nfaFunctions.cpp"
#include "regionsFunctions.cpp"


//==============================================================================
//====================================================================================================
//================================clustering==========================================================
//聚类
//求points中第i行与initializations中第j行里每个元素的平方差总和,每行维度都为nDims
inline double squaredDifference(int & nDims, double *& points, int & i, double *& initializations, int & j)
{
    double result = 0;
    for (int k = 0; k < nDims; ++k)
		result += pow(points[i*nDims+k] - initializations[j*nDims+k], 2);
    return result;
}
/**
 *输入
 *points: 待均值漂移的点集，总共有nPoints个点，每个点有nDims维度，是一维数组
 *initPoints: 均值漂移初始化位置，在nxd空间中找均值漂移初始时开始搜索的位置，总共有initLength个点，每个点有nDims维度
 *sigma = 1
 *window_size: window parameter = distance_tolerance或者window parameter = distance_tolerance/2
 *accuracy_tolerance: 收敛容忍误差1e-6
 *iter_times: 迭代次数50
 *输出
 *收敛的位置，位置个数与初始化搜索位置个数一样,我们将结果更新到initPoints,也就是它既是输入参数，也是输出参数，节省内存
 */
void meanShift( double * points, int nPoints, int nDims, double * & initPoints, int initLength, double sigma, double window_size, double accuracy_tolerance, int iter_times )
{
//	for (int i = 0; i<initLength; i++)
//		cout<<initPoints[2*i]<<'\t'<<initPoints[2*i+1]<<endl;
    int nQuerries = initLength;
    double * initializations = (double*)malloc(nQuerries * nDims * sizeof(double));
    memcpy(initializations, initPoints , nQuerries * nDims * sizeof(double));//copy

    double sigma2 = sigma*sigma;//sigma平方
    double radius2 = window_size *window_size;//平方
    double tolerance = accuracy_tolerance;
    int iters, maxiters = iter_times;//最大迭代次数
   //返回与初始搜索点集一样大小的最终定位点集
    double * finals = (double*)malloc(nQuerries * nDims * sizeof(double));;//最终定位点集的指针
    memcpy(finals, initializations, nQuerries * nDims * sizeof(double));
	double * distances = (double*)malloc(nPoints*sizeof(double));
    //printf("meanShift: nPoints:%d \tnDims: %d \tnQuerries:%d \n",nPoints,nDims,nQuerries);//打印
    for (int loop = 0; loop < nQuerries; ++loop)
    {
        iters = 0;
        while (iters < maxiters)
        {
            bool flag = false;
            double denominator = 0;//分母
            for (int i = 0; i < nPoints; ++i)//对所有的点集进行遍历，找到落在搜索圆域内的点
            {
                distances[i] = squaredDifference(nDims, points, i, initializations, loop);//求距离的平方
                if (distances[i] <= radius2)//在第loop个搜索中心的以sqrt(radius2)为半径的圆域内
                {
                    flag = true;
                    denominator += exp(-distances[i] / sigma2);
                }
            }
            if (!flag)
                break;
            for (int j = 0; j < nDims; ++j)
				finals[loop*nDims+j] = 0;//对最终定位点集中的第loop个点的向量赋值为0
            for (int i = 0; i < nPoints; ++i)
                if (distances[i] <= radius2)
                {
                    for (int j = 0; j < nDims; ++j)//每个内点向量的以一定权值累加
						finals[loop*nDims+j] += exp(-distances[i] / sigma2) * points[i*nDims+j];
                }
            for (int j = 0; j < nDims; ++j)//权值归一化
				finals[loop*nDims+j] /= denominator;
            if (sqrt(squaredDifference(nDims, finals, loop, initializations, loop)) < tolerance)//相继两次的迭代中心在误差内了，则认为已经收敛，没必要再继续迭代
                break;
            iters = iters + 1;
            for (int j = 0; j < nDims; ++j)//更新迭代的搜索中心
				initializations[loop*nDims+j] = finals[loop*nDims+j];
        }
    }
	memcpy(initPoints, finals, nQuerries * nDims * sizeof(double));
    free(distances);
    free(initializations);
	free(finals);
}

/***
 *输入
 *points,待聚类的点集,为一维数组,nPoints个点，每个点维度是nDims
 *distance_threshold 决定聚类的距离阈值
 *输出 outPoints
 *聚类后的点集 nOutPoints x nDims 
 *该函数要千万注意，当被调用后，函数内部会多申请nOutPoints个double型的数组内存，在外边使用完毕后，切记free(outPoints).
 */
void clusterByDistance(double * points, int nPoints, int nDims, double distance_threshold,int number_control, double * & outPoints, int * nOutPoints)
{ 
	double threshold2 = distance_threshold*distance_threshold;
    std::vector<double*> centers;
    std::vector<int> counts;
    centers.clear();
    counts.clear();
	char * labeled = (char*)malloc(sizeof(char)*nPoints);
    memset(labeled, 0, nPoints * sizeof(char));//初始化bool型标签为0
	if(nPoints == 1)
	{
		centers.push_back((double*)malloc(sizeof(double)*nDims));
		for (int k = 0; k < nDims; ++k)
			centers[centers.size() - 1][k] = points[k];
        counts.push_back(1);
	}
	else
	{
		for (int i = 0; i < nPoints; ++i)
		{
		    if (!labeled[i])
			{
		        labeled[i] = 1;
				centers.push_back((double*)malloc(sizeof(double)*nDims));
		        counts.push_back(1);
		        for (int k = 0; k < nDims; ++k)
				{
				   centers[centers.size() - 1][k] = points[i*nDims+k];  
				}
		        for (int j = i+1; j < nPoints; ++j)
		        {
		            if (!labeled[j])
		            {
		                double d = 0;
		                for (int k = 0; k < nDims; ++k)
				            d += pow(centers[centers.size() - 1][k] / counts[centers.size() - 1] - points[j*nDims+k], 2);
		                if (d <= threshold2)
		                {
		                    ++counts[centers.size() - 1];
		                    for (int k = 0; k < nDims; ++k)
								centers[centers.size() - 1][k] += points[j*nDims+k];
		                    labeled[j] = 1;
							if(counts[centers.size() - 1] >= number_control)//聚类数量控制，防止均值中心漂的太远  圆心聚类时20  半径聚类时10
								break;
		                }
		            }
		        }
		    }
		}
	}
    free(labeled);
    centers.shrink_to_fit();
    counts.shrink_to_fit();
    int m = (int) centers.size();
    outPoints = (double*)malloc(sizeof(double)*m*nDims);
	(*nOutPoints) = m;
    for (unsigned int i = 0; i < centers.size(); ++i)
    {
        for (int j = 0; j < nDims; ++j)
		{
			outPoints[i*nDims+j] = centers[i][j] / counts[i];
//			cout<<out[i*nDims+j]<<'\t';
		}
//		cout<<endl;
        free(centers[i]);
    }
    centers.resize(0);
    counts.resize(0);
//	vector<double*>().swap(centers);//释放回收vector内存
//	vector<int>().swap(counts);
}

//聚类算法，均值漂移
//三个步骤，一是选取初始迭代点，二是均值漂移，三是去除重复点，从而得到聚类中心
//获得候选圆心的聚类中心(xi,yi)
//输入：
//points，一维点数据,长度为points_num x 2
//distance_tolerance,数据点聚类的半径
//输出：
//二维数据点的聚类中心 centers是一维double数组， 大小为 centers_num x 2
//正确返回值为1，出现错误为0. 例如points为空
//切记切记！！！ centers为函数内部申请的内存，用来返回centers_num个点的聚类中心，使用完后一定要释放，记住free(centers)！！！
int  cluster2DPoints( double * points, int points_num, double distance_tolerance, double * & centers, int * centers_num)
{
	double xmax,xmin,ymax,ymin,xdelta,ydelta;
	int nbins_x,nbins_y;
	int x,y;
	int i;
	unsigned int addr,addr2;
	xmax = ymax = 0;
	xmin = ymin = DBL_MAX;
	for( i = 0; i< points_num; i++ )
	{
		addr = 2*i;
		if( points[addr] > xmax)
			xmax = points[addr];
		if( points[addr] < xmin)
			xmin = points[addr];
		if( points[addr+1] > ymax)
			ymax = points[addr+1];
		if( points[addr+1] < ymin)
			ymin = points[addr+1];
	}
	xmax += xmax*0.02;//避免xdelta、ydelta为0
	xmin -= xmin*0.02;
	ymax += ymax*0.02;
	ymin -= ymin*0.02;
	xdelta = (xmax-xmin);
	ydelta = (ymax-ymin);//有问题，假设所有数据一样大，此处为0
	nbins_x = (int)ceil(xdelta/distance_tolerance);
	nbins_y = (int)ceil(ydelta/distance_tolerance);
	if(nbins_x <= 0 )
	{
		nbins_x = 1;//至少保留1个bin
		//error("generateCircleCandidates,nbins_x,nbins_y error");
	}
	if(nbins_y <= 0)
	{
		nbins_y = 1;
	}
	point2d1i * center_bins;
	center_bins = (point2d1i *)calloc(nbins_y*nbins_x, sizeof(point2d1i));//(x,y,z),x用来记sum(xi),y用来记sum(yi),z用来记落在格子里的数量
	memset(center_bins,0,sizeof(point2d1i)*nbins_y*nbins_x);
	if(center_bins == NULL)
		error("cluster2DPoints, not enough memory");
//	cout<<"2D原始数据:"<<points_num<<endl;
	for ( i = 0; i< points_num; i++ )//将圆心记录到格子里面，同时落在相应格子里面的数量++
	{
		addr = 2*i;

//		cout<<points[addr]<<'\t'<<points[addr+1]<<endl;

		x = (int)((points[addr]   - xmin)/xdelta*nbins_x+0.5);//四舍五入
		y = (int)((points[addr+1] - ymin)/ydelta*nbins_y+0.5);
		if( x >= nbins_x)
			x = nbins_x-1;
		if( y >= nbins_y)
			y = nbins_y-1;
		addr2 = y*nbins_x+x;
		center_bins[addr2].x += points[addr];
		center_bins[addr2].y += points[addr+1];
		center_bins[addr2].z ++;
	}
	int initCentersLength = 0;
	for ( y = 0; y<nbins_y; y++)//将vote后非0的格子里面的point取均值，并按照顺序重写到center_bins里面，无内存消耗
		for ( x = 0; x<nbins_x; x++)
		{
			addr = y*nbins_x+x;
			if(center_bins[addr].z > 0)
			{
				center_bins[initCentersLength].x = center_bins[addr].x/center_bins[addr].z;
				center_bins[initCentersLength].y = center_bins[addr].y/center_bins[addr].z;
				initCentersLength++;
			}
		}
	if(initCentersLength == 0)
	{
		(*centers_num) = 0;
		centers = NULL;
		//cout<<"cluster2DPoints,points number:"<<points_num<<endl;
		//cout<<"cluster2DPoints,initCentersLength equals 0"<<endl;
		return 0;
		//error("generateCircleCandidates,initCentersLength equals 0");
	}
	double * initCenters; //initCentersLength x 2
	initCenters = (double*)malloc(sizeof(double)*initCentersLength*2); 
	//将记录在链表里面的分区后的圆心均值记录到数组里，便于作为初始点进行均值漂移
	for ( i = 0; i<initCentersLength; i++ )// initCenters 大小是 initCentersLength*2
	{
		int addr = 2*i;
		initCenters[addr]   = center_bins[i].x;
		initCenters[addr+1] = center_bins[i].y;
	}
	free((void*)center_bins);//赶紧释放该内存

//	cout<<"2D均值漂移前初始迭代点："<<endl;
//	for (int  i = 0; i<initCentersLength; i++)
//		cout<<initCenters[2*i]<<'\t'<<initCenters[2*i+1]<<endl;
	
	//均值漂移的结果会更新到initCenters里面
	meanShift(points,points_num,2,initCenters,initCentersLength,1,distance_tolerance,1e-6,50);//迭代20次

//	cout<<"2D均值漂移后的聚类中心:"<<endl;
//	for (int  i = 0; i<initCentersLength; i++)
//		cout<<initCenters[2*i]<<'\t'<<initCenters[2*i+1]<<endl;

	//聚类
	//千万要注意centers_num是int型指针，++--时要(*centers_num).
	clusterByDistance(initCenters,initCentersLength,2,distance_tolerance/2,40,centers, centers_num);//此处控制参数要改，要调节

//	cout<<"2D距离聚类，去除重复点后的点集:"<<endl;
//	for (int  i = 0; i<(*centers_num); i++)
//		cout<<centers[2*i]<<'\t'<<centers[2*i+1]<<endl;

	if((*centers_num) <= 0)//可无
	{
		return 0;  //不懂为什么，聚类中心的周围确没有最靠近它的点
		//system("pause");
		//error("cluster2DPoints,(*centers_num)<=0");
	}
	free(initCenters);
	//cout<<"2D聚类后数量:"<<(*centers_num)<<endl;
	return 1;
}

//聚类算法，均值漂移
//三个步骤，一是选取初始迭代点，二是均值漂移，三是去除重复点，从而得到聚类中心
//获得候选圆心的聚类中心(xi,yi)
//输入：
//datas，一维点数据,长度为datas_num x 1
//distance_tolerance,数据点聚类的半径
//输出：
//一维数据点的聚类中心 centers是一维double数组， 大小为 centers_num x 1
//正确返回值为1，出现错误为0. 例如points为空
//切记切记！！！ centers为函数内部申请的内存，用来返回centers_num个点的聚类中心，使用完后一定要释放，记住free(centers)！！！
int  cluster1DDatas( double * datas, int datas_num, double distance_tolerance, double * & centers, int * centers_num)
{
	double rmin,rmax,rdelta;
	int r;
	int i;
	rmin = DBL_MAX;
	rmax = 0;
	for( i  = 0; i < datas_num; i++)//将链表里的r集合复制到数组
	{
		if(datas[i] < rmin)//在这一次遍历中，记录最大最小值
			rmin = datas[i];
		if(datas[i] > rmax)
			rmax = datas[i];
	}
	int nbins_r = 0;
	point1d1i * center_bins;
	rmax += rmin*0.02;//避免rmax-rmin = 0
	rmin -= rmin*0.02;
	rdelta = rmax - rmin;
	nbins_r = (int)ceil((rdelta)/distance_tolerance);
	if(nbins_r <= 0)//至少有一个bin
		nbins_r = 1;
	center_bins = (point1d1i *)malloc(sizeof(point1d1i)*nbins_r);
	memset(center_bins,0,sizeof(point1d1i)*nbins_r);//初始化为0
//	cout<<"1D原始数据:"<<datas_num<<endl;
	for( i = 0; i<datas_num; i++)//对分区vote
	{
//		cout<<datas[i]<<endl;
		r = int((datas[i]-rmin)/rdelta*nbins_r+0.5);
		if(r>=nbins_r)
			r = nbins_r-1;
		center_bins[r].data += datas[i];
		center_bins[r].cnt  ++;			
	}
	int init_r_length = 0;
	for( i = 0; i<nbins_r; i++)
	{
		if(center_bins[i].cnt > 0)//统计非0分区,并且对每一个bin取均值，按照顺序重写到center_bins里面，无内存消耗
		{
			center_bins[init_r_length].data = center_bins[i].data/center_bins[i].cnt;
			init_r_length++;
		}
	}
	if(init_r_length == 0)
	{
		(*centers_num) = 0;
		centers = NULL;
		//cout<<"cluster1DDatas,points number:"<<datas_num<<endl;
		//cout<<"cluster2DDatas,init_r_length equals 0"<<endl;
		return 0;
		//error("generateCircleCandidates,initCentersLength equals 0");
	}
	double * initCenters; //init_r_length x 1
	initCenters = (double*)malloc(sizeof(double)*init_r_length); 
	//将记录在链表里面的分区后的圆心均值记录到数组里，便于作为初始点进行均值漂移
	for ( i = 0; i<init_r_length; i++ )// initCenters 大小是 init_r_length x 1
	{
		initCenters[i] = center_bins[i].data;
	}
	free(center_bins);//赶紧释放该内存

//	cout<<"1D均值漂移前初始迭代点："<<endl;
//	for (int  i = 0; i<init_r_length; i++)
//		cout<<initCenters[i]<<'\t';
//	cout<<endl;

	//至此，得到了均值漂移初始的initCenters，为一维double数组，长度是init_r_length
	meanShift(datas, datas_num, 1, initCenters, init_r_length, 1, distance_tolerance, 1e-6, 20);//迭代20次

//	cout<<"1D均值漂移后的聚类中心:"<<endl;
//	for (int  i = 0; i<init_r_length; i++)
//		cout<<initCenters[i]<<'\t';
//	cout<<endl;

	//聚类
	//千万要注意centers_num是int型指针，++--时要(*centers_num).
	clusterByDistance(initCenters, init_r_length, 1, distance_tolerance/2, 40, centers, centers_num);//控制参数40，最多40个点合成1个点
	
//	cout<<"1D距离聚类，去除重复点后的点集:"<<endl;
//	for (int  i = 0; i<(*centers_num); i++)
//		cout<<centers[i]<<'\t';
//	cout<<endl;

	if((*centers_num) <= 0)//可无
	{
		return 0;  //不懂为什么，聚类中心的周围确没有最靠近它的点
		//system("pause");
		//error("cluster1DDatas,(*centers_num)<=0");
	}
    free(initCenters);
//	cout<<"1D聚类后数量::"<<(*centers_num)<<endl;
	return 1;
}

