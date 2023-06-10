#include "generateEllipseCandidates.h"
#include "candidatesFunctions.cpp"
#include "clusteringFunctions.cpp"
#include "gradientFunctions.cpp"
#include "lsdFunctions.cpp"
#include "lsdInterface.cpp"
#include "miscFunctions.cpp"
#include "nfaFunctions.cpp"
#include "regionsFunctions.cpp"
//======================================MEX function==================================================================

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs!=3) 
      mexErrMsgIdAndTxt( "MATLAB:revord:invalidNumInputs","One input required.");
    else if(nlhs > 4) 
      mexErrMsgIdAndTxt( "MATLAB:revord:maxlhs","Too many output arguments.");
	uchar * inputimg = (uchar*)mxGetData(prhs[0]);
	int imgy,imgx;
	int edge_process_select = (int)mxGetScalar(prhs[1]);//边缘提取选择，1 canny; 2 sobel
	int specified_polarity  = (int)mxGetScalar(prhs[2]);//1,指定检测的椭圆极性要为正; -1指定极性为负; 0表示两种极性椭圆都检测
	imgy = (int)mxGetM(prhs[0]);
	imgx = (int)mxGetN(prhs[0]);
	double *data=(double*)malloc(imgy*imgx*sizeof(double));//将输入矩阵中的图像数据转存到一维数组中
    for(int c=0;c<imgx;c++)
    {
        for(int r=0;r<imgy;r++)
        {
           data[c+r*imgx]=inputimg[r+c*imgy];              
        }    
    }
	int n;//线段数量
	//int new_n;
	vector<vector<int>> groups;
	double * coverages;
	int * reg;
	int reg_x;
	int reg_y;
    double* out=mylsd(&n, data,imgx,imgy,&reg,&reg_x,&reg_y);
	groupLSs(out,n,reg,reg_x,reg_y,&groups);//分组
	free(reg); //释放内存
	calcuGroupCoverage(out,n,groups,coverages);//计算每个组的覆盖角度

    printf("The number of output arc-support line segments: %i\n",n);
	printf("The number of arc-support groups:%i\n",groups.size());
	/*int groups_t = 0;
	for (int i = 0; i<groups.size(); i++)
	{ 
		groups_t+= groups[i].size();
	}
	printf("Groups' total ls num:%i\n",groups_t);*/

	 image_double angles;
	 if(edge_process_select == 1)
		calculateGradient2(data,imgx,imgy,&angles); //version2, sobel; version 3 canny
	 else 
		 calculateGradient3(data,imgx,imgy,&angles); //version2, sobel; version 3 canny
	 PairGroupList * pairGroupList;
	 double distance_tolerance = 2;//max( 2.0, 0.005*min(angles->xsize,angles->ysize) ); // 0.005%*min(xsize,ysize)
	 double * candidates; //候选椭圆
	 double * candidates_out;//输出候选椭圆指针
	 int  candidates_num = 0;//候选椭圆数量
	 //rejectShortLines(out,n,&new_n);
	 pairGroupList = getValidInitialEllipseSet(out,n,&groups,coverages,angles,distance_tolerance,specified_polarity);
	 if(pairGroupList != NULL)
	 {
		printf("The number of initial ellipses：%i \n",pairGroupList->length);
		generateEllipseCandidates(pairGroupList, distance_tolerance, candidates, &candidates_num);
		printf("The number of ellipse candidates: %i \n",candidates_num);
		
		plhs[0] = mxCreateDoubleMatrix(5,candidates_num,mxREAL);
		candidates_out = (double*)mxGetPr(plhs[0]);
		//候选圆组合(xi,yi,ai,bi,phi_i)', 5 x candidates_num, 复制到矩阵candidates_out中
		memcpy(candidates_out,candidates,sizeof(double)*5*candidates_num);

		freePairGroupList(pairGroupList);
		free(candidates);
	 }
	 else
	 {
		 printf("The number of initial ellipses：%i \n",0);
		 double *candidates_out;
		 plhs[0] = mxCreateDoubleMatrix(5,1,mxREAL);
		 candidates_out = (double*)mxGetPr(plhs[0]);
		 candidates_out[0] = candidates_out[1] = candidates_out[2] = candidates_out[3] = candidates_out[4] = 0;
	 }
	 uchar *edgeimg_out;
	 unsigned long edge_pixels_total_num = 0;//边缘总像素
	 double *gradient_vec_out;
	 plhs[1] = mxCreateNumericMatrix(imgy,imgx,mxUINT8_CLASS,mxREAL);
	 edgeimg_out = (uchar*)mxGetData(plhs[1]);
	 //将边缘图复制到矩阵edgeimg_out中
	 //将梯度向量存到矩阵gradient_vec_out中
	 unsigned long addr,g_cnt = 0;
	 for ( int c = 0; c < imgx; c++ )
		 for ( int r = 0; r < imgy; r++)
		 {
			 addr = r*imgx+c;
			 if(angles->data[addr] == NOTDEF)
				 edgeimg_out[c*imgy+r] = 0;
			 else
			 {
				 edgeimg_out[c*imgy+r] = 255;//为边缘点，赋值为白色
				 //------------------------------------------------
				 edge_pixels_total_num++;
			 }
		 }
	 printf("edge pixel number: %i\n",edge_pixels_total_num);
	//申请edge_pixels_total_num x 2 来保存每一个边缘点的梯度向量，以列为优先，符合matlab的习惯
	 plhs[2] = mxCreateDoubleMatrix(2,edge_pixels_total_num,mxREAL);
	 gradient_vec_out = (double*)mxGetPr(plhs[2]);
	  for ( int c = 0; c < imgx; c++ )
		 for ( int r = 0; r < imgy; r++)
		 {
			 addr = r*imgx+c;
			 if(angles->data[addr] != NOTDEF)
			 {
				 gradient_vec_out[g_cnt++] = cos(angles->data[addr]);
				 gradient_vec_out[g_cnt++] = sin(angles->data[addr]);
			 }
		 }
	 //---------------------------------------------------------------------
	//输出线段检测的图像
	if(nlhs == 4)
	{
		Mat ls_mat = Mat::zeros(imgy,imgx,CV_8UC1);
		for ( int i = 0; i<n ; i++)//draw lines
		{
		  Point2d p1(out[8*i],out[8*i+1]),p2(out[8*i+2],out[8*i+3]);
		  line(ls_mat,p1,p2,Scalar(255,0,0));
		}
		if(candidates_num > 0)//draw ellipses
		{
			for ( int i = 0; i<candidates_num; i++)
				ellipse(ls_mat,cv::Point((int)candidates_out[i*5],(int)candidates_out[i*5+1]),cv::Size(candidates_out[i*5+2],candidates_out[i*5+3]),candidates_out[i*5+4]*180/M_PI,0,360,(Scalar(255,0,0)),1);
		}
		plhs[3] = mxCreateDoubleMatrix(imgy,imgx,mxREAL);
		double * ls_img_out = (double*)mxGetPr(plhs[3]);
		//memcpy(ls_out_mat,ls_mat.data ,sizeof(unsigned char)*M*N);
		for (int i = 0; i<imgx; i++)
			for (int j = 0; j<imgy;j++)
				ls_img_out[i*imgy+j]=ls_mat.data[j*imgx+i];
	}
	//---------------------------------------------------------------------
	//这里的free是释放程序中用于产生候选圆所用到的一系列内存
	free(data);
	free(coverages);
	free(out);
	free_image_double(angles);

}










/*
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	int m = mxGetM(prhs[0]);
	int n = mxGetN(prhs[0]);
	double * p = (double*)mxGetData(prhs[0]);
	int sum = 0;
	for (int c = 0; c<n; c++)
		for ( int r = 0; r<m; r++)
			sum += p[c*m+r];
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	double *pout = mxGetPr(plhs[0]);
	*pout = sum;

}
*
