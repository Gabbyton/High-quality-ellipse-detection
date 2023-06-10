#include "generateEllipseCandidates.h"
#include "clusteringFunctions.cpp"
#include "gradientFunctions.cpp"
#include "lsdFunctions.cpp"
#include "lsdInterface.cpp"
#include "mexFunction.cpp"
#include "miscFunctions.cpp"
#include "nfaFunctions.cpp"
#include "regionsFunctions.cpp"
//================================Generate Ellipse Candidates=========================================
//匹配组对，组对的索引参数，椭圆参数
typedef struct PairGroup_s
{
	point2i pairGroupInd;
	point2d center;  //(x0,y0)
	point2d axis;    //(a,b)
	double  phi;     //angle of orientation  
}PairGroup;

//匹配组对节点
typedef struct PairGroupNode_s
{
	point2i pairGroupInd;
	point2d center;  //(x0,y0)
	point2d axis;    //(a,b)
	double  phi;     //angle of orientation  
	PairGroupNode_s* next;
}PairGroupNode;

typedef struct  PairGroupList_s
{
	int length;
	PairGroup * pairGroup;
}PairGroupList;

typedef struct Point2dNode_s
{
	point2d point;
	Point2dNode_s * next;
}Point2dNode;

typedef struct Point3dNode_s
{
	point3d point;
	Point3dNode_s * next;
}Point3dNode;

typedef struct Point5dNode_s
{
	point2d center;
	point2d axis;
	double  phi;
	Point5dNode_s * next;
}Point5dNode;

typedef struct Point1dNode_s
{
	double data;
	Point1dNode_s * next;
}Point1dNode;

PairGroupList * pairGroupListInit( int length)
{
	if(length <= 0)
		error("paired groups length less equal than 0");
	PairGroupList * pairGroupList = (PairGroupList*)malloc(sizeof(PairGroupList));
	pairGroupList->length = length;
	pairGroupList->pairGroup = (PairGroup*)malloc(sizeof(PairGroup)*length);
	if(pairGroupList->pairGroup == NULL)
		error("pairGroupListInit,not enough memory");
	return pairGroupList;
}

void freePairGroupList( PairGroupList * list)
{
	if(list == NULL || list->pairGroup == NULL)
		error("freePairGroupList,invalidate free");
	free(list->pairGroup);
	free(list);
	list->pairGroup = NULL;
	list = NULL;
}

//计算梯度，返回模和角度，同时模值太小的像素点直接抑制掉，赋值为NOTDEF
//mod、angles为了传值，是二级指针
void calculateGradient( double * img_in, unsigned int imgx, unsigned int imgy,image_double * mod, image_double * angles)
{
	if(img_in == NULL || imgx == 0 || imgy == 0)
		error("calculateGradient error!");
	(*mod) = new_image_double(imgx,imgy);
	(*angles) = new_image_double(imgx,imgy);
	double threshold = 2/sin(22.5/180*M_PI);
	unsigned int x,y,adr;
	double com1,com2;
	double gx,gy;
	double norm,norm_square;
	double sum = 0;

	//double max_grad = 0.0;
	//边界初始为NOTDEF
	for ( x = 0; x<imgx; x++) 
	{
		//(*angles)->data[x]=NOTDEF;
		(*angles)->data[(imgy-1)*imgx+x]=NOTDEF;
		//(*mod)->data[x]=NOTDEF;
		(*mod)->data[(imgy-1)*imgx+x]=NOTDEF;
	}
	for ( y = 0; y<imgy; y++) 
	{
		//(*angles)->data[y*imgx] = NOTDEF;
		(*angles)->data[y*imgx+imgx-1] = NOTDEF;
		//(*mod)->data[y*imgx] = NOTDEF;
		(*mod)->data[y*imgx+imgx-1] = NOTDEF;
	}
	 /* compute gradient on the remaining pixels */
	for(x=0;x<imgx-1;x++)
		for(y=0;y<imgy-1;y++)
		{
			adr = y*imgx+x;
		  /*
		     Norm 2 computation using 2x2 pixel window:
		       A B
		       C D
		     and
		       com1 = D-A,  com2 = B-C.
		     Then
		       gx = B+D - (A+C)   horizontal difference
		       gy = C+D - (A+B)   vertical difference
		     com1 and com2 are just to avoid 2 additions.
		   */
		  com1 = img_in[adr+imgx+1] - img_in[adr];
		  com2 = img_in[adr+1]   - img_in[adr+imgx];

		  gx = com1+com2; /* gradient x component */
		  gy = com1-com2; /* gradient y component */
		  norm_square = gx*gx+gy*gy;

		  norm = sqrt( norm_square / 4.0 ); /* gradient norm */

		  (*mod)->data[adr] = norm; /* store gradient norm */

		  if( norm <= threshold ) /* norm too small, gradient no defined */
		  {
		    (*angles)->data[adr] = NOTDEF; /* gradient angle not defined */
			(*mod)->data[adr] = NOTDEF;
		  }
		  else
		    {
		      /* gradient angle computation */
		      (*angles)->data[adr] = atan2(gx,-gy);
		    }
		}
}

void calculateGradient2( double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles)
{
	if(img_in == NULL || imgx == 0 || imgy == 0)
		error("calculateGradient error!");
	image_double mod = new_image_double(imgx,imgy);
	(*angles) = new_image_double(imgx,imgy);
	unsigned int x,y,adr;
	double com1,com2;
	double gx,gy;
	double norm,norm_square;
	double threshold;
	double sum = 0;
	double value;  
	//double max_grad = 0.0;
	//边界初始为NOTDEF
	for ( x = 0; x<imgx; x++) 
	{
		(*angles)->data[x]=NOTDEF;
		(*angles)->data[(imgy-1)*imgx+x]=NOTDEF;
		(mod)->data[x]=NOTDEF;
		(mod)->data[(imgy-1)*imgx+x]=NOTDEF;
	}
	for ( y = 0; y<imgy; y++) 
	{
		(*angles)->data[y*imgx] = NOTDEF;
		(*angles)->data[y*imgx+imgx-1] = NOTDEF;
		(mod)->data[y*imgx] = NOTDEF;
		(mod)->data[y*imgx+imgx-1] = NOTDEF;
	}
	 /* compute gradient on the remaining pixels */
	for(x=1;x<imgx-1;x++)
		for(y=1;y<imgy-1;y++)
		{
			adr = y*imgx+x;
		  /*
		     Norm 2 computation using 2x2 pixel window:
		       A B C
		       D E F
			   G H I
		     and
		       com1 = C-G,  com2 = I-A.
		     Then
		       gx = C+2F+I - (A+2D+G)=com1+com2+2(F-D)   horizontal difference
		       gy = G+2H+I - (A+2B+C)=-com1+com2+2(H-B)   vertical difference
		     com1 and com2 are just to avoid 2 additions.
		   */
		  com1 = img_in[adr-imgx+1] - img_in[adr+imgx-1];
		  com2 = img_in[adr+imgx+1] - img_in[adr-imgx-1];

		  gx = (com1+com2+2*(img_in[adr+1] - img_in[adr-1]))/(8.0*255); /* gradient x component */
		  gy = (-com1+com2+2*(img_in[adr+imgx] - img_in[adr-imgx]))/(8.0*255); /* gradient y component */
		  norm_square = gx*gx+gy*gy;
		  sum+=norm_square;

		  norm = sqrt( norm_square); /* gradient norm */

		  (mod)->data[adr] = norm; /* store gradient norm */
		   /* gradient angle computation */
	     (*angles)->data[adr] = atan2(gy,gx);
		}
	threshold = 2*sqrt(sum/(imgx*imgy));//自动阈值
	//non maximum suppression
	for(x=1;x<imgx-1;x++)
		for(y=1;y<imgy-1;y++)
		{
			adr = y*imgx+x;
			value = (*angles)->data[adr];
			if((mod)->data[adr] < threshold )
			{
				(*angles)->data[adr] = NOTDEF;
				continue;
			}
			if( (value > -M_1_8_PI && value<=M_1_8_PI) || (value <= -M_7_8_PI ) || (value > M_7_8_PI))
			{
				if((mod)->data[adr] <= (mod)->data[adr+1] || (mod)->data[adr] <= (mod)->data[adr-1])
					(*angles)->data[adr] = NOTDEF;
			}
			else if( (value> M_1_8_PI && value<= M_3_8_PI) || (value> -M_7_8_PI && value<= -M_5_8_PI) )
			{
				if((mod)->data[adr] <= (mod)->data[adr-imgx-1] || (mod)->data[adr] <= (mod)->data[adr+imgx+1])
					(*angles)->data[adr] = NOTDEF;
			}
			else if((value> M_3_8_PI && value<= M_5_8_PI) || (value> -M_5_8_PI && value<= -M_3_8_PI))
			{
				if((mod)->data[adr] <= (mod)->data[adr-imgx] || (mod)->data[adr] <= (mod)->data[adr+imgx])
					(*angles)->data[adr] = NOTDEF;
			}
			else 
			{
				if((mod)->data[adr] <= (mod)->data[adr-imgx+1] || (mod)->data[adr] <= (mod)->data[adr+imgx-1])
					(*angles)->data[adr] = NOTDEF;
			}
		}
    //也标记到mod图上面
	//for(x=1;x<imgx-1;x++)
	//	for(y=1;y<imgy-1;y++)
	//	{
	//		if((*angles)->data[y*imgx+x] == NOTDEF)
	//			(mod)->data[y*imgx+x] = NOTDEF;
	//	}
		free_image_double(mod);
}

//=============================================================================
//需要包含如下头文件
//#include <opencv2\opencv.hpp>
//using namespace cv;
void cvCanny3(	const void* srcarr, void* dstarr,
				void* dxarr, void* dyarr,
                int aperture_size )
{
    //cv::Ptr<CvMat> dx, dy;
    cv::AutoBuffer<char> buffer;
    std::vector<uchar*> stack;
    uchar **stack_top = 0, **stack_bottom = 0;

    CvMat srcstub, *src = cvGetMat( srcarr, &srcstub );
    CvMat dststub, *dst = cvGetMat( dstarr, &dststub );

	CvMat dxstub, *dx = cvGetMat( dxarr, &dxstub );
	CvMat dystub, *dy = cvGetMat( dyarr, &dystub );


    CvSize size;
    int flags = aperture_size;
    int low, high;
    int* mag_buf[3];
    uchar* map;
    ptrdiff_t mapstep;
    int maxsize;
    int i, j;
    CvMat mag_row;

    if( CV_MAT_TYPE( src->type ) != CV_8UC1 ||
        CV_MAT_TYPE( dst->type ) != CV_8UC1 ||
		CV_MAT_TYPE( dx->type  ) != CV_16SC1 ||
		CV_MAT_TYPE( dy->type  ) != CV_16SC1 )
        CV_Error( CV_StsUnsupportedFormat, "" );

    if( !CV_ARE_SIZES_EQ( src, dst ))
        CV_Error( CV_StsUnmatchedSizes, "" );
	
    aperture_size &= INT_MAX;
    if( (aperture_size & 1) == 0 || aperture_size < 3 || aperture_size > 7 )
        CV_Error( CV_StsBadFlag, "" );


	size.width = src->cols;
    size.height = src->rows;

	//aperture_size = -1; //SCHARR
    cvSobel( src, dx, 1, 0, aperture_size );
    cvSobel( src, dy, 0, 1, aperture_size );

	Mat1f magGrad(size.height, size.width, 0.f);
	float maxGrad(0);
	float val(0);
	for(i=0; i<size.height; ++i)
	{
		float* _pmag = magGrad.ptr<float>(i);
		const short* _dx = (short*)(dx->data.ptr + dx->step*i);
        const short* _dy = (short*)(dy->data.ptr + dy->step*i);
		for(j=0; j<size.width; ++j)
		{
			val = float(abs(_dx[j]) + abs(_dy[j]));
			_pmag[j] = val;
			maxGrad = (val > maxGrad) ? val : maxGrad;
		}
	}
	
	//% Normalize for threshold selection
	//normalize(magGrad, magGrad, 0.0, 1.0, NORM_MINMAX);

	//% Determine Hysteresis Thresholds
	
	//set magic numbers
	const int NUM_BINS = 64;	
	const double percent_of_pixels_not_edges = 0.9;
	const double threshold_ratio = 0.3;

	//compute histogram
	int bin_size = cvFloor(maxGrad / float(NUM_BINS) + 0.5f) + 1;
	if (bin_size < 1) bin_size = 1;
	int bins[NUM_BINS] = { 0 }; 
	for (i=0; i<size.height; ++i) 
	{
		float *_pmag = magGrad.ptr<float>(i);
		for(j=0; j<size.width; ++j)
		{
			int hgf = int(_pmag[j]);
			bins[int(_pmag[j]) / bin_size]++;
		}
	}	

	
	

	//% Select the thresholds
	float total(0.f);	
	float target = float(size.height * size.width * percent_of_pixels_not_edges);
	int low_thresh, high_thresh(0);
	
	while(total < target)
	{
		total+= bins[high_thresh];
		high_thresh++;
	}
	high_thresh *= bin_size;
	low_thresh = cvFloor(threshold_ratio * float(high_thresh));
	
    if( flags & CV_CANNY_L2_GRADIENT )
    {
        Cv32suf ul, uh;
        ul.f = (float)low_thresh;
        uh.f = (float)high_thresh;

        low = ul.i;
        high = uh.i;
    }
    else
    {
        low = cvFloor( low_thresh );
        high = cvFloor( high_thresh );
    }

    
	buffer.allocate( (size.width+2)*(size.height+2) + (size.width+2)*3*sizeof(int) );
    mag_buf[0] = (int*)(char*)buffer;
    mag_buf[1] = mag_buf[0] + size.width + 2;
    mag_buf[2] = mag_buf[1] + size.width + 2;
    map = (uchar*)(mag_buf[2] + size.width + 2);
    mapstep = size.width + 2;

    maxsize = MAX( 1 << 10, size.width*size.height/10 );
    stack.resize( maxsize );
    stack_top = stack_bottom = &stack[0];

    memset( mag_buf[0], 0, (size.width+2)*sizeof(int) );
    memset( map, 1, mapstep );
    memset( map + mapstep*(size.height + 1), 1, mapstep );

    /* sector numbers
       (Top-Left Origin)

        1   2   3
         *  *  *
          * * *
        0*******0
          * * *
         *  *  *
        3   2   1
    */

    #define CANNY_PUSH(d)    *(d) = (uchar)2, *stack_top++ = (d)
    #define CANNY_POP(d)     (d) = *--stack_top

    mag_row = cvMat( 1, size.width, CV_32F );

    // calculate magnitude and angle of gradient, perform non-maxima supression.
    // fill the map with one of the following values:
    //   0 - the pixel might belong to an edge
    //   1 - the pixel can not belong to an edge
    //   2 - the pixel does belong to an edge
    for( i = 0; i <= size.height; i++ )
    {
        int* _mag = mag_buf[(i > 0) + 1] + 1;
        float* _magf = (float*)_mag;
        const short* _dx = (short*)(dx->data.ptr + dx->step*i);
        const short* _dy = (short*)(dy->data.ptr + dy->step*i);
        uchar* _map;
        int x, y;
        ptrdiff_t magstep1, magstep2;
        int prev_flag = 0;

        if( i < size.height )
        {
            _mag[-1] = _mag[size.width] = 0;

            if( !(flags & CV_CANNY_L2_GRADIENT) )
                for( j = 0; j < size.width; j++ )
                    _mag[j] = abs(_dx[j]) + abs(_dy[j]);

            else
            {
                for( j = 0; j < size.width; j++ )
                {
                    x = _dx[j]; y = _dy[j];
                    _magf[j] = (float)std::sqrt((double)x*x + (double)y*y);
                }
            }
        }
        else
            memset( _mag-1, 0, (size.width + 2)*sizeof(int) );

        // at the very beginning we do not have a complete ring
        // buffer of 3 magnitude rows for non-maxima suppression
        if( i == 0 )
            continue;

        _map = map + mapstep*i + 1;
        _map[-1] = _map[size.width] = 1;

        _mag = mag_buf[1] + 1; // take the central row
        _dx = (short*)(dx->data.ptr + dx->step*(i-1));
        _dy = (short*)(dy->data.ptr + dy->step*(i-1));

        magstep1 = mag_buf[2] - mag_buf[1];
        magstep2 = mag_buf[0] - mag_buf[1];

        if( (stack_top - stack_bottom) + size.width > maxsize )
        {
            int sz = (int)(stack_top - stack_bottom);
            maxsize = MAX( maxsize * 3/2, maxsize + 8 );
            stack.resize(maxsize);
            stack_bottom = &stack[0];
            stack_top = stack_bottom + sz;
        }

        for( j = 0; j < size.width; j++ )
        {
            #define CANNY_SHIFT 15
            #define TG22  (int)(0.4142135623730950488016887242097*(1<<CANNY_SHIFT) + 0.5)

            x = _dx[j];
            y = _dy[j];
            int s = x ^ y;
            int m = _mag[j];

            x = abs(x);
            y = abs(y);
            if( m > low )
            {
                int tg22x = x * TG22;
                int tg67x = tg22x + ((x + x) << CANNY_SHIFT);

                y <<= CANNY_SHIFT;

                if( y < tg22x )
                {
                    if( m > _mag[j-1] && m >= _mag[j+1] )
                    {
                        if( m > high && !prev_flag && _map[j-mapstep] != 2 )
                        {
                            CANNY_PUSH( _map + j );
                            prev_flag = 1;
                        }
                        else
                            _map[j] = (uchar)0;
                        continue;
                    }
                }
                else if( y > tg67x )
                {
                    if( m > _mag[j+magstep2] && m >= _mag[j+magstep1] )
                    {
                        if( m > high && !prev_flag && _map[j-mapstep] != 2 )
                        {
                            CANNY_PUSH( _map + j );
                            prev_flag = 1;
                        }
                        else
                            _map[j] = (uchar)0;
                        continue;
                    }
                }
                else
                {
                    s = s < 0 ? -1 : 1;
                    if( m > _mag[j+magstep2-s] && m > _mag[j+magstep1+s] )
                    {
                        if( m > high && !prev_flag && _map[j-mapstep] != 2 )
                        {
                            CANNY_PUSH( _map + j );
                            prev_flag = 1;
                        }
                        else
                            _map[j] = (uchar)0;
                        continue;
                    }
                }
            }
            prev_flag = 0;
            _map[j] = (uchar)1;
        }

        // scroll the ring buffer
        _mag = mag_buf[0];
        mag_buf[0] = mag_buf[1];
        mag_buf[1] = mag_buf[2];
        mag_buf[2] = _mag;
    }

    // now track the edges (hysteresis thresholding)
    while( stack_top > stack_bottom )
    {
        uchar* m;
        if( (stack_top - stack_bottom) + 8 > maxsize )
        {
            int sz = (int)(stack_top - stack_bottom);
            maxsize = MAX( maxsize * 3/2, maxsize + 8 );
            stack.resize(maxsize);
            stack_bottom = &stack[0];
            stack_top = stack_bottom + sz;
        }

        CANNY_POP(m);

        if( !m[-1] )
            CANNY_PUSH( m - 1 );
        if( !m[1] )
            CANNY_PUSH( m + 1 );
        if( !m[-mapstep-1] )
            CANNY_PUSH( m - mapstep - 1 );
        if( !m[-mapstep] )
            CANNY_PUSH( m - mapstep );
        if( !m[-mapstep+1] )
            CANNY_PUSH( m - mapstep + 1 );
        if( !m[mapstep-1] )
            CANNY_PUSH( m + mapstep - 1 );
        if( !m[mapstep] )
            CANNY_PUSH( m + mapstep );
        if( !m[mapstep+1] )
            CANNY_PUSH( m + mapstep + 1 );
    }

    // the final pass, form the final image
    for( i = 0; i < size.height; i++ )
    {
        const uchar* _map = map + mapstep*(i+1) + 1;
        uchar* _dst = dst->data.ptr + dst->step*i;

        for( j = 0; j < size.width; j++ )
		{
            _dst[j] = (uchar)-(_map[j] >> 1);
		}
	}
};

void Canny3(	InputArray image, OutputArray _edges,
				OutputArray _sobel_x, OutputArray _sobel_y,
                int apertureSize, bool L2gradient )
{
    Mat src = image.getMat();
    _edges.create(src.size(), CV_8U);
	_sobel_x.create(src.size(), CV_16S);
	_sobel_y.create(src.size(), CV_16S);


    CvMat c_src = src, c_dst = _edges.getMat();
	CvMat c_dx = _sobel_x.getMat();
	CvMat c_dy = _sobel_y.getMat();


    cvCanny3(	&c_src, &c_dst, 
				&c_dx, &c_dy,
				apertureSize + (L2gradient ? CV_CANNY_L2_GRADIENT : 0));
};

//canny
void calculateGradient3( double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles)
{
	Mat1b edge;
	Mat1s DX,DY;
	Mat1b gray = Mat::zeros(imgy,imgx,CV_8UC1);
	unsigned int x,y,addr;
	(*angles) = new_image_double(imgx,imgy);
	//copy to gray image
	for ( y = 0; y<imgy; y++)
		for ( x = 0; x<imgx; x++)
		{
			addr = y*imgx+x;
			gray.data[addr] = (uchar)(img_in[addr]);
		}
	//canny
   Canny3(gray,edge,DX,DY,3,false);
   for ( y = 0; y<imgy; y++)
   {
	    short* _dx = DX.ptr<short>(y);
		short* _dy = DY.ptr<short>(y);
		uchar* _e = edge.ptr<uchar>(y);
		for ( x = 0; x<imgx; x++)
		{
			if(_e[x] > 0)//0 or 255
			{
				(*angles)->data[y*imgx+x]  = atan2((double)_dy[x],(double)_dx[x]);//calculate gradient 
			}
			else
				(*angles)->data[y*imgx+x] = NOTDEF;
		}
   }
   edge.release();
   DX.release();
   DY.release();
   gray.release();
}


//=============================================================================
/** Convert ellipse from matrix form to common form:
    ellipse = (centrex,centrey,ax,ay,orientation).
 */
int ellipse2Param(double *p,double param[])
{
	// ax^2 + bxy + cy^2 + dx + ey + f = 0 
  double a,b,c,d,e,f;
  double thetarad,cost,sint,cos_squared,sin_squared,cos_sin,Ao,Au,Av,Auu,Avv,tuCentre,tvCentre,wCentre,uCentre,vCentre,Ru,Rv;
  a = p[0];
  b = p[1];
  c = p[2];
  d = p[3];
  e = p[4];
  f = p[5]; 

  thetarad=0.5*atan2(b,a-c); 
  cost=cos(thetarad);
  sint=sin(thetarad);
  sin_squared=sint*sint;
  cos_squared=cost*cost;
  cos_sin=sint*cost;
  Ao=f;
  Au=d*cost+e* sint;
  Av=-d*sint+e* cost;
  Auu=a*cos_squared+c*sin_squared+b*cos_sin;
  Avv=a*sin_squared+c*cos_squared-b*cos_sin;

  if(Auu==0 || Avv==0){ param[0]=0;param[1]=0;param[2]=0;param[3]=0;param[4]=0;return 0;}
  else
    {
      tuCentre=-Au/(2.*Auu);
      tvCentre=-Av/(2.*Avv);
      wCentre=Ao-Auu*tuCentre*tuCentre-Avv*tvCentre*tvCentre;
      uCentre=tuCentre*cost-tvCentre*sint;
      vCentre=tuCentre*sint+tvCentre*cost;
      Ru=-wCentre/Auu;
      Rv=-wCentre/Avv;
 //     if (Ru>0) Ru=pow(Ru,0.5);
 //     else Ru=-pow(-Ru,0.5);
 //     if (Rv>0) Rv=pow(Rv,0.5);
 //     else Rv=-pow(-Rv,0.5);
	  if (Ru <= 0 || Rv <= 0)//长短轴小于0的情况？？？
		  return 0;
	  Ru = sqrt(Ru);
	  Rv = sqrt(Rv);
      param[0]=uCentre;param[1]=vCentre;
      param[2]=Ru;param[3]=Rv;param[4]=thetarad;
	  //会出现Ru < Rv情况，对调一下
	  if(Ru < Rv )
	  {
		  param[2] = Rv;
		  param[3] = Ru;
		  if(thetarad < 0)//调换长短轴，使得第三个参数为长轴，第四个为短轴
			  param[4] += M_1_2_PI;
		  else
			  param[4] -= M_1_2_PI;
		  if(thetarad < - M_1_2_PI)//长轴倾角限定在-pi/2 ~ pi/2，具备唯一性
			  param[4] += M_PI;
		  if(thetarad > M_1_2_PI)
			  param[4] -= M_PI;
	  }
    }
  return 1;
}
//input : (xi,yi)
//output: x0,y0,a,b,phi,ellipara需要事先申请内存
//successfull, return 1; else return 0
int fitEllipse(point2d* dataxy, int datanum, double* ellipara)
{
	double* D = (double*)malloc(datanum*6*sizeof(double));
	double S[36]; 
	double C[36];
	memset(D,0,sizeof(double)*datanum);
	memset(S,0,sizeof(double)*36);
	memset(C,0,sizeof(double)*36);
	for ( int i = 0; i<datanum; i++)
	{
		D[i*6]  = dataxy[i].x*dataxy[i].x;
		D[i*6+1]= dataxy[i].x*dataxy[i].y;
		D[i*6+2]= dataxy[i].y*dataxy[i].y;
		D[i*6+3]= dataxy[i].x;
		D[i*6+4]= dataxy[i].y;
		D[i*6+5]= 1;
	}
	for ( int i = 0; i<6; i++)
		for ( int j = i; j<6; j++)
		{
			//S[i*6+j]
			for ( int k = 0; k<datanum; k++)
				S[i*6+j] += D[k*6+i]*D[k*6+j];
		}
	free(D);//释放内存
	//对称矩阵赋值
	for ( int i = 0; i<6; i++)
		for ( int j = 0; j<i; j++)
			S[i*6+j]=S[j*6+i];
	C[0*6+2] = 2;
	C[1*6+1] = -1;
	C[2*6+0] = 2;
	// eig(S,C) eig(inv(S)*C)
	double alphar[6],alphai[6],beta[6];
	double vl[36] = {0};//此处不用
	double vr[36] = {0};
	char JOBVL = 'N';
	char JOBVR = 'V';
	ptrdiff_t fitN = 6;
	double fitWork[64];
	ptrdiff_t workLen = 64;
	ptrdiff_t info;
	//info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
	//注意S为对称矩阵，故转置后等于本身，变成列优先，S可以不变
	dggev(&JOBVL,&JOBVR,&fitN,S,&fitN,C,&fitN,alphar,alphai,beta,vl,&fitN,vr,&fitN,fitWork,&workLen,&info);
	if(info == 0)
	{
		int index = -1;
		for ( int i = 0; i<6; i++)
			if( (alphar[i]>=-(2.2204460492503131e-014)) && (alphai[i] == 0) && (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
				index = i;//vr[:,i],vr第i列对应的特征向量则为拟合参数
		if(index == -1)//再试一次，放宽对实部>0的约束，放宽到>-0.005
		{
			double temp = -0.005;//这个参数很关键
			for ( int i = 0; i<6; i++)
			if( (alphar[i]>=temp) && (alphai[i] == 0) && (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
			{
				temp = alphar[i];
				index = i;//vr[:,i],vr第i列对应的特征向量则为拟合参数
			}
		}
		if(index != -1)
		{
			//此处借用beta来传递参数
		    //beta[0] = vr[6*0+index];
		    //beta[1] = vr[6*1+index];
		    //beta[2] = vr[6*2+index];
		    //beta[3] = vr[6*3+index];
		    //beta[4] = vr[6*4+index];
		    //beta[5] = vr[6*5+index];
			  beta[0] = vr[6*index+0];
		      beta[1] = vr[6*index+1];
		      beta[2] = vr[6*index+2];
		      beta[3] = vr[6*index+3];
		      beta[4] = vr[6*index+4];
		      beta[5] = vr[6*index+5];
			ellipse2Param(beta,ellipara);//ax^2 + bxy + cy^2 + dx + ey + f = 0, transform to (x0,y0,a,b,phi)
			return 1;
		}
	}
	return 0;
}

//input: dataxy为数据点(xi,yi),总共有datanum个
//output: 拟合矩阵S. 注意：S需要事先申请内存，double S[36].
inline void calcuFitMatrix(point2d* dataxy, int datanum, double * S)
{
	double* D = (double*)malloc(datanum*6*sizeof(double));
	memset(D,0,sizeof(double)*datanum);
	for ( int i = 0; i<datanum; i++)
	{
		D[i*6]  = dataxy[i].x*dataxy[i].x;
		D[i*6+1]= dataxy[i].x*dataxy[i].y;
		D[i*6+2]= dataxy[i].y*dataxy[i].y;
		D[i*6+3]= dataxy[i].x;
		D[i*6+4]= dataxy[i].y;
		D[i*6+5]= 1;
	}
	for ( int i = 0; i<6; i++)
	{
		for ( int j = i; j<6; j++)
		{
			//S[i*6+j]
			for ( int k = 0; k<datanum; k++)
				S[i*6+j] += D[k*6+i]*D[k*6+j];
		}
	}
    free(D);//释放内存
	//对称矩阵赋值
	for ( int i = 0; i<6; i++)
		for ( int j = 0; j<i; j++)
			S[i*6+j]=S[j*6+i];
}
//input: fit matrixes S1,S2. length is 36.
//output: fit matrix S_out. S_out = S1 + S2.
//S_out事先需要申请内存
inline void addFitMatrix(double * S1, double * S2, double * S_out)
{
	int ind;
	for ( int i = 0; i<6; i++ )
		for ( int j = i; j<6; j++)
		{
			ind = i*6+j;
			S_out[ind] = S1[ind]+S2[ind];
		}
	//对称矩阵赋值
	for ( int i = 0; i<6; i++)
		for ( int j = 0; j<i; j++)
			S_out[i*6+j]=S_out[j*6+i];
}
//input : S矩阵，6 x 6 = 36
//output: (A,B,C,D,E,F)且A>0, ellicoeff需要事先申请内存. 当要转换成(x0,y0,a,b,phi)时，则要用
//ellipse2Param(ellicoeff,ellipara); ax^2 + bxy + cy^2 + dx + ey + f = 0, transform to (x0,y0,a,b,phi)
//successfull, return 1; else return 0
int fitEllipse2(double * S, double* ellicoeff)
{
	double C[36];
	memset(C,0,sizeof(double)*36);
	
	C[0*6+2] = 2;
	C[1*6+1] = -1;
	C[2*6+0] = 2;
	// eig(S,C) eig(inv(S)*C)
	double alphar[6],alphai[6],beta[6];
	double vl[36] = {0};//此处不用
	double vr[36] = {0};
	char JOBVL = 'N';
	char JOBVR = 'V';
	ptrdiff_t fitN = 6;
	double fitWork[64];
	ptrdiff_t workLen = 64;
	ptrdiff_t info;
	//info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
	dggev(&JOBVL,&JOBVR,&fitN,S,&fitN,C,&fitN,alphar,alphai,beta,vl,&fitN,vr,&fitN,fitWork,&workLen,&info);
	if(info == 0)
	{
		int index = -1;
		for ( int i = 0; i<6; i++)
			if( (alphar[i]>=-(2.2204460492503131e-014)) && (alphai[i] == 0) && (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
				index = i;//vr[:,i],vr第i列对应的特征向量则为拟合参数
		if(index == -1)//再试一次，放宽对实部>0的约束，放宽到>-0.005
		{
			double temp = -0.005;//这个参数很关键
			for ( int i = 0; i<6; i++)
			if( (alphar[i]>=temp) && (alphai[i] == 0) && (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
			{
				temp = alphar[i];
				index = i;//vr[:,i],vr第i列对应的特征向量则为拟合参数
			}
		}
		if(index != -1)
		{
			//此处借用beta来传递参数
	        if(vr[6*index+0] < 0)//注意列优先
			{
				ellicoeff[0] = -vr[6*index+0]; //-vr[6*0+index];
				ellicoeff[1] = -vr[6*index+1]; //-vr[6*1+index];
				ellicoeff[2] = -vr[6*index+2]; //-vr[6*2+index];
				ellicoeff[3] = -vr[6*index+3]; //-vr[6*3+index];
				ellicoeff[4] = -vr[6*index+4]; //-vr[6*4+index];
				ellicoeff[5] = -vr[6*index+5]; //-vr[6*5+index];
			}
			else
			{
				ellicoeff[0] = vr[6*index+0];//vr[6*0+index];
				ellicoeff[1] = vr[6*index+1];//vr[6*1+index];
				ellicoeff[2] = vr[6*index+2];//vr[6*2+index];
				ellicoeff[3] = vr[6*index+3];//vr[6*3+index];
				ellicoeff[4] = vr[6*index+4];//vr[6*4+index];
				ellicoeff[5] = vr[6*index+5];//vr[6*5+index];
			}
			return 1;
		}
	}
	return 0;
}

//入参：e1 = (x1,y1,a1,b1,phi1), e2 = (x2,y2,a2,b2,phi2)
//输出：相等为1，否则为0
inline bool isEllipseEqual(double * ellipse1, double * ellipse2, double centers_distance_threshold, double semimajor_errorratio, double semiminor_errorratio, double angle_errorratio, double iscircle_ratio)
{
	bool con1 = ( abs(ellipse1[0] - ellipse2[0]) < centers_distance_threshold && abs(ellipse1[1] - ellipse2[1]) < centers_distance_threshold &&
		abs(ellipse1[2] - ellipse2[2])/MAX(ellipse1[2],ellipse2[2]) < semimajor_errorratio && abs(ellipse1[3] - ellipse2[3])/MIN(ellipse1[3],ellipse2[3]) < semiminor_errorratio );
	bool con2 = ( ellipse1[3]/ellipse1[2] >= iscircle_ratio );//0.9 0.85
	bool con3 = ( ellipse2[3]/ellipse2[2] >= iscircle_ratio );
	bool con4 = ( (con2 && con3) || (con2 == false && con3 == false && abs(ellipse1[4]-ellipse2[4])<= angle_errorratio*M_PI) );
	return (con1 && con4);
}

inline bool regionLimitation( point2d point_g1s, point2d g1s_ls_dir, point2d point_g1e, point2d g1e_ls_dir, point2d point_g2s, point2d g2s_ls_dir, point2d point_g2e, point2d g2e_ls_dir, double polarity, double region_limitation_dis_tolerance)
{
	point2d g1m_ls_dir, g2m_ls_dir;
	point2d g1s_arc_dir,g1e_arc_dir,g1m_arc_dir,g2s_arc_dir,g2e_arc_dir,g2m_arc_dir;
	point2d test_vec1,test_vec2,test_vec3; //弧指向圆心的向量和测试向量
	//组的pend<-pstart构成的向量为gim_arc_dir
	double xdelta, ydelta, theta;
	xdelta = point_g1e.x - point_g1s.x;
	ydelta = point_g1e.y - point_g1s.y;
	theta = atan2(ydelta,xdelta);
	g1m_ls_dir.x = cos(theta);
	g1m_ls_dir.y = sin(theta);
	xdelta = point_g2e.x - point_g2s.x;
	ydelta = point_g2e.y - point_g2s.y;
	theta = atan2(ydelta,xdelta);
	g2m_ls_dir.x = cos(theta);
	g2m_ls_dir.y = sin(theta);

	if( polarity == 1)// polarity is equal 1, arc vector = (dy,-dx)
	{
		g1s_arc_dir.x = g1s_ls_dir.y;
		g1s_arc_dir.y = -g1s_ls_dir.x;
		g1e_arc_dir.x = g1e_ls_dir.y;
		g1e_arc_dir.y = -g1e_ls_dir.x;
		g1m_arc_dir.x = g1m_ls_dir.y;
		g1m_arc_dir.y = -g1m_ls_dir.x;
		g2s_arc_dir.x = g2s_ls_dir.y;
		g2s_arc_dir.y = -g2s_ls_dir.x;
		g2e_arc_dir.x = g2e_ls_dir.y;
		g2e_arc_dir.y = -g2e_ls_dir.x;
		g2m_arc_dir.x = g2m_ls_dir.y;
		g2m_arc_dir.y = -g2m_ls_dir.x;
	}
	else// polarity is equal -1, arc vector = (-dy,dx)
	{
		g1s_arc_dir.x = -g1s_ls_dir.y;
		g1s_arc_dir.y = g1s_ls_dir.x;
		g1e_arc_dir.x = -g1e_ls_dir.y;
		g1e_arc_dir.y = g1e_ls_dir.x;
		g1m_arc_dir.x = -g1m_ls_dir.y;
		g1m_arc_dir.y = g1m_ls_dir.x;
		g2s_arc_dir.x = -g2s_ls_dir.y;
		g2s_arc_dir.y = g2s_ls_dir.x;
		g2e_arc_dir.x = -g2e_ls_dir.y;
		g2e_arc_dir.y = g2e_ls_dir.x;
		g2m_arc_dir.x = -g2m_ls_dir.y;
		g2m_arc_dir.y = g2m_ls_dir.x;
	}
	test_vec1.x = (point_g2e.x - point_g1s.x);
	test_vec1.y = (point_g2e.y - point_g1s.y);
	test_vec2.x = (point_g2s.x - point_g1e.x);
	test_vec2.y = (point_g2s.y - point_g1e.y);
	test_vec3.x = (test_vec1.x + test_vec2.x)/2;
	test_vec3.y = (test_vec1.y + test_vec2.y)/2;
	double t1,t2,t3,t4,t5,t6;
	t1 = dotProduct(g1s_arc_dir,test_vec1);
	t2 = dotProduct(g1e_arc_dir,test_vec2);
	t3 = dotProduct(g1m_arc_dir,test_vec3);
	t4 = -dotProduct(g2e_arc_dir,test_vec1);
	t5 = -dotProduct(g2s_arc_dir,test_vec2);
	t6 = -dotProduct(g2m_arc_dir,test_vec3);

	if(  dotProduct(g1s_arc_dir,test_vec1)  >= region_limitation_dis_tolerance && \
		 dotProduct(g1e_arc_dir,test_vec2)  >= region_limitation_dis_tolerance && \
		 dotProduct(g1m_arc_dir,test_vec3)  >= region_limitation_dis_tolerance && \
		-dotProduct(g2e_arc_dir,test_vec1) >= region_limitation_dis_tolerance && \
	    -dotProduct(g2s_arc_dir,test_vec2) >= region_limitation_dis_tolerance && \
		-dotProduct(g2m_arc_dir,test_vec3) >= region_limitation_dis_tolerance
		)
		return TRUE;
	return FALSE;
}

/*
void drawEllipse(Mat img, double * ellipara)
{
  Point peliicenter(ellipara[0],ellipara[1]);
  Size  saxis(ellipara[2],ellipara[3]);
  //Mat ellimat = Mat::zeros(img.rows,img.cols,CV_8UC3);
  //ellimat.setTo(255);
  static int ccc = 0;
  static unsigned int cnt = 0;
  if(cnt % 2 == 0 )
	  ccc = 0;
  else
  {
	  ccc = 255;
	  cout<<cnt/2<<'\t'<<ellipara[0]<<'\t'<<ellipara[1]<<"\t"<<ellipara[2]<<'\t'<<ellipara[3]<<'\t'<<ellipara[4]<<endl;
  }
  cnt++;

  Mat imgtemp = img.clone();
  ellipse(imgtemp,peliicenter,saxis,ellipara[4]*180/M_PI,0,360,(Scalar(0,255,ccc)),2);
  namedWindow("w1");
  imshow("w1",imgtemp);
  //waitKey(0);
}
void drawEdge(Mat img, point2d * dataxy, int num)
{
	 static int ccc = 0;
     static int cnt = 0;
     cnt++;
     if(cnt % 2 == 0 )
	     ccc = 0;
     else
	  ccc = 255;
	Mat imgtemp = img.clone();
	for (int i = 0; i<num; i++)
	{
		imgtemp.at<Vec3b>(dataxy[i].y,dataxy[i].x) = (Vec3b(ccc,255,0));
	}
	namedWindow("w2");
    imshow("w2",imgtemp);
}
*/

/*----------------------------------------------------------------------------*/
/** Approximate the distance between a point and an ellipse using Rosin distance.
 */
inline double d_rosin (double *param, double x, double y)
{ 
  double ae2 = param[2]*param[2];
  double be2 = param[3]*param[3];
  x = x - param[0];
  y = y - param[1];
  double xp = x*cos(-param[4])-y*sin(-param[4]);
  double yp = x*sin(-param[4])+y*cos(-param[4]);
  double fe2;
  fe2 = ae2-be2;
  double X = xp*xp;
  double Y = yp*yp;
  double delta = (X+Y+fe2)*(X+Y+fe2)-4*X*fe2;
  double A = (X + Y + fe2 - sqrt(delta))/2.0; 
  double ah = sqrt(A);
  double bh2 = fe2-A;
  double term = (A*be2+ae2*bh2);
  double xi = ah*sqrt(ae2*(be2+bh2)/term);
  double yi = param[3]*sqrt(bh2*(ae2-A)/term);
  double d[4],dmin;


  d[0] = dist(xp,yp,xi,yi);
  d[1] = dist(xp,yp,xi,-yi);
  d[2] = dist(xp,yp,-xi,yi);
  d[3] = dist(xp,yp,-xi,-yi);
  dmin = DBL_MAX;
  for ( int i = 0; i<4; i++)
  {
	  if( d[i] <= dmin)
		  dmin = d[i];
  }
//  if (X+Y>xi*xi+yi*yi)
//    return dmin;
//  else return -dmin; 
  return dmin;
}
/*----------------------------------------------------------------------------*/

//输入
//lsd算法检测得到的线段集合lines的数量line_num，return的返回值是line_nums条线段，为一维double型数组lines，长度为8*n，每8个为一组
//存着x1,y1,x2,y2,dx,dy,length,polarity
//groups: 线段分组，每个组存按照几何分布顺序顺时针或者逆时针存储着线段索引，线段索引范围是0~line_num-1. 这里由于是指针，使用时要注意(*group)
//first_group_ind、second_group_ind是匹配组队的索引，当提取salient hypothesis时，second_group_ind = -1, fit_matrix2 = NULL.
//fit_matrix1, fit_matrix2, 分别是组队的对应的拟合矩阵
//angles, 是边缘点图+梯度方向。 无边缘点时是NODEF
//distance_tolerance:
//group_inliers_num:记录着各个组的支持内点数量的数组，实时更新，初始时为0
//输出
//ellipara
bool calcEllipseParametersAndValidate( double * lines, int line_num, vector<vector<int>> * groups, int first_group_ind,int second_group_ind, double * fit_matrix1, double * fit_matrix2, image_double angles, double distance_tolerance, unsigned int * group_inliers_num, point5d *ellipara)
{
	double S[36]; //拟合矩阵S
	double Coefficients[6] = {0,0,0,0,0,0};// ax^2 + bxy + cy^2 + dx + ey + f = 0 
	double param[5], param2[5];
	int info,addr;
	rect rec;
	rect_iter* iter;
	int rec_support_cnt,rec_inliers_cnt;
	bool flag1 = TRUE, flag2 = TRUE;
	double point_normalx, point_normaly, point_normal, temp;
	vector<point2i> first_group_inliers, second_group_inliers;
	point2i pixel_temp;
	double semimajor_errorratio,semiminor_errorratio,iscircle_ratio;
	if( fit_matrix2 == NULL || second_group_ind == -1)//只对一个覆盖度较大的组进行拟合
	{
		for ( int i  = 0; i < 36; i++)
			S[i] = fit_matrix1[i];
	}
	else
	{
		addFitMatrix(fit_matrix1,fit_matrix2,S);//对组对进行拟合， S = fit_matrix1 + fit_matrix2
	}
	info = fitEllipse2(S, Coefficients);// ax^2 + bxy + cy^2 + dx + ey + f = 0, a > 0
	if ( info == 0 )//拟合失败
	{
		ellipara = NULL;
		return FALSE;
	}
	ellipse2Param(Coefficients, param);// (x0,y0,a,b,phi)
	if ( min(param[2],param[3]) < 3*distance_tolerance || max(param[2],param[3]) > min(angles->xsize,angles->ysize) ||  param[0] < 0 || param[0] > angles->xsize || param[1] < 0 || param[1] > angles->ysize )
	{
		ellipara = NULL;
		return FALSE;
	}
	//if ( first_group_ind == 2 && second_group_ind == 8)
	//drawEllipse(img,param);
	//组队中的 first group先进行内点准则验证，并且更新组的支持内点数量
	for ( unsigned int i = 0; i<(*groups)[first_group_ind].size(); i++)
	{
		addr = (*groups)[first_group_ind][i] * 8; //第first_group_ind分组的第i条线段索引*8
		rec.x1 = lines[addr];
		rec.y1 = lines[addr+1];
		rec.x2 = lines[addr+2];
		rec.y2 = lines[addr+3];
		rec.x  = (rec.x1 + rec.x2)/2;
		rec.y  = (rec.y1 + rec.y2)/2;
		rec.dx = lines[addr+4];
		rec.dy = lines[addr+5];
		rec.width = 3*distance_tolerance;
		//line_length[i] = (int)lines[addr+6];//记录线段长度到数组line_length[i]
		rec_support_cnt = rec_inliers_cnt = 0;//清零很重要
		if ( lines[addr+7] == 1) //极性一致
		{
			for(iter = ri_ini(&rec);!ri_end(iter);ri_inc(iter))//线段1
			{
				//外接矩形可能会越界
				if(iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize)
				{
					temp  = angles->data[iter->y*angles->xsize+iter->x] ;//内点的梯度方向
					if(temp!= NOTDEF )
					{
						//test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
						point_normalx = Coefficients[0]*iter->x + (Coefficients[1]*iter->y + Coefficients[3])/2;
						point_normaly = Coefficients[2]*iter->y + (Coefficients[1]*iter->x + Coefficients[4])/2;
						point_normal = atan2(-point_normaly,-point_normalx); //边缘点的法线方向,指向椭圆内侧
						rec_inliers_cnt++;
						if(angle_diff(point_normal,temp) <= M_1_8_PI ) //+- 22.5°内 且 || d - r || < 3 dis_t
						{
							rec_support_cnt++;
							pixel_temp.x = iter->x; pixel_temp.y = iter->y;
							first_group_inliers.push_back(pixel_temp);//添加该线段对应的内点
						}
					} 
				}
			}
		}
		else//极性相反
		{
			for(iter = ri_ini(&rec);!ri_end(iter);ri_inc(iter))//线段1
			{
				//外接矩形可能会越界
				if(iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize)
				{
					temp  = angles->data[iter->y*angles->xsize+iter->x] ;//内点的梯度方向
					if(temp!= NOTDEF )
					{
						//test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
						point_normalx = Coefficients[0]*iter->x + (Coefficients[1]*iter->y + Coefficients[3])/2;
						point_normaly = Coefficients[2]*iter->y + (Coefficients[1]*iter->x + Coefficients[4])/2;
						point_normal = atan2(point_normaly,point_normalx); //边缘点的法线方向,指向椭圆外侧
						rec_inliers_cnt++;
						if(angle_diff(point_normal,temp) <= M_1_8_PI ) //+- 22.5°内 且 || d - r || < 3 dis_t
						{
							rec_support_cnt++;
							pixel_temp.x = iter->x; pixel_temp.y = iter->y;
							first_group_inliers.push_back(pixel_temp);//添加该线段对应的内点
						}
					} 
				}
			}
		}
		if( !( rec_support_cnt > 0 && ( rec_support_cnt >= 0.8*lines[addr+6] || rec_support_cnt*1.0/rec_inliers_cnt >= 0.6) ) )
		{
			flag1 = FALSE; //flag1 初始化时为TRUE, 一旦组内有一条线段不满足要求，直接false, 内点准则验证不通过
			break;
		}
	}
	if ( flag1 == TRUE && first_group_inliers.size() >= 0.8*group_inliers_num[first_group_ind] )//靠近最大统计过的内点,通过验证
	{
		if( first_group_inliers.size() >= group_inliers_num[first_group_ind])//更新组出现过的最大内点数
			group_inliers_num[first_group_ind] =  first_group_inliers.size();
	}
	else 
		flag1 = FALSE;
	//第一个组完成验证
	if ( second_group_ind == -1 || fit_matrix2 == NULL)//只对一个覆盖度较大的组进行拟合
	{
		ellipara->x = param[0];//因为无论如何，都需要返回显著性强的椭圆
	    ellipara->y = param[1];
	    ellipara->a = param[2];
	    ellipara->b = param[3];
	    ellipara->phi = param[4];
		if ( flag1 == TRUE)//通过内点再次拟合，提高质量
		{
			point2d * dataxy = (point2d*)malloc(sizeof(point2d)*first_group_inliers.size());
			for ( unsigned int i = 0; i<first_group_inliers.size(); i++)
			{
				dataxy[i].x = first_group_inliers[i].x;
				dataxy[i].y = first_group_inliers[i].y;
			}
			info = fitEllipse(dataxy,first_group_inliers.size(), param2);
			free(dataxy); //释放内存
			if ( info == 1  && isEllipseEqual(param2,param,3*distance_tolerance,0.1,0.1,0.1,0.9) )
			{
				ellipara->x = param2[0];//更新椭圆，提高品质
			    ellipara->y = param2[1];
			    ellipara->a = param2[2];
			    ellipara->b = param2[3];
			    ellipara->phi = param2[4];
			    //drawEllipse(img,param2);
			}
		}
		return TRUE;//对于只有一个组的提取椭圆，此时直接返回
	}
	//接下来，对组队中的 second group进行内点准则验证，并且更新组的支持内点数量
	if (flag1 == FALSE)//在组队运算中，如果第一个组都无法满足内点要求，直接返回false
		return FALSE;
	for ( unsigned int i = 0; i<(*groups)[second_group_ind].size(); i++)
	{
		addr = (*groups)[second_group_ind][i] * 8; //第first_group_ind分组的第i条线段索引*8
		rec.x1 = lines[addr];
		rec.y1 = lines[addr+1];
		rec.x2 = lines[addr+2];
		rec.y2 = lines[addr+3];
		rec.x  = (rec.x1 + rec.x2)/2;
		rec.y  = (rec.y1 + rec.y2)/2;
		rec.dx = lines[addr+4];
		rec.dy = lines[addr+5];
		rec.width = 3*distance_tolerance;
		//line_length[i] = (int)lines[addr+6];//记录线段长度到数组line_length[i]
		rec_support_cnt = rec_inliers_cnt = 0;//清零很重要
		if ( lines[addr+7] == 1) //极性一致
		{
			for(iter = ri_ini(&rec);!ri_end(iter);ri_inc(iter))//线段1
			{
				//外接矩形可能会越界
				if(iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize)
				{
					temp  = angles->data[iter->y*angles->xsize+iter->x] ;//内点的梯度方向
					if(temp!= NOTDEF )
					{
						//test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
						point_normalx = Coefficients[0]*iter->x + (Coefficients[1]*iter->y + Coefficients[3])/2;
						point_normaly = Coefficients[2]*iter->y + (Coefficients[1]*iter->x + Coefficients[4])/2;
						point_normal = atan2(-point_normaly,-point_normalx); //边缘点的法线方向,指向椭圆内侧
						rec_inliers_cnt++;
						if(angle_diff(point_normal,temp) <= M_1_8_PI ) //+- 22.5°内 且 || d - r || < 3 dis_t
						{
							rec_support_cnt++;
							pixel_temp.x = iter->x; pixel_temp.y = iter->y;
							second_group_inliers.push_back(pixel_temp);//添加该线段对应的内点
						}
					} 
				}
			}
		}
		else//极性相反
		{
			for(iter = ri_ini(&rec);!ri_end(iter);ri_inc(iter))//线段1
			{
				//外接矩形可能会越界
				if(iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize)
				{
					temp  = angles->data[iter->y*angles->xsize+iter->x] ;//内点的梯度方向
					if(temp!= NOTDEF )
					{
						//test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
						point_normalx = Coefficients[0]*iter->x + (Coefficients[1]*iter->y + Coefficients[3])/2;
						point_normaly = Coefficients[2]*iter->y + (Coefficients[1]*iter->x + Coefficients[4])/2;
						point_normal = atan2(point_normaly,point_normalx); //边缘点的法线方向,指向椭圆外侧
						rec_inliers_cnt++;
						if(angle_diff(point_normal,temp) <= M_1_8_PI ) //+- 22.5°内 且 || d - r || < 3 dis_t
						{
							rec_support_cnt++;
							pixel_temp.x = iter->x; pixel_temp.y = iter->y;
							second_group_inliers.push_back(pixel_temp);//添加该线段对应的内点
						}
					} 
				}
			}
		}
		if( !(rec_support_cnt > 0 && ( rec_support_cnt >= 0.8*lines[addr+6] || rec_support_cnt*1.0/rec_inliers_cnt >= 0.6) ) )
		{
			flag2 = FALSE; //flag1 初始化时为TRUE, 一旦组内有一条线段不满足要求，直接false, 内点准则验证不通过
			break;
		}
	}
	if ( flag2 == TRUE && second_group_inliers.size() >= 0.8*group_inliers_num[second_group_ind] )//靠近最大统计过的内点,通过验证
	{
		if(second_group_inliers.size() >= group_inliers_num[second_group_ind])//更新组出现过的最大内点数
			group_inliers_num[second_group_ind] = second_group_inliers.size();
	}
	else 
		flag2 = FALSE;
	//第二个组完成验证
	if ( flag1 == TRUE && flag2 == TRUE)
	{
		point2d * dataxy = (point2d*)malloc(sizeof(point2d)*(first_group_inliers.size() + second_group_inliers.size()));
		for ( unsigned int i = 0; i<first_group_inliers.size(); i++)
		{
			dataxy[i].x = first_group_inliers[i].x;
			dataxy[i].y = first_group_inliers[i].y;
		}
		addr = first_group_inliers.size();
		for ( unsigned int i = 0; i<second_group_inliers.size(); i++)//连接两个数组时一定要注意索引范围
		{
			dataxy[addr+i].x = second_group_inliers[i].x;
			dataxy[addr+i].y = second_group_inliers[i].y;
		}
//		drawEdge(img,dataxy,(first_group_inliers.size() + second_group_inliers.size()));
		info = fitEllipse(dataxy,(first_group_inliers.size() + second_group_inliers.size()), param2);
		free(dataxy); //释放内存
		//小长短轴的椭圆需要放宽参数
		if ( param[2] <= 50 )
			semimajor_errorratio = 0.25;
		else if (param[2] <= 100 )
			semimajor_errorratio = 0.15;
		else
			semimajor_errorratio = 0.1;
		if ( param[3] <= 50 )
			semiminor_errorratio = 0.25;
		else if ( param[3] <= 100)
			semiminor_errorratio = 0.15;
		else
			semiminor_errorratio = 0.1;
		if (param[2] <= 50 && param[3] <= 50 )
			iscircle_ratio = 0.75;
		else if (param[2] >= 50 && param[2] <= 100 &&  param[3] >= 50 && param[3] <= 100 )
			iscircle_ratio = 0.85;
		else
			iscircle_ratio = 0.9;
		if ( info == 1  && isEllipseEqual(param2,param,3*distance_tolerance,semimajor_errorratio,semiminor_errorratio,0.1, iscircle_ratio) )
		{
			ellipara->x = param2[0];//更新椭圆，提高品质
		    ellipara->y = param2[1];
		    ellipara->a = param2[2];
		    ellipara->b = param2[3];
		    ellipara->phi = param2[4];
		    //drawEllipse(img,param2);
			return TRUE;
		}
	}
	return FALSE;
}


//输入
//lsd算法检测得到的线段集合lines的数量line_num，return的返回值是line_nums条线段，为一维double型数组lines，长度为8*n，每8个为一组
//存着x1,y1,x2,y2,dx,dy,length,polarity
//groups: 线段分组，每个组存按照几何分布顺序顺时针或者逆时针存储着线段索引，线段索引范围是0~line_num-1
//coverages: 每个分组的角度覆盖范围0~2pi，如果组里只有1条线段，覆盖角度为0。数组长度等于分组的数量。
//angles 存边缘点的梯度方向gradient direction, 无边缘点位NOTDEF
//返回值 PairedGroupList* list 返回的是初始椭圆集合的数组，长度list->length. 
//切记，该内存在函数内申请，用完该函数记得释放内存，调用函数freePairedSegmentList()进行释放

PairGroupList * getValidInitialEllipseSet( double * lines, int line_num, vector<vector<int>> * groups, double * coverages, image_double angles, double distance_tolerance, int specified_polarity)
{
	//加速计算
	//int* lineInliersIndex = (int*)malloc(sizeof(int)*line_num);//如果第i条线段找到了内点，则记录其索引为j = length(supportInliers),即supportInliers.at(j)存着该线段的支持内点,没找到内点的线段对应索引为初始值-1.
    //vector<vector<point2d>> supportInliers;//保存相应线段的支持内点
	//memset(lineInliersIndex,-1,sizeof(int)*line_num);//此处要实践确实可行，对于整数可以初始化为0，-1.对于浮点数则只可以为0.

	PairGroupList * pairGroupList = NULL;
	PairGroupNode *head, *tail;
	int pairlength = 0;
	point2d pointG1s,pointG1e,pointG2s,pointG2e,g1s_ls_dir,g1e_ls_dir,g2s_ls_dir,g2e_ls_dir;
	double polarity;
	point5d ellipara;
    int groupsNum = (*groups).size();//组的数量
	double * fitMatrixes = (double*)malloc(sizeof(double)*groupsNum*36);//定义拟合矩阵S_{6 x 6}. 每个组都有一个拟合矩阵
	unsigned int * supportInliersNum = (unsigned int*)malloc(sizeof(int)*groupsNum);//用于存储每个组曾经最大出现的支持内点数量
	memset(fitMatrixes,0,sizeof(double)*groupsNum*36);
	memset(supportInliersNum, 0, sizeof(unsigned int)*groupsNum);//初始化为0.
	//double distance_tolerance = max( 2.0, 0.005*min(angles->xsize,angles->ysize) ); // 0.005%*min(xsize,ysize)
    int i,j;
	int cnt_temp,ind_start,ind_end;
	bool info;
    
	//实例化拟合矩阵Si
	point2d * dataxy = (point2d*)malloc(sizeof(point2d)*line_num*2);//申请足够大内存, line_num条线段，共有2line_num个端点
	for ( i = 0; i<groupsNum; i++)
	{
		cnt_temp = 0;//千万注意要清0
		for ( j = 0; j<(*groups)[i].size(); j++)
		{
			//每一条线段有2个端点
			dataxy[cnt_temp].x = lines[(*groups)[i][j]*8];
			dataxy[cnt_temp++].y = lines[(*groups)[i][j]*8+1];
			dataxy[cnt_temp].x = lines[(*groups)[i][j]*8+2];
			dataxy[cnt_temp++].y = lines[(*groups)[i][j]*8+3];
		}
		calcuFitMatrix(dataxy,cnt_temp, fitMatrixes+i*36);
	}
	free(dataxy);//释放内存

	head = tail = NULL;//将初始椭圆集合存储到链表中
	//selection of salient elliptic hypothesis
	for ( i = 0; i<groupsNum; i++)
	{
		if(coverages[i] >= M_4_9_PI )//当组的覆盖角度>= 4pi/9 = 80°, 我们认为具有很大的显著性，可直接拟合提取
		{
			//加入极性判断,只提取指定极性的椭圆
			if (specified_polarity == 0 || (lines[(*groups)[i][0]*8+7] == specified_polarity))
			{
				//显著性大的初始椭圆提取，一定会返回TRUE，因此没必要再判断
				info = calcEllipseParametersAndValidate(lines,line_num,groups,i,-1,(fitMatrixes+i*36),NULL,angles,distance_tolerance,supportInliersNum,&ellipara);
				if (info == FALSE) 
				{
					continue;
					error("getValidInitialEllipseSet, selection of salient ellipses failed!");//这种情况会出现？？,跑54.jpg出现该问题
				}
				PairGroupNode * node = (PairGroupNode*)malloc(sizeof(PairGroupNode));
				node->center.x = ellipara.x;
				node->center.y = ellipara.y;
				node->axis.x   = ellipara.a;
				node->axis.y   = ellipara.b;
				node->phi      = ellipara.phi;
				node->pairGroupInd.x = i;
				node->pairGroupInd.y = -1;//无
				if(head != NULL)
				{
					tail->next = node;
					tail = node;
				}
				else
				{
					head = tail = node;
				}
				pairlength++;
			}
		}
	}
    //selection of pair group hypothesis
	for ( i = 0; i<groupsNum-1; i++)
		for ( j = i+1; j<groupsNum; j++)
			{
				//加入极性判断,只提取指定极性的椭圆
			   if (specified_polarity == 0 || (lines[(*groups)[i][0]*8+7] == specified_polarity))
			    {
					//group i 's polarity is the same as group j; and the number of two paired groups should be >= 3.
					if( lines[(*groups)[i][0]*8+7] == lines[(*groups)[j][0]*8+7] && ((*groups)[i].size() + (*groups)[j].size()) >= 3)
					{
						ind_start = (*groups)[i][0];//第i组的最开始一条线段索引
						ind_end   = (*groups)[i][(*groups)[i].size()-1];//第i组的最后一条线段索引
						pointG1s.x = lines[ind_start*8];
						pointG1s.y = lines[ind_start*8+1];
						g1s_ls_dir.x = lines[ind_start*8+4];
						g1s_ls_dir.y = lines[ind_start*8+5];
						pointG1e.x = lines[ind_end*8+2];
						pointG1e.y = lines[ind_end*8+3];
						g1e_ls_dir.x = lines[ind_end*8+4];
						g1e_ls_dir.y = lines[ind_end*8+5];
						ind_start = (*groups)[j][0];//第j组的最开始一条线段索引
						ind_end   = (*groups)[j][(*groups)[j].size()-1];//第j组的最后一条线段索引
						pointG2s.x = lines[ind_start*8];
						pointG2s.y = lines[ind_start*8+1];
						g2s_ls_dir.x = lines[ind_start*8+4];
						g2s_ls_dir.y = lines[ind_start*8+5];
						pointG2e.x = lines[ind_end*8+2];
						pointG2e.y = lines[ind_end*8+3];
						g2e_ls_dir.x = lines[ind_end*8+4];
						g2e_ls_dir.y = lines[ind_end*8+5];
						polarity = lines[ind_start*8+7]; //i,j两组的极性
						if(regionLimitation(pointG1s,g1s_ls_dir,pointG1e,g1e_ls_dir,pointG2s,g2s_ls_dir,pointG2e,g2e_ls_dir,polarity,-3*distance_tolerance))//都在彼此的线性区域内
						{
							//if ( i == 2)
							//	drawPairGroup(img,lines,(*groups),i,j);

							if(calcEllipseParametersAndValidate(lines,line_num,groups,i,j,(fitMatrixes+i*36),(fitMatrixes+j*36),angles,distance_tolerance,supportInliersNum,&ellipara))//二次一般方程线性求解，线段的内点支持比例
							{
								PairGroupNode * node = (PairGroupNode*)malloc(sizeof(PairGroupNode));
								node->center.x = ellipara.x;
								node->center.y = ellipara.y;
								node->axis.x   = ellipara.a;
								node->axis.y   = ellipara.b;
								node->phi      = ellipara.phi;
								node->pairGroupInd.x = i;
								node->pairGroupInd.y = -1;//无
								if(head != NULL)
								{
									tail->next = node;
									tail = node;
								}
								else
								{
									head = tail = node;
								}
								pairlength++;
							}
						}
						
					}
			   }
			}
	if(pairlength > 0)
	{
		PairGroupNode *p;
		p = head;
		pairGroupList = pairGroupListInit(pairlength);
		for( int i = 0; i<pairGroupList->length; i++)
		{
			pairGroupList->pairGroup[i].center.x = p->center.x;
			pairGroupList->pairGroup[i].center.y = p->center.y;
			pairGroupList->pairGroup[i].axis.x = p->axis.x;
			pairGroupList->pairGroup[i].axis.y = p->axis.y;
			pairGroupList->pairGroup[i].phi = p->phi;
			pairGroupList->pairGroup[i].pairGroupInd.x = p->pairGroupInd.x;//记录组对(i,j),由groups中的第i个组和第j个组构成的匹配组产生该有效椭圆参数
			pairGroupList->pairGroup[i].pairGroupInd.y = p->pairGroupInd.y;
			p = p->next;
		}
		tail->next = NULL;
		while (head != NULL)
		{
			p = head;
			head = head->next;
			free(p);
		}
	}
	//supportInliers.resize(0);
	//free(lineInliersIndex);//释放线段内点的索引
	free(supportInliersNum);//释放存储各个组的支持内点数量的数组
	free(fitMatrixes);//释放存储各个组的拟合矩阵
	return pairGroupList;
}


void generateEllipseCandidates( PairGroupList * pairGroupList, double distance_tolerance, double * & ellipse_candidates, int * candidates_num)
{
	if( pairGroupList->length <= 0 )//检测，至少要有1个样本用来产生候选
	{
		ellipse_candidates = NULL;
		(*candidates_num) = 0;
		return;
	}
	double * centers;
	int center_num; //椭圆中心(xi,yi)的聚类数量
	double * phis;
	int phi_num;    //针对每一个椭圆中心(xi,yi)，倾斜角度phi的聚类数量
	double * axises;
	int axis_num;   //针对每一个椭圆中心和倾角(xi,yi,phi),长短半轴(a,b)的聚类数量
	double * bufferXY = (double*)calloc(pairGroupList->length*2,sizeof(double));
	double * bufferPhi = (double*)calloc(pairGroupList->length,sizeof(double));
	double * bufferAB = (double*)calloc(pairGroupList->length*2,sizeof(double));
	point2i * bufferIndexes = (point2i *)calloc(pairGroupList->length,sizeof(point2i));//point[i].x记录第i个分类在bufferXX中的起始索引位置，point[i].y记录第i个分类在bufferXX中的长度
	double  * buffer2AB = (double*)calloc(pairGroupList->length*2,sizeof(double));
	point2i * buffer2Indexes = (point2i *)calloc(pairGroupList->length,sizeof(point2i));//point[i].x记录第i个分类在bufferXX中的起始索引位置，point[i].y记录第i个分类在bufferXX中的长度
	int     * buffer_temp = (int*)calloc(pairGroupList->length,sizeof(int));
	int addr,addr2,info,ind;
	double dis_min,dis_temp;
	if ( bufferXY == NULL || bufferPhi == NULL || bufferAB == NULL || bufferIndexes == NULL ||
		 buffer2AB == NULL || buffer2Indexes == NULL || buffer_temp == NULL
		)
	{
		ellipse_candidates = NULL;
		(*candidates_num) = 0;
		error("generateEllipseCandidates, not enough memory");
	}
	(*candidates_num) = 0; //候选椭圆数量，初始化为0,非常重要
	//copy
	for ( int i = 0; i<pairGroupList->length; i++)
	{
		addr = 2*i;
		bufferXY[addr] = pairGroupList->pairGroup[i].center.x;
		bufferXY[addr+1] = pairGroupList->pairGroup[i].center.y;
	}
	//cluster the ellipses' centers
	info = cluster2DPoints(bufferXY,pairGroupList->length,distance_tolerance,centers,&center_num);
	if( info == 0)
	{
		ellipse_candidates = NULL;
		(*candidates_num) = 0;
		error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic centers");
	}
	//classification,寻找每个点归属的聚类中心
	for ( int i = 0; i<pairGroupList->length; i++)
	{
		dis_min = DBL_MAX;
		ind = -1;
		for ( int j = 0; j<center_num; j++)
		{
			addr = 2*j;
			dis_temp = (pairGroupList->pairGroup[i].center.x - centers[addr])*(pairGroupList->pairGroup[i].center.x - centers[addr]) + (pairGroupList->pairGroup[i].center.y - centers[addr+1])*(pairGroupList->pairGroup[i].center.y - centers[addr+1]);
			if(dis_temp < dis_min)
			{
				dis_min = dis_temp;
				ind = j; //record the nearest center's index
			}
		}
		buffer_temp[i] = ind; //此处借用buffer2来记下第i个初始椭圆对应第ind个椭圆聚类中心
	}
	//将分类结果按顺序存到bufferXY,bufferPhi,bufferAB中，且bufferIndexes[i]存着第i个聚类中心的起始索引位置和长度
	memset(bufferIndexes,0,sizeof(point2i)*pairGroupList->length);
	ind = 0;//清零，样本点起始位置，索引位置是ind*2,分区的基址
	for ( int i = 0; i<center_num; i++)
	{
		bufferIndexes[i].x = ind; 
		for ( int j = 0; j<pairGroupList->length; j++)
		{
			if ( buffer_temp[j] == i)
			{
				addr = ind*2;//切记长短半轴是一组一组寸储的，需要 x 2
				addr2 = bufferIndexes[i].y*2;
				bufferPhi[ind+bufferIndexes[i].y] = pairGroupList->pairGroup[j].phi;
				bufferAB[addr+addr2] = pairGroupList->pairGroup[j].axis.x;
				bufferAB[addr+addr2+1] = pairGroupList->pairGroup[j].axis.y;
				bufferIndexes[i].y++;//第i个聚类中心周围的点数量加1
			}
		}
		if(bufferIndexes[i].y == 0)//聚类中心周围没有靠近的点
		{
			error("generateEllipseCandidates, no XY points near to the clustering center");
		}
		ind += bufferIndexes[i].y;
	}
	//cout<<"2D cluster centers over"<<endl;
	//对每一个椭圆中心的周围的点进行倾角聚类
	//第i个椭圆聚类中心，其邻近点的索引范围是：bufferIndexs[i].x ~ (bufferIndex[i].x + bufferIndex[i].y-1)
	for ( int i = 0; i<center_num; i++)
	{
		

		double * phi_pointer_temp = bufferPhi+bufferIndexes[i].x;//倾角指针
		double * ab_pointer_temp = bufferAB+bufferIndexes[i].x*2;//长短半轴的指针,记住 x 2
		info = cluster1DDatas(phi_pointer_temp, bufferIndexes[i].y, 0.0873, phis, &phi_num);//对phi聚类, pi/180*5 = 0.0873, 5°误差
		if (info == 0) //不懂为什么，聚类中心centers[i]的周围可能没有最靠近它的点,数量bufferIndexes[i].y = 0
		{ 
			//cout<<"generateEllipseCandidates, cluster2DPoints, error in clustering elliptic phis"<<endl;
			continue;
			//error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic phis");
		}
		//classification,寻找每个点归属的聚类中心
		for ( int j = 0; j<bufferIndexes[i].y; j++ )
		{
			dis_min = DBL_MAX;
			ind = -1;
			for ( int k = 0; k<phi_num; k++)
			{
				dis_temp = (*(phi_pointer_temp+j)-phis[k]) * (*(phi_pointer_temp+j)-phis[k]);
				if(dis_temp < dis_min)
				{
					dis_min = dis_temp;
					ind = k;//record the nearest phi's index
				}
			}
			buffer_temp[j] = ind;
		}
		//将分类结果按顺序存储到buffer2AB中，且buffer2Indexes[j].x对应第i个phi的聚类中心起始点，buffer2Indexes[j].y对应数量(长度)
		memset(buffer2Indexes,0,sizeof(point2i)*bufferIndexes[i].y);
		ind = 0;
		for ( int j = 0; j<phi_num; j++)
		{
			buffer2Indexes[j].x = ind;//起始点
			for ( int k = 0; k<bufferIndexes[i].y; k++)
			{
				if ( buffer_temp[k] == j)
				{
					addr = ind*2;
					addr2 = buffer2Indexes[j].y*2;
					buffer2AB[addr+addr2] = *(ab_pointer_temp+k*2);
					buffer2AB[addr+addr2+1] = *(ab_pointer_temp+k*2+1);
					buffer2Indexes[j].y++;//长度加1
				}
			}
			ind += buffer2Indexes[j].y;
		}
		for ( int j = 0; j<phi_num; j++ )
		{
			double * ab_pointer_temp2 = buffer2AB+buffer2Indexes[j].x*2; //长短半轴的指针,记住 x 2
			info = cluster2DPoints(ab_pointer_temp2, buffer2Indexes[j].y, distance_tolerance, axises, &axis_num);
			if (info == 0) //不懂为什么，聚类中心phi_j的周围可能没有最靠近它的点,数量buffer2Indexes[j].y = 0
			{   
				//cout<<"generateEllipseCandidates, cluster2DPoints, error in clustering elliptic axises"<<endl;
				continue;
				//error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic axises");
			}
			//将候选椭圆重写到bufferXY,bufferPhi,bufferAB里面, 候选椭圆数量(*candidates_num)++
			for ( int k = 0; k<axis_num; k++)
			{
				addr = (*candidates_num)*2;
				bufferXY[addr] = centers[i*2];
				bufferXY[addr+1] = centers[i*2+1];
				bufferPhi[(*candidates_num)] = phis[j];
				bufferAB[addr] = axises[k*2];
				bufferAB[addr+1] = axises[k*2+1];
				(*candidates_num)++;
			}
			free(axises);//cluster2DPoints严格要求，用完axises后，需要释放函数内部申请的内存
		}
		free(phis);//cluster1DDatas严格要求，用完phis后，需要释放函数内部申请的内存
	}
	free(centers);//cluster2DPoints严格要求，用完centers后，需要释放函数内部申请的内存
	//释放在函数开头申请的部分内存
	free(buffer_temp); //此处释放出问题
	free(buffer2Indexes);
	free(buffer2AB);
	free(bufferIndexes);
	ellipse_candidates = (double*)malloc(sizeof(double)*(*candidates_num)*5);
	for ( int i = 0; i < (*candidates_num); i++ )
	{
		addr = 2*i;
		ellipse_candidates[i*5]  = bufferXY[addr];
		ellipse_candidates[i*5+1]= bufferXY[addr+1];
		ellipse_candidates[i*5+2]= bufferAB[addr];
		ellipse_candidates[i*5+3]= bufferAB[addr+1];
		ellipse_candidates[i*5+4]= bufferPhi[i];
	}
	//释放在函数开头申请的内存
	free(bufferAB);
	free(bufferPhi);
	free(bufferXY);
	if((*candidates_num)<= 0)
	{
		*candidates_num = 0;
		ellipse_candidates = NULL;
		//cout<<"no any candidates generated!"<<endl;
	}
}







//==========================================END=======================================================================
/**
输入：
prhs[0]: 输入的灰度图像，单通道，大小是imgy x imgx
prhs[1]: 边缘提取选择，1 canny; 2 sobel
prhs[2]: 检测指定的椭圆极性
输出：
plhs[0]: 候选椭圆组合(xi,yi,ai,bi,phi_i)', 5 x m
plhs[1]: 边缘图，大小是imgy x imgx，设边缘点总数为 edgepix_n. 二值化，0 或者 255
plhs[2]: 边缘点的梯度向量矩阵，大小是 2 x edgepix_n, (cos(theta_rad),sin(theta_rad))'...
plhs[3]: 线段图，大小是imgy x imgx 
*/
/*
compile：
mex generateEllipseCandidates.cpp -IF:\OpenCV\opencv2.4.9\build\include -IF:\OpenCV\opencv2.4.9\build\include\opencv -IF:\OpenCV\opencv2.4.9\build\include\opencv2 -LF:\OpenCV\opencv2.4.9\build\x64\vc11\lib -IF:\Matlab\settlein\extern\include -LF:\Matlab\settlein\extern\lib\win64\microsoft -lopencv_core249 -lopencv_highgui249 -lopencv_imgproc249 -llibmwlapack.lib
*/
/
