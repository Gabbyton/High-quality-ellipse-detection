#include "generateEllipseCandidates.h"
#include "candidatesFunctions.cpp"
#include "clusteringFunctions.cpp"
#include "gradientFunctions.cpp"
#include "lsdFunctions.cpp"
#include "mexFunction.cpp"
#include "miscFunctions.cpp"
#include "nfaFunctions.cpp"
#include "regionsFunctions.cpp"
/*----------------------------------------------------------------------------*/
/*-------------------------- Line Segment Detector ---------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** LSD full interface.
 */
double * LineSegmentDetection( int * n_out,
                               double * img, int X, int Y,
                               double scale, double sigma_scale, double quant,
                               double ang_th, double log_eps, double density_th,
                               int n_bins,
                               int ** reg_img, int * reg_x, int * reg_y )
{
  image_double image;
  ntuple_list out = new_ntuple_list(8);
  double * return_value;
  image_double scaled_image,angles,modgrad;
  image_char used;
  image_char pol;  //对于构成圆弧的像素标记极性，如果梯度的方向和弧的方向指向一致，则为SAME_POLE,否则为OPP_POLE,该标记初始是为0
  image_int region = NULL;
  struct coorlist * list_p;
  struct coorlist * list_p_temp;
//  struct coorlist * mem_p;
  struct rect main_rect;//main rect
  struct rect rect_up,rect_down;//divide the rect into 2 rects:rect_up and rect_down
  struct point2i * reg;
  int reg_size,min_reg_size,i;
  unsigned int xsize,ysize;
  double rho,reg_angle,prec,p;
  double log_nfa = -1,logNT;
//  double log_nfa1,log_nfa2;
  int ls_count = 0;                   /* line segments are numbered 1,2,3,... */
  int seed_cnt = 0;
  int refine_cnt = 0;
  int reg_size_toosmall_cnt=0;

  /* check parameters */
  if( img == NULL || X <= 0 || Y <= 0 ) error("invalid image input.");
  if( scale <= 0.0 ) error("'scale' value must be positive.");
  if( sigma_scale <= 0.0 ) error("'sigma_scale' value must be positive.");
  if( quant < 0.0 ) error("'quant' value must be positive.");
  if( ang_th <= 0.0 || ang_th >= 180.0 )
    error("'ang_th' value must be in the range (0,180).");
  if( density_th < 0.0 || density_th > 1.0 )
    error("'density_th' value must be in the range [0,1].");
  if( n_bins <= 0 ) error("'n_bins' value must be positive.");


  /* angle tolerance */
  prec = M_PI * ang_th / 180.0;
  p = ang_th / 180.0;

  rho = quant / sin(prec); /* gradient magnitude threshold */


  /* load and scale image (if necessary) and compute angle at each pixel */
  image = new_image_double_ptr( (unsigned int) X, (unsigned int) Y, img );
  if( scale != 1.0 )
    {
	  //按照scale进行高斯降采样的图像，注意宽高是上取整，设采样后高宽为imgN*imgM
      scaled_image = gaussian_sampler( image, scale, sigma_scale );
	  //返回一张梯度角度顺时针旋转90°后的align角度图angles，如果梯度角度是(gx,gy)->(-gy,gx)，
	  //和梯度的模的图modgrad,然后按照n_bins进行伪排序返回链表的头指针list_p,里面存的是坐标
	  angles = ll_angle( scaled_image, rho, &list_p,&modgrad, (unsigned int) n_bins );
      free_image_double(scaled_image);
    }
  else
    angles = ll_angle( image, rho, &list_p,&modgrad,(unsigned int) n_bins );
  xsize = angles->xsize;//降采样后的图像的x size，宽度imgM
  ysize = angles->ysize;//降采样后的图像的y size，高度imgN

  /* Number of Tests - NT

     The theoretical number of tests is Np.(XY)^(5/2)
     where X and Y are number of columns and rows of the image.
     Np corresponds to the number of angle precisions considered.
     As the procedure 'rect_improve' tests 5 times to halve the
     angle precision, and 5 more times after improving other factors,
     11 different precision values are potentially tested. Thus,
     the number of tests is
       11 * (X*Y)^(5/2)
     whose logarithm value is
       log10(11) + 5/2 * (log10(X) + log10(Y)).
  */
  logNT = 5.0 * ( log10( (double) xsize ) + log10( (double) ysize ) ) / 2.0
          + log10(11.0);
  min_reg_size = (int) (-logNT/log10(p)); /* minimal number of point2is in region that can give a meaningful event，每个矩形区域内align point2i最小数量*/
  /* initialize some structures */
  if( reg_img != NULL && reg_x != NULL && reg_y != NULL ) /* save region data */
    region = new_image_int_ini(angles->xsize,angles->ysize,0);//申请与降采样后图像一样大小的int类型的内存，该内存的作用是将检测到的线段序号标到相应的图像格子里，该部分可有可无
  used = new_image_char_ini(xsize,ysize,NOTUSED);//申请与降采样后图像一样大小的char类型的内存
  pol  = new_image_char_ini(xsize,ysize,NOTDEF_POL);//像素点处的梯度和弧指向的方向的极性标记
  reg = (struct point2i *) calloc( (size_t) (xsize*ysize), sizeof(struct point2i) );
  if( reg == NULL ) error("not enough memory!");

  list_p_temp = list_p;//记录头链表的头指针，后面需要利用该头指针进行内存释放
  /* search for line segments */
  for(; list_p_temp != NULL; list_p_temp = list_p_temp->next )
    if( used->data[ list_p_temp->x + list_p_temp->y * used->xsize ] == NOTUSED &&
        angles->data[ list_p_temp->x + list_p_temp->y * angles->xsize ] != NOTDEF )
       /* there is no risk of double comparison problems here
          because we are only interested in the exact NOTDEF value */
      {
        /* find the region of connected point2i and ~equal angle */
		//reg是长度为imgN*imgM的一维point2i型数组，有足够大的空间存储生长的区域，reg_size是里面存储了数据的数量，记录的是区域的point2i
		//reg_angle是该区域的主方向的double型变量，存的角度是弧度制
		  seed_cnt ++;
        region_grow( list_p_temp->x, list_p_temp->y, angles, reg, &reg_size,&reg_angle, used, prec );

        /* reject small regions */
        if( reg_size < min_reg_size ) 
		{
			reg_size_toosmall_cnt++;
			continue;
		}

        /* construct rectangular approximation for the region */
		//根据生长的区域得到近似外接矩阵的参数，矩形参数包括:起点，终点，方向theta，宽度等
        region2rect(reg,reg_size,modgrad,reg_angle,prec,p,&main_rect);
		if( FALSE == isArcSegment(reg,reg_size,&main_rect,angles,used,pol,prec,p,&rect_up,&rect_down))
			continue;
        /* Check if the rectangle exceeds the minimal density of
           region point2is. If not, try to improve the region.
           The rectangle will be rejected if the final one does
           not fulfill the minimal density condition.
           This is an addition to the original LSD algorithm published in
           "LSD: A Fast Line Segment Detector with a False Detection Control"
           by R. Grompone von Gioi, J. Jakubowicz, J.M. Morel, and G. Randall.
           The original algorithm is obtained with density_th = 0.0.
         */

        //提纯，通过重新生长区域来达到期望的密度阈值 
        if( !refine( reg, &reg_size, modgrad, reg_angle,
                     prec, p, &main_rect, used, angles, density_th ) ) continue;

		refine_cnt++;
        // compute NFA value 
        log_nfa = rect_improve(&main_rect,angles,logNT,log_eps);//通过改善矩形区域以尝试得到期望的nfa值
        if( log_nfa <= log_eps ) //错误控制
			continue;
        // A New Line Segment was found! 
        ++ls_count;  // increase line segment counter 

        //
        //  The gradient was computed with a 2x2 mask, its value corresponds to
        //  point2is with an offset of (0.5,0.5), that should be added to output.
        //  The coordinates origin is at the center of pixel (0,0).
        //
        main_rect.x1 += 0.5; main_rect.y1 += 0.5;
        main_rect.x2 += 0.5; main_rect.y2 += 0.5;

        // scale the result values if a subsampling was performed */
        if( scale != 1.0 )
          {
            main_rect.x1 /= scale; main_rect.y1 /= scale;
            main_rect.x2 /= scale; main_rect.y2 /= scale;
          //  main_rect.width /= scale;
          }

        /* add line segment found to output */
		add_8tuple( out, main_rect.x1, main_rect.y1, main_rect.x2, main_rect.y2,main_rect.dx,main_rect.dy,
			dist(main_rect.x1, main_rect.y1, main_rect.x2, main_rect.y2), main_rect.polarity);

		//------------------------------------------------------------------------------------------------- 
		/*
		cout<<ls_count<<'\t'<<main_rect.theta<<'\t'<<main_rect.theta*180/M_PI<<"\t polarity:"<<main_rect.polarity<<endl;//打印theta
		
			fstream file1,file2;
			if(ls_count == 1)//清空内容
			{
				file1.open("D:\\Graduate Design\\picture\\sp\\coor.txt",ios::out | ios::trunc);
				file1.close();
				file2.open("D:\\Graduate Design\\picture\\sp\\reg.txt",ios::out | ios::trunc);
				file2.close();
			}
			
			file1.open("D:\\Graduate Design\\picture\\sp\\coor.txt",ios::app);
			file1<<main_rect.x1<<'\t'<<main_rect.y1<<'\t'<<main_rect.x2<<'\t'<<main_rect.y2<<'\t'<<(main_rect.theta*180/M_PI)<<endl;
			file1.close();
			
			if(ls_count == 1)//保持第1根线段的区域
			{
				file2.open("D:\\Graduate Design\\picture\\sp\\reg.txt",ios::app);
				for(i=0; i<reg_size; i++)
					file2<<angles->data[ reg[i].x + reg[i].y * angles->xsize ]*180/M_PI<<endl;
				file2.close();
			}
			*/
		//-------------------------------------------------------------------------------------------------------
        /* add region number to 'region' image if needed */ //将检测到的线段序号标到相应的图像格子里，该部分可有可无
        if( region != NULL )
          for(i=0; i<reg_size; i++)
            region->data[ reg[i].x + reg[i].y * region->xsize ] = ls_count;
      }


  /* free memory */
  free( (void *) image );   /* only the double_image structure should be freed,
                               the data point2ier was provided to this functions
                               and should not be destroyed.                 */
  free_image_double(angles);
  free_image_double(modgrad);
  free_image_char(used);
  free_image_char(pol);
  free( (void *) reg );
//  free( (void *) mem_p );
  //释放分成1024区的存储梯度从大到小的链表,mycode
  //---------------------------------------
  list_p_temp = list_p->next;
  while(list_p_temp != NULL)
  {
	  free(list_p);
	  list_p = list_p_temp;
	  list_p_temp = list_p->next;
  }
  free(list_p);

  //cout<<"seed cnt:"<<seed_cnt<<endl;
  //cout<<"refine cnt:"<<refine_cnt<<endl;
  //cout<<"reg_size_toosmall cnt:"<<reg_size_toosmall_cnt<<endl;
  //----------------------------------------
  /* return the result */
  if( reg_img != NULL && reg_x != NULL && reg_y != NULL )
    {
      if( region == NULL ) error("'region' should be a valid image.");
      *reg_img = region->data;
      if( region->xsize > (unsigned int) INT_MAX ||
          region->xsize > (unsigned int) INT_MAX )
        error("region image to big to fit in INT sizes.");
      *reg_x = (int) (region->xsize);
      *reg_y = (int) (region->ysize);

      /* free the 'region' structure.
         we cannot use the function 'free_image_int' because we need to keep
         the memory with the image data to be returned by this function. */
      free( (void *) region );
    }
  if( out->size > (unsigned int) INT_MAX )
    error("too many detections to fit in an INT.");
  *n_out = (int) (out->size);

  return_value = out->values;
  free( (void *) out );  /* only the 'ntuple_list' structure must be freed,
                            but the 'values' point2ier must be keep to return
                            as a result. */
  return return_value;
}

/*------------------------------------------------------------------------------------------------*/
/**
my code,Alan Lu
输入
img  : 输入图像的一维double型数组,大小为Y*X，按照行优先存储，传入前需要拥有内存
X    : 输入图像的columns
Y    ：输入图像的rows
输出
n_out: lsd算法检测得到的线段的数量n，return的返回值是n条线段，为一维double型数组，长度为8*n，每8个为一组，存着x1,y1,x2,y2,dx,dy,width,polarity
reg_img: 输出标记区域，是一维的int型数组，大小reg_y*reg_x,在相应的像素位置标记着它属于的线段(1,2,3,...n),如果值为0表示不属于任何线段.
         假如外部是int * region_img,则只需要 &region_img,就可以得到标记区域的返回，不需要时直接NULL传入
reg_x  : 输出标记区域的columns,不需要时直接NULL传入
reg_y  : 输出标记区域的rows,不需要时直接NULL传入
*/
double * mylsd(int * n_out, double * img, int X, int Y, int ** reg_img, int * reg_x, int * reg_y)
{
	 /* LSD parameters */
  double scale = 0.8;       /* Scale the image by Gaussian filter to 'scale'. */
  double sigma_scale = 0.6; /* Sigma for Gaussian filter is computed as
                                sigma = sigma_scale/scale.                    */
  double quant = 2.0;       /* Bound to the quantization error on the
                                gradient norm.                                */
  double ang_th = 22.5;     /* Gradient angle tolerance in degrees.           */
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  double density_th = 0.7;  /* Minimal density of region point2is in rectangle. */
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */ 

  return LineSegmentDetection( n_out, img, X, Y, scale, sigma_scale, quant,
                               ang_th, log_eps, density_th, n_bins,
                               reg_img, reg_x, reg_y );
}
//lines: 输入的lines_num条线段，每条线段8个值，存着x1,y1,x2,y2,dx,dy,width,polarity
//lines_num:
//new_lines_num: 拒绝短线段后的new_lines_num条线段，存在lines的前面，而短的线段会放到尾巴处
//此处长度限制参数很重要：目前取8^2, 14^2
void     rejectShortLines(double * lines, int lines_num, int * new_lines_num )
{
	int    new_num = 0;
	int    shor_lines_num = 0;
	double temp;
	new_num = lines_num - shor_lines_num;
	for ( int i = 0; i< new_num; i++)
	{
		if( lines[i*8+6] < 10)//reject short lines, the length threshold is important: 8,14 最后需要调节
		{
			for ( int j = 0; j<8; j++)
			{
				temp = lines[i*8+j];
				lines[i*8+j] = lines[(new_num-1)*8+j];
				lines[(new_num-1)*8+j] = temp;
			}
			i--; //调换后需要检查调换来的线段长度，需要回退
			shor_lines_num++;
			new_num = lines_num - shor_lines_num;
		}
	}
	*new_lines_num = new_num;
}

/*----------------------------------------------------------------------------*/
//输入：
//start_angle,end_angle, 角度方位是(-pi,pi).  
//  pi    ------->x  0
//        |
//        |
//       y\/ pi/2
//polarity: 当polarity为1时，表示的是从start_angle按照逆时针方向旋转到end_angle的角度;当polarity为-1时，表示的是从start_angle按照顺时针方向旋转到end_angle的角度;
//返回值： 旋转角度coverage
inline double rotateAngle(double start_angle, double end_angle, int polarity)
{
	double coverage;
	//首先需要将angle1和angle2转换到 0 ~ 2pi
	if(start_angle < 0) start_angle += M_2__PI;//限制角度在0~2pi之间
	if(end_angle < 0 ) end_angle += M_2__PI;
	if(polarity == 1)//极性为1
	{
		coverage = start_angle - end_angle;
	}
	else //极性为-1
	{ 
		coverage = end_angle - start_angle;
	}
	if(coverage < 0) coverage += M_2__PI;
	return coverage;
}
//对线段按照凸性和距离进行分组
//lines: 输入的lines_num条线段，每条线段8个值，存着x1,y1,x2,y2,dx,dy,length,polarity
//lines_num:
//输出分组groups. 每个组是一个vector<int>
//注意：切记用完region,需要在函数外面手动释放region
void groupLSs(double *lines, int line_num, int * region, int imgx, int imgy, vector<vector<int>> * groups)
{
	if(line_num == 0)
	{
		groups = NULL;
		return;
	}
	unsigned char isEnd = 0;//是否还可以继续搜寻
	int currentLine; //当前线段
	char * label = (char*)calloc(line_num, sizeof(char));
	memset(label,0,sizeof(char)*line_num); //init the label all to be zero
	int * group_up = (int*)malloc(sizeof(int)*line_num);//申请足够内存，存储延线段方向得到的分组的线段
	int * group_down = (int*)malloc(sizeof(int)*line_num);//存储线段反方向分组的线段
	int group_up_cnt,group_down_cnt;
	//coorlist * head,*tail;
	vector<int> group_temp;
	point2d dir_vec1,dir_vec2;
	point2i *votebin = (point2i*)calloc(line_num,sizeof(point2i));//申请足够内存，用来投票. x记录线段索引，y记录票数
	int bincnt = 0;
	int xx,yy,temp;
	double start_angle,end_angle,angle_delta;
	for ( int i = 0; i<line_num; i++)
	{
		if( label[i] == 0)//未被分组过
		{
			group_up_cnt = group_down_cnt = 0;//每开始寻找一组，需要置零
			//先从第i条线段的头部开始搜索，进行分组,结果存在group_up里面
			group_up[group_up_cnt++] = i;//记录线段i,注意线段是0~line_num-1
			isEnd = 0;//置零，表示还可以从当前线段开始搜索，还未结束
	     	currentLine = i;
			while(isEnd == 0)
			{
				label[currentLine] = 1; //标记该线段已经被分组
				//head = tail = NULL;
		        bincnt = 0;
				dir_vec1.x = lines[currentLine*8+4];
				dir_vec1.y = lines[currentLine*8+5];
				if ( lines[currentLine*8+7] == 1)//极性为正
				{
					//将dir_vec1逆时针旋转45°
					dir_vec2.x = (dir_vec1.x + dir_vec1.y)*0.707106781186548; // sqrt(2)/2 = 0.707106781186548
				    dir_vec2.y = (-dir_vec1.x + dir_vec1.y)*0.707106781186548;
				}
				else
				{
					//将dir_vec1顺时针旋转45°
					dir_vec2.x = (dir_vec1.x - dir_vec1.y)*0.707106781186548; // sqrt(2)/2 = 0.707106781186548
				    dir_vec2.y = (dir_vec1.x + dir_vec1.y)*0.707106781186548;
				}
				for ( int j = 1; j<=4; j++)
					for ( int k = 1; k<=4; k++)//在4x4邻域内搜索
					{
						xx = (int)(lines[currentLine*8+2]*0.8+j*dir_vec1.x+k*dir_vec2.x);
						yy = (int)(lines[currentLine*8+3]*0.8+j*dir_vec1.y+k*dir_vec2.y);
						if(xx < 0 || xx >= imgx || yy < 0 || yy >= imgy)//越界
							continue;
						temp = region[yy*imgx+xx];
						if(temp>0)//表示有线段的支持区域，在1~line_num
						{
							region[yy*imgx+xx] = -temp;//取负数标记
							for (xx = 0; xx<bincnt; xx++)
							{
								if(votebin[xx].x == temp - 1)//如果以前投票过，直接在相应的bin的票数上加1
								{
									votebin[xx].y++;
									break;
								}
							}
							if(xx == bincnt)//如果以前没有投票过，增加该线段，并记录票数为1
							{
								if(bincnt == line_num)
									error("group ls error1");
								votebin[bincnt].x = temp - 1;
								votebin[bincnt].y = 1;
								bincnt++; //bin的总数加1
							}
						}
					}
			    //寻找投票最多的线段，并且需要满足数量大于一定值
			    temp = 0;
				for ( int j = 0; j<bincnt; j++)
				{
					if(votebin[j].y>temp)
					{
						temp = votebin[j].y;
						xx = votebin[j].x;//借用xx变量
					}
				}
				if ( temp >= 5 && label[xx] == 0 && lines[8*xx+7] == lines[8*i+7] )//待实验调整参数值
				{
					if(group_up_cnt == line_num)
					   error("group ls error2");
					yy = group_up_cnt-1;//借用yy变量
					start_angle = atan2(lines[8*group_up[yy]+5],lines[8*group_up[yy]+4]);
					end_angle = atan2(lines[8*xx+5],lines[8*xx+4]);
					angle_delta = rotateAngle(start_angle,end_angle,(int)lines[8*i+7]);
					if(angle_delta <= M_3_8_PI)//相邻两线段的旋转夹角也需要满足在pi/4内
					{
						group_up[group_up_cnt++] = xx;//压入线段
						currentLine = xx; //更新当前搜索线段
					}
					else
						isEnd = 1;
				}
				else
					isEnd = 1;//结束，已经找不到可以分组的线段了
			}
			//先从第i条线段的尾部开始搜索，进行分组,结果存在group_down里面。记住，第i条线段在group_up和group_down中的0索引处都储存了
			group_down[group_down_cnt++] = i; 
			isEnd = 0;//置零，表示还可以从当前线段开始搜索，还未结束
	     	currentLine = i;
			while(isEnd == 0)
			{
				label[currentLine] = 1; //标记该线段已经被分组
				//head = tail = NULL;
		        bincnt = 0;
				dir_vec1.x = -lines[currentLine*8+4];
				dir_vec1.y = -lines[currentLine*8+5];
				if ( lines[currentLine*8+7] == 1)//极性相同
				{
					//将dir_vec1顺时针旋转45°
					dir_vec2.x = (dir_vec1.x - dir_vec1.y)*0.707106781186548; // sqrt(2)/2 = 0.707106781186548
				    dir_vec2.y = (dir_vec1.x + dir_vec1.y)*0.707106781186548;
				}
				else
				{
					//将dir_vec1顺时针旋转45°
					dir_vec2.x = (dir_vec1.x + dir_vec1.y)*0.707106781186548; // sqrt(2)/2 = 0.707106781186548
				    dir_vec2.y = (-dir_vec1.x + dir_vec1.y)*0.707106781186548;
				}
				for ( int j = 1; j<=4; j++)
					for ( int k = 1; k<=4; k++)//在4x4邻域内搜索
					{
						xx = (int)(lines[currentLine*8+0]*0.8+j*dir_vec1.x+k*dir_vec2.x);
						yy = (int)(lines[currentLine*8+1]*0.8+j*dir_vec1.y+k*dir_vec2.y);
						if(xx < 0 || xx >= imgx || yy < 0 || yy >= imgy)//越界
							continue;
						temp = region[yy*imgx+xx];
						if(temp>0)//表示有线段的支持区域，在1~line_num
						{
							region[yy*imgx+xx] = -temp;//取负数标记
							for (xx = 0; xx<bincnt; xx++)
							{
								if(votebin[xx].x == temp - 1)//如果以前投票过，直接在相应的bin的票数上加1
								{
									votebin[xx].y++;
									break;
								}
							}
							if(xx == bincnt)//如果以前没有投票过，增加该线段，并记录票数为1
							{
								if(bincnt == line_num)
									error("group ls error3");
								votebin[bincnt].x = temp - 1;
								votebin[bincnt].y = 1;
								bincnt++; //bin的总数加1
							}
						}
					}
			    //寻找投票最多的线段，并且需要满足数量大于一定值
			    temp = 0;
				for ( int j = 0; j<bincnt; j++)
				{
					if(votebin[j].y>temp)
					{
						temp = votebin[j].y;
						xx = votebin[j].x;//借用xx变量
					}
				}
				if ( temp >= 5 && label[xx] == 0 && lines[8*xx+7] == lines[8*i+7])//待实验调整参数值
				{
					if(group_down_cnt == line_num)
					   error("group ls error2");
					yy = group_down_cnt-1;//借用yy变量
					start_angle = atan2(lines[8*group_down[yy]+5],lines[8*group_down[yy]+4]);
					end_angle = atan2(lines[8*xx+5],lines[8*xx+4]);
					angle_delta = rotateAngle(end_angle,start_angle,(int)lines[8*i+7]);//注意此时需要调换一下，因为是从尾部开始搜索
					if(angle_delta < M_3_8_PI)//相邻两线段的旋转夹角也需要满足在pi/4内,pi*3/8 = 66.5°
					{
						group_down[group_down_cnt++] = xx; //压入线段
						currentLine = xx; //更新当前搜索线段
					}
					else
						isEnd = 1;
				}
				else
					isEnd = 1;//结束，已经找不到可以分组的线段了
			}
			(*groups).push_back(group_temp); //添加线段分组
			temp = (*groups).size()-1;
			for (int j = group_down_cnt-1; j>= 0; j--)
			{
				(*groups)[temp].push_back(group_down[j]);
			}
			for (int j = 1; j<group_up_cnt; j++)//由于第i条线段在group_up和group_down都储存了，所以就从索引1开始
			{
				(*groups)[temp].push_back(group_up[j]);
			}
		}
	}
	free(label);
	free(group_up);
	free(group_down);
	free(votebin);
}
//计算groups中每个组的跨度
//输入：
//lines: 输入的lines_num条线段，每条线段8个值，存着x1,y1,x2,y2,dx,dy,length,polarity
//lines_num:
//groups: 分组，每个分组都存着线段的索引
//输出:
//coverages: 每个组的跨度，当组内线段只有1条时，跨度为0. coverages的长度等于组的数量 = groups.size()
//注意，coverages用前不需要申请内存，coverages用完后，需要在函数外手动释放内存，长度等于分组数量
void calcuGroupCoverage(double * lines, int line_num, vector<vector<int>> groups, double * &coverages)
{
	int groups_num = groups.size();
	int temp;
	double start_angle,end_angle;
	coverages = (double*)malloc(sizeof(double)*groups_num);
	for ( int i = 0; i<groups_num; i++)
	{
		temp = groups[i].size()-1;
		if(groups[i].size() == 0)//第i个分组只有1条线段，则跨度为0
		{
			coverages[i] = 0;
		}
		else
		{
			start_angle = atan2(lines[8*groups[i][0]+5],lines[8*groups[i][0]+4]);
			end_angle = atan2(lines[8*groups[i][temp]+5],lines[8*groups[i][temp]+4]);
			coverages[i] = rotateAngle(start_angle,end_angle,(int)lines[8*groups[i][0]+7]);
		}
	}
}

