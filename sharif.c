#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

int main( ){

char line[255];
int e,i,j,l,n,m,index;
double x,y;
double E;
//E=29.5e6;
double A;
//A=1; 
 int h;  
double data1,data2;
 int data3,data4,data5;  
      
//coordinate of X and Y direction

/*printf("Enter the NODE number :\n");
scanf("%d",&n);

printf("Enter the ELEMENTS number :\n");
scanf("%d",&m);


printf("Enter the fiexd DOF number :\n");
scanf("%d",&h);

printf("Enter the young modulus :\n");
scanf("%lf",&E);

printf("Enter the each ELEMENTS cross section area :\n");
scanf("%lf",&A);*/


 FILE *fp;
 fp = fopen("truss.txt","r");
                            
    while (strcmp (line,"BASIC DATA START\r\n"))  // coordinate reading
 
fgets(line,255,fp);

     while (strcmp (line,"BASIC DATA END\r\n")) 
      {
      fgets(line,255,fp);
      if (!strcmp(line,"BASIC DATA END\r\n")) break;

sscanf(line," %lf %lf %d %d %d",&data1,&data2,&data3,&data4,&data5);
  A=data1;
  E=data2;
  n=data3;
  m=data4;
  h=data5;

printf("A=%lf \n",A);
printf("E=%lf \n",E);
printf("n=%d \n",n);
printf("m=%d \n",m);
printf("h=%d \n",h);

}

double coor_xy[n][2];// coordintate 
int index1,start,end;
double p[m][2]; // coordintate diffrent between two nodes




    while (strcmp (line,"BEGIN NODE \r\n"))  // node and coordinate reading
 
fgets(line,255,fp);

     while (strcmp (line,"END NODE \r\n")) 
      {
      fgets(line,255,fp);
      if (!strcmp(line,"END NODE \r\n")) break;

sscanf(line,"%d %lf %lf",&index,&x,&y);

coor_xy[index][0]=x;
coor_xy[index][1]=y;

printf("coor_xy[%d][0]=%lf\n",index,coor_xy[index][0]);
printf("coor_xy[%d][1]=%lf\n",index,coor_xy[index][1]);

}



while (strcmp (line,"BEGIN ELEMENTS\r\n"))  // element and node reading
 
fgets(line,255,fp);

     while (strcmp (line,"END ELEMENTS \r\n")) 
      {
      fgets(line,255,fp);

 if (!strcmp(line,"END ELEMENTS \r\n")) break;

sscanf(line,"%d %d %d\n",&index1,&start,&end);

p[index1][0]=coor_xy[end][0]-coor_xy[start][0];// coordintate diffrent between two nodes of element in x direction
p[index1][1]=coor_xy[end][1]-coor_xy[start][1];// coordintate diffrent between two nodes of element in y direction

printf("p[%d][0]=%lf\n",index1,p[index1][0]);
printf("p[%d][1]=%lf\n",index1,p[index1][1]);

printf("p[%d][0]=coor_xy[%d][0]-coor_xy[%d][0]\n",index1,end,start);
printf("p[%d][1]=coor_xy[%d][1]-coor_xy[%d][1]\n",index1,end,start);
}


//Now calculation of the length of 9 elements from 1 to 9 respective
double L[m];
double c[m];
double s[m];
double d[m];
for(i=0;i<m;i++){      //distance between two nodes of element
                  
L[i] = sqrt((p[i][1]*p[i][1])+(p[i][0]*p[i][0])); // sqrt(y2-y1)²+(x2-x1)² formula

c[i]=(p[i][0])/L[i]; // for c=(x2-x1)/L formula 
s[i]=(p[i][1])/L[i];  // for s=(y2-y1)/L formula

d[i]=(E*A)/L[i];       //EA/L calculation 

printf("L[%d]=%lf \n",i,L[i]);
printf("c[%d]=%lf \n",i,c[i]);
printf("s[%d]=%lf \n",i,s[i]);
printf("d[%d]=%lf \n",i,d[i]);
printf("E=%lf \n",E);
printf("A=%lf \n",A);
} 

printf("\n");
printf("stiffness matrix calculation for each element \n");
printf("\n");

double k[i][4][4];  //4x4   element stifness matrix formula for 2D structure

for(i=0;i<m;i++){   // calculation of all element matrix value 

k[i][0][0]=(d[i])*(c[i]*c[i]);    
k[i][0][1]=(d[i])*(c[i]*s[i]);
k[i][0][2]=(d[i])*(-(c[i]*c[i]));   
k[i][0][3]=(d[i])*(-(c[i]*s[i]));
k[i][1][0]=(d[i])*(c[i]*s[i]);
k[i][1][1]=(d[i])*(s[i]*s[i]);
k[i][1][2]=(d[i])*(-(c[i]*s[i]));
k[i][1][3]=(d[i])*(-(s[i]*s[i]));
k[i][2][0]=(d[i])*(-(c[i]*c[i]));    
k[i][2][1]=(d[i])*(-(c[i]*s[i]));
k[i][2][2]=(d[i])*(c[i]*c[i]);    
k[i][2][3]=(d[i])*(c[i]*s[i]);
k[i][3][0]=(d[i])*(-(c[i]*s[i]));
k[i][3][1]=(d[i])*(-(s[i]*s[i]));
k[i][3][2]=(d[i])*(c[i]*s[i]);
k[i][3][3]=(d[i])*(s[i]*s[i]);
}
for(i=0;i<m;i++){
     for(j=0;j<=3;j++){
          for(l=0;l<=3;l++){

  printf("k[%d][%d][%d]=%lf \n",i,j,l,(k[i][j][l]));

         //printf(" \t %lf",(k[i][j][l]));
          //printf("\n");   
                           }
                      }
                 }
printf("\n");
printf("Now Identity matrix begin from file \n");
printf("\n");

//double global[12][12]={0}; //global striffness matrix

double global[2*n][2*n]; //global striffness matrix
memset(global,0,(2*n)*(2*n)*sizeof(double));

int ID[4][m]; // indentity matrix

 FILE *fp1;  // file open for  ID matrix data reading
 fp1 = fopen("ID.txt","r");


 for(i=0;i<4;i++){  //row number  all time 4 for ID matrix because of 2D truss
          for(j=0;j<m;j++){
 fscanf(fp1, "%d",&ID[i][j]); 

}
}
fclose(fp1);

 for(i=0;i<4;i++){ // 
          for(j=0;j<m;j++){
printf("ID[%d][%d]=%d \n",i,j,ID[i][j]);
}
printf("\n");
}

// loop for  element stiffness matrix to global stiffness matrix
printf("\n");
printf("Global stiffness matrix calculation \n");
printf("\n");

for(e=0;e<m;e++){
     for(i=0;i<=3;i++){
          for(j=0;j<=3;j++){

global[ID[i][e]][ID[j][e]]=global[ID[i][e]][ID[j][e]]+k[e][i][j];

                        }
                      }
                 }

/*for(e=0;e<m;e++){
     for(i=0;i<4;i++){
          for(j=0;j<4;j++){

 printf("global[%d][%d]=%lf \n",ID[i][e],ID[j][e],(global[ID[i][e]][ID[j][e]]));
                         }
                      }
                 }
*/
for(i=0;i<(n*2);i++){
          for(j=0;j<(n*2);j++){
printf("global[%d][%d]=%lf \n",i,j,(global[i][j]));         
                        
                     }
                 }
printf("DOF with serial number which is not fixed \n");

// now we will delete row and column for boundary condition by file reading
/*int h;
printf("Enter the fiexd DOF number :\n");
scanf("%d",&h);*/

double intermediate[(n*2)-h][n*2];
double final[(n*2)-h][(n*2)-h];
int pos[(n*2)-h]; // FREE DOF  
int serial3,index3;
//int pos[12-4];


    while (strcmp (line,"FREE DOF BEGIN\r\n"))  // DOF reading
 
fgets(line,255,fp);

     while (strcmp (line,"FREE DOF END \r\n")) 
      {

      fgets(line,255,fp);

      if (!strcmp(line,"FREE DOF END\r\n")) break;                          
sscanf(line,"%d %d \n",&serial3,&index3);

pos[serial3]=index3;

printf("pos[%d]=%d \n",serial3,index3);
}

printf("\n");
printf("Final global matrix after row column eleminate\n");
printf("\n");
/*
pos[0]=0;
pos[1]=1;
pos[2]=4;
pos[3]=5;
pos[4]=6;
pos[5]=7;
pos[6]=8;
pos[7]=9;*/
/*deleting 3rd, 4th, 11th and 12th row*/
for(i=0;i<(n*2)-h;i++){
for (j=0;j<(n*2);j++){
    intermediate[i][j] = global[pos[i]][j];
}
}

/*deleting 3rd, 4th, 10th and 11th column*/
for (i=0;i<(n*2)-h;i++){
for (j=0;j<(n*2)-h;j++){
   final[i][j] = intermediate[i][pos[j]];
     }
}
for ( i = 0 ; i <(n*2)-h ; i++ ){
for ( j =0 ; j <(n*2)-h; j++ ){
printf("final[%d][%d]=%lf \n ",i,j,final[i][j]);
}
}
printf("\n");
/*printf("\n");
printf("final \n");
for ( i = 0 ; i <(n*2)-4  ; i++ ){
for ( j =0 ; j < (n*2)-4 ; j++ ){
printf("%lf  ",final[i][j]);
}
printf("\r\n");
}
*/
//  now solution of system of linear equation Ku=F

int serial5,dof,force1;
double xold[(n*2)-h];
int force[(n*2)-h];

printf("Force boundary condition  \n");
printf("\n");

while (strcmp (line,"FORCE ON DOF\r\n"))  // force boundary condition reading
 
fgets(line,255,fp);

     while (strcmp (line,"END FORCE\r\n")) 
      {

      fgets(line,255,fp);

      if (!strcmp(line,"END FORCE\r\n")) break;                          
sscanf(line,"%d %d %d \n",&serial5,&dof,&force1);

force[serial5]=force1;

printf("force[%d]=%d \n",serial5,force1);
}
fclose(fp);

/*force[0]=0;
force[1]=0;
force[2]=0;
force[3]=-25000;
force[4]=0;
force[5]=0;
force[6]=0;
force[7]=0;*/

//printf("Force boundary condition with DOF \n");

// Now we will solve our system of linear equation by jacobi iteration method
for(i=0;i<(n*2)-h;i++){
    xold[i]=force[i]/final[i][i];
}

    double xnew[(n*2)-h];

    int t;
    double sum;
    //for(i=0;i<n;i++){
        //sum=0;
    //}
    t=0;

    while(t<1000000){

    for(i=0;i<(n*2)-h;i++){
        sum=force[i];
        for(j=0;j<(n*2)-h;j++){
            if(i!=j)
            //continue;
                //sum=0;
            sum=sum-(final[i][j]*xold[j]);

            }
            xnew[i]=sum/final[i][i];
        }
    //for(i=0;i<n;i++){
      //  xnew[i]=b[i]/a[i][i]-sum;
    //}
    for(i=0;i<(n*2)-h;i++){
        xold[i]=xnew[i];
    }
    t++;
    }
printf("\n");
 printf("Displacement(Ku=F,u=FK-¹) \n\n");
  for(i=0;i<(n*2)-h;i++){
        printf("u_DOF[%d]=%lf\n",pos[i],xold[i]);
    }
printf("\n");

return 0;
}
